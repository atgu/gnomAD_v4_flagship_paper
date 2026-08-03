"""LLM service for making calls to Anthropic and Google Vertex AI APIs."""
import os
import time
import anthropic
import traceback
from typing import Optional, Dict, Any
from config import ANTHROPIC_API_KEY, LLM_MAX_RETRIES, LLM_INITIAL_DELAY
from services.token_tracker import get_tracker

# No module-level Anthropic client. The one that used to live here was never
# read, but constructing it raised at import time when ANTHROPIC_API_KEY was
# unset, which made a Vertex-only run impossible. Each backend builds its own
# client on demand; see _call_anthropic_direct and _get_anthropic_vertex_client.

# Global flag to control verbose output
VERBOSE_MODE = False

# Vertex AI configuration (lazy initialization)
VERTEX_PROJECT = os.environ.get("GOOGLE_CLOUD_PROJECT", "guez-sandbox-aedc")
# us-central1 is the historical region for gemini-2.x. gemini-3.x models
# (3.5-flash, 3.1-flash-lite, ...) are only served from the special 'global'
# endpoint. We resolve the location per-model in _resolve_gemini_location.
VERTEX_LOCATION = os.environ.get("VERTEX_LOCATION", "us-central1")

# Claude 4.5+ on Vertex AI is recommended via the 'global' endpoint (best
# capacity, no regional premium). Older Claude variants only live in regional
# endpoints, so this can be overridden per-model later if needed.
ANTHROPIC_VERTEX_REGION = os.environ.get("ANTHROPIC_VERTEX_REGION", "global")

# Llama 4 MaaS models are only served from us-east5 on Vertex (as of 2026).
LLAMA_VERTEX_REGION = os.environ.get("LLAMA_VERTEX_REGION", "us-east5")

# Gemini per-run options. These can be tuned via set_gemini_options() at
# process startup. Workers in a ProcessPoolExecutor MUST call
# set_gemini_options() themselves (do not rely on inheriting global state
# through fork — it works on Linux but not on spawn-based platforms).
GEMINI_THINKING_BUDGET: int = 0  # 0 = no reasoning (cheapest, fastest)
GEMINI_FLEX: bool = False        # True = use Flex PayGo (-50% price, slower)
GEMINI_MAX_OUTPUT_TOKENS: int = 8192

# Lazily-built google.genai client (one per process). Built when first
# needed, after set_gemini_options() has had a chance to update the flags.
_genai_client = None


def set_gemini_options(
    thinking_budget: Optional[int] = None,
    flex: Optional[bool] = None,
    max_output_tokens: Optional[int] = None,
):
    """
    Configure Gemini call options for the current process.
    
    Must be called BEFORE the first Gemini API call. Subsequent calls
    rebuild the genai client only if the flex flag changes (because Flex
    is selected via HTTP headers baked into the client).
    """
    global GEMINI_THINKING_BUDGET, GEMINI_FLEX, GEMINI_MAX_OUTPUT_TOKENS, _genai_client
    rebuild = False
    if thinking_budget is not None:
        GEMINI_THINKING_BUDGET = int(thinking_budget)
    if flex is not None:
        if bool(flex) != GEMINI_FLEX:
            rebuild = True
        GEMINI_FLEX = bool(flex)
    if max_output_tokens is not None:
        GEMINI_MAX_OUTPUT_TOKENS = int(max_output_tokens)
    if rebuild:
        _genai_client = None


def _is_gemini_model(model: str) -> bool:
    """Check if the model is a Gemini model."""
    return model.lower().startswith("gemini")


def _is_anthropic_vertex_model(model: str) -> bool:
    """A Claude model addressed via Vertex AI rather than the Anthropic API.

    Convention: append '@vertex' to any Anthropic model id to route through
    Vertex. Plain 'claude-haiku-4-5' stays on the direct Anthropic API; the
    explicit 'claude-haiku-4-5@vertex' uses Vertex with GCP credits.
    """
    m = model.lower()
    return "@vertex" in m and "claude" in m


def _is_llama_vertex_model(model: str) -> bool:
    """A Llama model addressed via the Vertex MaaS OpenAI-compatible endpoint."""
    m = model.lower()
    return m.startswith("llama") or m.startswith("meta/llama")


def _normalize_anthropic_vertex_model(model: str) -> str:
    """Strip the '@vertex' suffix so the Vertex SDK gets the bare model id.

    Example: 'claude-haiku-4-5@vertex' -> 'claude-haiku-4-5'.
    """
    return model.replace("@vertex", "").strip()


# Map short aliases to the canonical Vertex MaaS model id required by the
# OpenAI-compatible endpoint. Keys must be checked AFTER an exact match in
# the canonical form.
_LLAMA_MAAS_ALIASES = {
    "llama-4-scout": "meta/llama-4-scout-17b-16e-instruct-maas",
    "llama-4-maverick": "meta/llama-4-maverick-17b-128e-instruct-maas",
    "llama-3.3-70b": "meta/llama-3.3-70b-instruct-maas",
    "llama-3-3-70b": "meta/llama-3.3-70b-instruct-maas",
}


def _normalize_llama_model(model: str) -> str:
    """Resolve a user-friendly Llama alias to the canonical Vertex MaaS id."""
    m = model.lower()
    if m.startswith("meta/"):
        # Already canonical (e.g. meta/llama-4-scout-17b-16e-instruct-maas).
        return model
    if m in _LLAMA_MAAS_ALIASES:
        return _LLAMA_MAAS_ALIASES[m]
    # Best-effort partial match (e.g. user passes 'llama-4-scout-17b-16e').
    for alias, canonical in _LLAMA_MAAS_ALIASES.items():
        if alias in m:
            return canonical
    # Fall back to whatever the caller passed; the Vertex API will 404 if
    # it's invalid, which is the clearest error message anyway.
    return model


def _resolve_gemini_location(model: str) -> str:
    """
    Resolve the Vertex AI location for a given Gemini model.
    
    Gemini 3.x models are served exclusively on the 'global' endpoint.
    Older Gemini 2.x / 1.x models stay on the legacy regional endpoint
    (default us-central1) for backward compatibility with run_005.
    """
    m = model.lower()
    if m.startswith("gemini-3"):
        return "global"
    return VERTEX_LOCATION


def _get_genai_client(location: str):
    """
    Build or return the cached google.genai client for the given location.
    
    The client is cached across calls within the same process. If the Flex
    flag is on, the client is configured with the Vertex shared/flex HTTP
    headers so that ALL requests through it are routed to Flex PayGo.
    """
    global _genai_client
    if _genai_client is not None and getattr(_genai_client, "_location", None) == location:
        return _genai_client
    
    from google import genai
    from google.genai.types import HttpOptions
    
    http_options = None
    if GEMINI_FLEX:
        http_options = HttpOptions(
            api_version="v1",
            headers={
                "X-Vertex-AI-LLM-Request-Type": "shared",
                "X-Vertex-AI-LLM-Shared-Request-Type": "flex",
            },
        )
    
    client = genai.Client(
        vertexai=True,
        project=VERTEX_PROJECT,
        location=location,
        http_options=http_options,
    )
    # Stash the location so we know when to rebuild
    client._location = location
    _genai_client = client
    return client


def _call_gemini(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """
    Call Gemini model via Vertex AI (google.genai SDK).
    
    Honours module-level GEMINI_THINKING_BUDGET, GEMINI_FLEX, and
    GEMINI_MAX_OUTPUT_TOKENS. The returned 'flex' boolean reflects whether
    Vertex actually used Flex PayGo (from traffic_type), not just our
    request preference — Vertex may fall back to Standard on its own.
    
    Returns:
        Dict with 'text', 'input_tokens', 'output_tokens', 'flex'
    """
    from google.genai.types import GenerateContentConfig, ThinkingConfig
    from google.api_core import exceptions as gcp_exceptions
    
    location = _resolve_gemini_location(model)
    client = _get_genai_client(location)
    
    config = GenerateContentConfig(
        temperature=temperature,
        max_output_tokens=GEMINI_MAX_OUTPUT_TOKENS,
        thinking_config=ThinkingConfig(thinking_budget=GEMINI_THINKING_BUDGET),
    )
    
    delay = initial_delay
    
    for attempt in range(max_retries):
        try:
            if VERBOSE_MODE:
                print(
                    f"INFO: Sending request to Gemini (model: {model}, temp: {temperature}, "
                    f"thinking_budget: {GEMINI_THINKING_BUDGET}, flex: {GEMINI_FLEX}, "
                    f"attempt: {attempt + 1}/{max_retries})..."
                )
            
            response = client.models.generate_content(
                model=model,
                contents=prompt,
                config=config,
            )
            
            input_tokens = 0
            output_tokens = 0
            actual_flex = False
            if response.usage_metadata is not None:
                input_tokens = response.usage_metadata.prompt_token_count or 0
                # candidates_token_count = visible response tokens.
                # thoughts_token_count  = hidden reasoning tokens, also billed
                # as output. We sum them so the cost report reflects what
                # Google actually charges.
                candidates = response.usage_metadata.candidates_token_count or 0
                thoughts = getattr(response.usage_metadata, "thoughts_token_count", 0) or 0
                output_tokens = candidates + thoughts
                # Vertex reports the real traffic class used for this request.
                traffic_type = getattr(response.usage_metadata, "traffic_type", None)
                # The field may be a str or an enum depending on SDK version.
                actual_flex = str(traffic_type).endswith("ON_DEMAND_FLEX")
            
            return {
                "text": response.text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": actual_flex,
            }
            
        except gcp_exceptions.ResourceExhausted as e:
            # Quota / 429 — common with Flex during peak load. Retry.
            print(f"WARNING: Gemini ResourceExhausted (likely quota / Flex saturation). Retrying in {delay}s... {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
                delay *= 2
            else:
                raise LLMOverloadedError(f"Gemini API exhausted after {max_retries} attempts: {e}")
        except gcp_exceptions.ServiceUnavailable as e:
            print(f"WARNING: Gemini ServiceUnavailable. Retrying in {delay}s... {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
                delay *= 2
            else:
                raise LLMOverloadedError(f"Gemini API unavailable after {max_retries} attempts: {e}")
        except Exception as e:
            error_str = str(e).lower()
            if any(err in error_str for err in ['429', '500', '503', 'overloaded', 'quota', 'resource exhausted', 'unavailable']):
                print(f"WARNING: Gemini API error. Retrying in {delay} seconds... Error: {e}")
                if attempt < max_retries - 1:
                    time.sleep(delay)
                    delay *= 2
                else:
                    raise LLMOverloadedError(f"Gemini API call failed after {max_retries} attempts. Last error: {e}")
            else:
                print(f"ERROR: Gemini API error: {e}")
                raise e
    
    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


# ---------------------------------------------------------------------------
# Anthropic Claude via Vertex AI (uses GCP credits instead of Anthropic API)
# ---------------------------------------------------------------------------

# One cached client per process keyed by region. AnthropicVertex handles ADC
# token refresh internally, so we don't need a manual refresh loop.
_anthropic_vertex_clients: Dict[str, Any] = {}


def _get_anthropic_vertex_client(region: str):
    """Build or return a cached AnthropicVertex client for the given region."""
    if region in _anthropic_vertex_clients:
        return _anthropic_vertex_clients[region]
    from anthropic import AnthropicVertex
    client_ = AnthropicVertex(region=region, project_id=VERTEX_PROJECT)
    _anthropic_vertex_clients[region] = client_
    return client_


def _call_anthropic_vertex(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """Call a Claude model via Vertex AI. Returns dict with text + token counts.

    Token accounting and response shape mirror the direct Anthropic SDK so the
    caller doesn't need to branch.
    """
    real_model = _normalize_anthropic_vertex_model(model)
    client_ = _get_anthropic_vertex_client(ANTHROPIC_VERTEX_REGION)
    delay = initial_delay

    # claude-opus-4.x and later models deprecate the temperature parameter.
    # Detect them by the presence of "-4-" followed by a version >= 6, or
    # simply by checking the known deprecation pattern.
    _no_temp_models = ("claude-opus-4-",)
    _omit_temperature = any(real_model.startswith(p) for p in _no_temp_models)

    for attempt in range(max_retries):
        try:
            if VERBOSE_MODE:
                print(f"INFO: Sending request to Anthropic Vertex (model: {real_model}, region: {ANTHROPIC_VERTEX_REGION}, attempt: {attempt + 1}/{max_retries})...")
            create_kwargs = dict(
                model=real_model,
                max_tokens=4096,
                messages=[{"role": "user", "content": prompt}],
            )
            if not _omit_temperature:
                create_kwargs["temperature"] = temperature
            message = client_.messages.create(**create_kwargs)
            input_tokens = message.usage.input_tokens if hasattr(message, "usage") else 0
            output_tokens = message.usage.output_tokens if hasattr(message, "usage") else 0
            return {
                "text": message.content[0].text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": False,
            }
        except anthropic.APIStatusError as e:
            # Vertex shapes its errors via the Anthropic SDK so we get the
            # same status_code semantics as the direct API.
            if e.status_code in [429, 500, 503, 529]:
                if VERBOSE_MODE:
                    print(f"WARNING: Anthropic Vertex returned {e.status_code}. Retrying in {delay}s...")
                if attempt < max_retries - 1:
                    time.sleep(delay)
                    delay *= 2
                else:
                    raise LLMOverloadedError(f"Anthropic Vertex call failed after {max_retries} attempts. Last error: {e}")
            else:
                print(f"ERROR: Anthropic Vertex API error: {e}")
                raise e
        except Exception as e:
            print(f"ERROR: Unexpected error during Anthropic Vertex call: {e}")
            raise e

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


# ---------------------------------------------------------------------------
# Meta Llama via Vertex AI MaaS (OpenAI-compatible endpoint)
# ---------------------------------------------------------------------------

# GCP access tokens issued by ADC live ~1h. We refresh proactively 5 minutes
# before expiry to avoid mid-call 401s in long-running workers.
_GCP_TOKEN_TTL_SECONDS = 3600
_GCP_TOKEN_REFRESH_MARGIN = 300

_llama_vertex_client = None
_llama_vertex_token_expires_at: float = 0.0


def _get_llama_vertex_client(force_refresh: bool = False):
    """Return an OpenAI client pointed at the Vertex MaaS endpoint.

    The OpenAI client expects a static api_key, but GCP access tokens expire,
    so we rebuild the client whenever the cached token is near expiry. This
    is safe to call on every request.
    """
    global _llama_vertex_client, _llama_vertex_token_expires_at

    now = time.time()
    needs_refresh = (
        force_refresh
        or _llama_vertex_client is None
        or now >= (_llama_vertex_token_expires_at - _GCP_TOKEN_REFRESH_MARGIN)
    )
    if not needs_refresh:
        return _llama_vertex_client

    from google.auth import default as gauth_default
    from google.auth.transport.requests import Request as GAuthRequest
    from openai import OpenAI

    credentials, _ = gauth_default(
        scopes=["https://www.googleapis.com/auth/cloud-platform"]
    )
    credentials.refresh(GAuthRequest())
    _llama_vertex_token_expires_at = now + _GCP_TOKEN_TTL_SECONDS

    endpoint_url = (
        f"https://{LLAMA_VERTEX_REGION}-aiplatform.googleapis.com/v1beta1/"
        f"projects/{VERTEX_PROJECT}/locations/{LLAMA_VERTEX_REGION}/endpoints/openapi"
    )
    _llama_vertex_client = OpenAI(base_url=endpoint_url, api_key=credentials.token)
    return _llama_vertex_client


def _call_llama_vertex(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """Call a Llama MaaS model via Vertex. Returns dict with text + token counts.

    Llama MaaS speaks the OpenAI Chat Completions protocol, so the response
    shape is OpenAI-style (`choices[0].message.content`, `usage.prompt_tokens`,
    `usage.completion_tokens`). We map it to our internal dict.
    """
    real_model = _normalize_llama_model(model)
    delay = initial_delay

    for attempt in range(max_retries):
        try:
            client_ = _get_llama_vertex_client()
            if VERBOSE_MODE:
                print(f"INFO: Sending request to Llama Vertex (model: {real_model}, region: {LLAMA_VERTEX_REGION}, attempt: {attempt + 1}/{max_retries})...")
            response = client_.chat.completions.create(
                model=real_model,
                max_tokens=4096,
                temperature=temperature,
                messages=[{"role": "user", "content": prompt}],
            )
            usage = getattr(response, "usage", None)
            input_tokens = getattr(usage, "prompt_tokens", 0) if usage else 0
            output_tokens = getattr(usage, "completion_tokens", 0) if usage else 0
            text = response.choices[0].message.content or ""
            return {
                "text": text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": False,
            }
        except Exception as e:
            err_str = str(e)
            err_lower = err_str.lower()
            # 401 = expired/invalid token. Force a refresh before the next try.
            if "401" in err_str or "unauthorized" in err_lower:
                if VERBOSE_MODE:
                    print("WARNING: Llama Vertex 401, forcing token refresh.")
                try:
                    _get_llama_vertex_client(force_refresh=True)
                except Exception as refresh_err:
                    print(f"ERROR: Failed to refresh GCP token: {refresh_err}")
            # Retry on transient errors. We classify by HTTP status string
            # because the openai SDK wraps a variety of exception types and
            # we don't want to be brittle about which one we caught.
            is_retryable = any(
                code in err_str for code in ("429", "500", "502", "503", "504", "529")
            ) or "timeout" in err_lower
            if attempt < max_retries - 1 and (is_retryable or "401" in err_str):
                if VERBOSE_MODE:
                    print(f"WARNING: Llama Vertex error: {e}. Retrying in {delay}s...")
                time.sleep(delay)
                delay *= 2
            else:
                print(f"ERROR: Llama Vertex API error: {e}")
                if attempt == max_retries - 1 and is_retryable:
                    raise LLMOverloadedError(f"Llama Vertex call failed after {max_retries} attempts. Last error: {e}")
                raise

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


# ---------------------------------------------------------------------------
# Google MedGemma via a Vertex AI dedicated endpoint (vLLM, OpenAI-compatible)
# ---------------------------------------------------------------------------
# Unlike Llama MaaS (a shared OpenAI route), MedGemma is self-hosted on a
# dedicated endpoint we deploy ourselves. Its URL and served-model id are not
# known until deploy time, so we read them from data/medgemma_endpoint.json
# (written by benchmark/scripts/deploy_medgemma.py).

_MEDGEMMA_INFO_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "medgemma_endpoint.json",
)

_medgemma_client = None
_medgemma_token_expires_at: float = 0.0
_medgemma_info: Optional[Dict[str, Any]] = None


def _is_medgemma_model(model: str) -> bool:
    """A MedGemma model addressed via our Vertex dedicated endpoint."""
    return model.lower().startswith("medgemma")


def _load_medgemma_info() -> Dict[str, Any]:
    """Load (and cache) the deployed-endpoint metadata JSON."""
    global _medgemma_info
    if _medgemma_info is not None:
        return _medgemma_info
    import json
    if not os.path.exists(_MEDGEMMA_INFO_PATH):
        raise RuntimeError(
            f"MedGemma endpoint info not found at {_MEDGEMMA_INFO_PATH}. "
            "Deploy it first: python3 benchmark/scripts/deploy_medgemma.py"
        )
    with open(_MEDGEMMA_INFO_PATH) as f:
        _medgemma_info = json.load(f)
    return _medgemma_info


def _get_medgemma_client(force_refresh: bool = False):
    """Return an OpenAI client pointed at the MedGemma dedicated endpoint.

    Same expiring-token handling as the Llama client: GCP ADC tokens last ~1h,
    so we rebuild the client when the cached token nears expiry.
    """
    global _medgemma_client, _medgemma_token_expires_at

    now = time.time()
    needs_refresh = (
        force_refresh
        or _medgemma_client is None
        or now >= (_medgemma_token_expires_at - _GCP_TOKEN_REFRESH_MARGIN)
    )
    if not needs_refresh:
        return _medgemma_client

    from google.auth import default as gauth_default
    from google.auth.transport.requests import Request as GAuthRequest
    from openai import OpenAI

    info = _load_medgemma_info()
    credentials, _ = gauth_default(
        scopes=["https://www.googleapis.com/auth/cloud-platform"]
    )
    credentials.refresh(GAuthRequest())
    _medgemma_token_expires_at = now + _GCP_TOKEN_TTL_SECONDS

    dns = info["dedicated_endpoint_dns"]
    base_url = (
        f"https://{dns}/v1beta1/"
        f"projects/{info['project']}/locations/{info['location']}/"
        f"endpoints/{info['endpoint_id']}"
    )
    _medgemma_client = OpenAI(base_url=base_url, api_key=credentials.token)
    return _medgemma_client


def _call_medgemma_vertex(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """Call MedGemma on its dedicated vLLM endpoint (OpenAI Chat Completions)."""
    info = _load_medgemma_info()
    served_model = info.get("served_model_name", "medgemma-27b-text-it")
    delay = initial_delay

    for attempt in range(max_retries):
        try:
            client_ = _get_medgemma_client()
            if VERBOSE_MODE:
                print(f"INFO: Sending request to MedGemma (model: {served_model}, attempt: {attempt + 1}/{max_retries})...")
            response = client_.chat.completions.create(
                model=served_model,
                max_tokens=4096,
                temperature=temperature,
                messages=[{"role": "user", "content": prompt}],
                # vLLM-native repetition penalty. Greedy decoding (temp=0) without
                # an anti-repeat term makes MedGemma fall into a whitespace-padding
                # loop on the wide Markdown tables (esp. the A1 association table),
                # emitting up to 4096 blank tokens and zero parseable rows. 1.3
                # breaks the loop and stays near-greedy. NB: stronger anti-repeat
                # (frequency_penalty) instead pushes it into multilingual gibberish,
                # so keep repetition_penalty alone and run at low concurrency (the
                # degeneration is sensitive to vLLM continuous-batching nondeterminism).
                extra_body={"repetition_penalty": 1.3},
            )
            usage = getattr(response, "usage", None)
            input_tokens = getattr(usage, "prompt_tokens", 0) if usage else 0
            output_tokens = getattr(usage, "completion_tokens", 0) if usage else 0
            text = response.choices[0].message.content or ""
            return {
                "text": text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": False,
            }
        except Exception as e:
            err_str = str(e)
            err_lower = err_str.lower()
            if "401" in err_str or "unauthorized" in err_lower:
                if VERBOSE_MODE:
                    print("WARNING: MedGemma 401, forcing token refresh.")
                try:
                    _get_medgemma_client(force_refresh=True)
                except Exception as refresh_err:
                    print(f"ERROR: Failed to refresh GCP token: {refresh_err}")
            is_retryable = any(
                code in err_str for code in ("429", "500", "502", "503", "504", "529")
            ) or "timeout" in err_lower
            if attempt < max_retries - 1 and (is_retryable or "401" in err_str):
                if VERBOSE_MODE:
                    print(f"WARNING: MedGemma error: {e}. Retrying in {delay}s...")
                time.sleep(delay)
                delay *= 2
            else:
                print(f"ERROR: MedGemma endpoint error: {e}")
                if attempt == max_retries - 1 and is_retryable:
                    raise LLMOverloadedError(f"MedGemma call failed after {max_retries} attempts. Last error: {e}")
                raise

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


# ---------------------------------------------------------------------------
# Qwen3.6-27B via a Vertex AI dedicated endpoint (vLLM, OpenAI-compatible)
# ---------------------------------------------------------------------------
# Self-hosted on a dedicated endpoint we deploy ourselves (same pattern as
# MedGemma). URL + served-model id are written at deploy time to
# data/qwen_endpoint.json (by benchmark/scripts/deploy_qwen.py).
# Qwen3.6 is a reasoning model; we disable thinking at request time via the
# chat-template kwarg so outputs are short, parseable and cheap.

_QWEN_INFO_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "qwen_endpoint.json",
)

_qwen_client = None
_qwen_token_expires_at: float = 0.0
_qwen_info: Optional[Dict[str, Any]] = None


def _is_qwen_model(model: str) -> bool:
    """A Qwen model addressed via our Vertex dedicated endpoint."""
    return model.lower().startswith("qwen")


def _load_qwen_info() -> Dict[str, Any]:
    """Load (and cache) the deployed-endpoint metadata JSON."""
    global _qwen_info
    if _qwen_info is not None:
        return _qwen_info
    import json
    if not os.path.exists(_QWEN_INFO_PATH):
        raise RuntimeError(
            f"Qwen endpoint info not found at {_QWEN_INFO_PATH}. "
            "Deploy it first: python3 benchmark/scripts/deploy_qwen.py"
        )
    with open(_QWEN_INFO_PATH) as f:
        _qwen_info = json.load(f)
    return _qwen_info


def _get_qwen_client(force_refresh: bool = False):
    """Return an OpenAI client pointed at the Qwen dedicated endpoint.

    Same expiring-token handling as MedGemma: GCP ADC tokens last ~1h, so we
    rebuild the client when the cached token nears expiry.
    """
    global _qwen_client, _qwen_token_expires_at

    now = time.time()
    needs_refresh = (
        force_refresh
        or _qwen_client is None
        or now >= (_qwen_token_expires_at - _GCP_TOKEN_REFRESH_MARGIN)
    )
    if not needs_refresh:
        return _qwen_client

    from google.auth import default as gauth_default
    from google.auth.transport.requests import Request as GAuthRequest
    from openai import OpenAI

    info = _load_qwen_info()
    credentials, _ = gauth_default(
        scopes=["https://www.googleapis.com/auth/cloud-platform"]
    )
    credentials.refresh(GAuthRequest())
    _qwen_token_expires_at = now + _GCP_TOKEN_TTL_SECONDS

    dns = info["dedicated_endpoint_dns"]
    base_url = (
        f"https://{dns}/v1beta1/"
        f"projects/{info['project']}/locations/{info['location']}/"
        f"endpoints/{info['endpoint_id']}"
    )
    _qwen_client = OpenAI(base_url=base_url, api_key=credentials.token)
    return _qwen_client


def _call_qwen_vertex(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """Call Qwen3.6-27B on its dedicated vLLM endpoint (OpenAI Chat Completions).

    Thinking is disabled via chat_template_kwargs so the model returns only the
    final answer (no <think> trace): keeps outputs short, parseable and cheap.
    """
    info = _load_qwen_info()
    served_model = info.get("served_model_name", "qwen3.6-27b")
    delay = initial_delay

    for attempt in range(max_retries):
        try:
            client_ = _get_qwen_client()
            if VERBOSE_MODE:
                print(f"INFO: Sending request to Qwen (model: {served_model}, attempt: {attempt + 1}/{max_retries})...")
            response = client_.chat.completions.create(
                model=served_model,
                max_tokens=4096,
                temperature=temperature,
                messages=[{"role": "user", "content": prompt}],
                # Disable Qwen's reasoning/thinking mode so the response is just
                # the final answer (no <think> trace). Mirrors how we use the
                # other agents (non-reasoning) and avoids token blow-up.
                extra_body={"chat_template_kwargs": {"enable_thinking": False}},
            )
            usage = getattr(response, "usage", None)
            input_tokens = getattr(usage, "prompt_tokens", 0) if usage else 0
            output_tokens = getattr(usage, "completion_tokens", 0) if usage else 0
            text = response.choices[0].message.content or ""
            return {
                "text": text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": False,
            }
        except Exception as e:
            err_str = str(e)
            err_lower = err_str.lower()
            if "401" in err_str or "unauthorized" in err_lower:
                if VERBOSE_MODE:
                    print("WARNING: Qwen 401, forcing token refresh.")
                try:
                    _get_qwen_client(force_refresh=True)
                except Exception as refresh_err:
                    print(f"ERROR: Failed to refresh GCP token: {refresh_err}")
            is_retryable = any(
                code in err_str for code in ("429", "500", "502", "503", "504", "529")
            ) or "timeout" in err_lower
            if attempt < max_retries - 1 and (is_retryable or "401" in err_str):
                if VERBOSE_MODE:
                    print(f"WARNING: Qwen error: {e}. Retrying in {delay}s...")
                time.sleep(delay)
                delay *= 2
            else:
                print(f"ERROR: Qwen endpoint error: {e}")
                if attempt == max_retries - 1 and is_retryable:
                    raise LLMOverloadedError(f"Qwen call failed after {max_retries} attempts. Last error: {e}")
                raise

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


# ---------------------------------------------------------------------------
# DeepSeek via Vertex AI Model-as-a-Service (managed, OpenAI-compatible)
# ---------------------------------------------------------------------------
# Fully managed/serverless MaaS (like Llama): no GPU to deploy, pay-per-token.
# Served on the 'global' location. DeepSeek-V3.2 returns its answer directly in
# `content` (no separate <think>/reasoning channel on Vertex MaaS); we still
# pass thinking=False to keep it in non-reasoning mode and avoid token bloat.

_DEEPSEEK_MAAS_LOCATION = "global"

_deepseek_client = None
_deepseek_token_expires_at: float = 0.0

_DEEPSEEK_ALIASES = {
    "deepseek-v3.2": "deepseek-ai/deepseek-v3.2-maas",
    "deepseek-v3.2-maas": "deepseek-ai/deepseek-v3.2-maas",
    "deepseek-v3.1": "deepseek-ai/deepseek-v3.1-maas",
    "deepseek-v3.1-maas": "deepseek-ai/deepseek-v3.1-maas",
    "deepseek-r1": "deepseek-ai/deepseek-r1-0528-maas",
    "deepseek": "deepseek-ai/deepseek-v3.2-maas",
}


def _is_deepseek_model(model: str) -> bool:
    """A DeepSeek model addressed via the Vertex MaaS OpenAI-compatible endpoint."""
    return model.lower().startswith("deepseek")


def _normalize_deepseek_model(model: str) -> str:
    """Resolve a friendly DeepSeek alias to the canonical Vertex MaaS id."""
    m = model.lower()
    if m.startswith("deepseek-ai/"):
        return model
    if m in _DEEPSEEK_ALIASES:
        return _DEEPSEEK_ALIASES[m]
    for alias, canonical in _DEEPSEEK_ALIASES.items():
        if alias in m:
            return canonical
    return model


def _get_deepseek_maas_client(force_refresh: bool = False):
    """Return an OpenAI client pointed at the DeepSeek MaaS global endpoint."""
    global _deepseek_client, _deepseek_token_expires_at

    now = time.time()
    needs_refresh = (
        force_refresh
        or _deepseek_client is None
        or now >= (_deepseek_token_expires_at - _GCP_TOKEN_REFRESH_MARGIN)
    )
    if not needs_refresh:
        return _deepseek_client

    from google.auth import default as gauth_default
    from google.auth.transport.requests import Request as GAuthRequest
    from openai import OpenAI

    credentials, _ = gauth_default(
        scopes=["https://www.googleapis.com/auth/cloud-platform"]
    )
    credentials.refresh(GAuthRequest())
    _deepseek_token_expires_at = now + _GCP_TOKEN_TTL_SECONDS

    endpoint_url = (
        f"https://aiplatform.googleapis.com/v1/"
        f"projects/{VERTEX_PROJECT}/locations/{_DEEPSEEK_MAAS_LOCATION}/endpoints/openapi"
    )
    _deepseek_client = OpenAI(base_url=endpoint_url, api_key=credentials.token)
    return _deepseek_client


# DeepSeek MaaS has a low per-project concurrent-request quota. When exceeded
# it returns HTTP 200 with choices=None and error="The request is throttled due
# to too many concurrent requests." (NOT a 429). We must detect this soft
# throttle and back off generously, independently of the small global retry
# budget used for hard HTTP errors.
_DEEPSEEK_THROTTLE_MAX_ATTEMPTS = 10
_DEEPSEEK_THROTTLE_BASE_DELAY = 2.0
_DEEPSEEK_THROTTLE_MAX_DELAY = 30.0


def _deepseek_throttle_reason(response) -> Optional[str]:
    """Return the throttle/error message if the response carries no usable
    choices (soft throttle or server-side error wrapped in a 200), else None."""
    choices = getattr(response, "choices", None)
    if choices:
        return None
    # The error string lives in the raw payload, not a typed attribute.
    err = getattr(response, "error", None)
    if err is None:
        try:
            err = response.model_dump().get("error")
        except Exception:
            err = None
    return str(err) if err else "empty response (no choices)"


def _call_deepseek_maas(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """Call a DeepSeek MaaS model via Vertex (OpenAI Chat Completions protocol).

    Handles DeepSeek MaaS's soft concurrency throttle (200 + choices=None) with
    a dedicated, jittered exponential backoff so individual calls survive bursts
    rather than crashing on `choices[0]`.
    """
    import random

    real_model = _normalize_deepseek_model(model)
    delay = initial_delay
    throttle_delay = _DEEPSEEK_THROTTLE_BASE_DELAY
    throttle_attempts = 0

    attempt = 0
    while attempt < max_retries:
        try:
            client_ = _get_deepseek_maas_client()
            if VERBOSE_MODE:
                print(f"INFO: Sending request to DeepSeek MaaS (model: {real_model}, attempt: {attempt + 1}/{max_retries})...")
            response = client_.chat.completions.create(
                model=real_model,
                max_tokens=4096,
                temperature=temperature,
                messages=[{"role": "user", "content": prompt}],
                # Keep DeepSeek's hybrid thinking off: answer goes straight to
                # `content`, no reasoning trace, no token bloat.
                extra_body={"chat_template_kwargs": {"thinking": False}},
            )

            # Soft throttle / empty response: retry generously without consuming
            # the (small) hard-error retry budget.
            throttle_reason = _deepseek_throttle_reason(response)
            if throttle_reason is not None:
                throttle_attempts += 1
                if throttle_attempts >= _DEEPSEEK_THROTTLE_MAX_ATTEMPTS:
                    raise LLMOverloadedError(
                        f"DeepSeek MaaS throttled {throttle_attempts}x: {throttle_reason}"
                    )
                sleep_s = min(throttle_delay, _DEEPSEEK_THROTTLE_MAX_DELAY) * (0.5 + random.random())
                if VERBOSE_MODE:
                    print(f"WARNING: DeepSeek MaaS throttled ({throttle_reason}). Backing off {sleep_s:.1f}s (throttle attempt {throttle_attempts}).")
                time.sleep(sleep_s)
                throttle_delay = min(throttle_delay * 2, _DEEPSEEK_THROTTLE_MAX_DELAY)
                continue  # do NOT increment `attempt`; throttle has its own budget

            usage = getattr(response, "usage", None)
            input_tokens = getattr(usage, "prompt_tokens", 0) if usage else 0
            output_tokens = getattr(usage, "completion_tokens", 0) if usage else 0
            text = response.choices[0].message.content or ""
            return {
                "text": text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": False,
            }
        except LLMOverloadedError:
            raise
        except Exception as e:
            err_str = str(e)
            err_lower = err_str.lower()
            if "401" in err_str or "unauthorized" in err_lower:
                if VERBOSE_MODE:
                    print("WARNING: DeepSeek MaaS 401, forcing token refresh.")
                try:
                    _get_deepseek_maas_client(force_refresh=True)
                except Exception as refresh_err:
                    print(f"ERROR: Failed to refresh GCP token: {refresh_err}")
            is_retryable = any(
                code in err_str for code in ("429", "500", "502", "503", "504", "529")
            ) or "timeout" in err_lower
            if attempt < max_retries - 1 and (is_retryable or "401" in err_str):
                if VERBOSE_MODE:
                    print(f"WARNING: DeepSeek MaaS error: {e}. Retrying in {delay}s...")
                time.sleep(delay)
                delay *= 2
                attempt += 1
            else:
                print(f"ERROR: DeepSeek MaaS API error: {e}")
                if attempt == max_retries - 1 and is_retryable:
                    raise LLMOverloadedError(f"DeepSeek MaaS call failed after {max_retries} attempts. Last error: {e}")
                raise
        else:
            attempt += 1

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


# Define a retry decorator
def retry_with_backoff(retries=LLM_MAX_RETRIES, initial_delay=LLM_INITIAL_DELAY):
    # ... decorator implementation ...
    pass

class LLMOverloadedError(Exception):
    """Exception raised when LLM API is overloaded after retries."""
    pass


# ---------------------------------------------------------------------------
# Anthropic Claude via the direct Anthropic API (default backend)
# ---------------------------------------------------------------------------

def _call_anthropic_direct(
    prompt: str,
    model: str,
    temperature: float,
    max_retries: int,
    initial_delay: float,
) -> Dict[str, Any]:
    """Call a Claude model via the official Anthropic API (uses ANTHROPIC_API_KEY).

    This is the legacy path. Kept as a dedicated function so the routing in
    `_dispatch_llm` is symmetric across all backends.
    """
    client_ = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    delay = initial_delay

    for attempt in range(max_retries):
        try:
            if VERBOSE_MODE:
                print(f"INFO: Sending request to Anthropic API (model: {model}, temp: {temperature}, attempt: {attempt + 1}/{max_retries})...")
            message = client_.messages.create(
                model=model,
                max_tokens=4096,
                temperature=temperature,
                messages=[{"role": "user", "content": prompt}],
            )
            input_tokens = message.usage.input_tokens if hasattr(message, "usage") else 0
            output_tokens = message.usage.output_tokens if hasattr(message, "usage") else 0
            return {
                "text": message.content[0].text,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "flex": False,
            }
        except anthropic.APIStatusError as e:
            if e.status_code in [529, 500, 503]:
                print(f"WARNING: Anthropic API returned status {e.status_code}. Retrying in {delay} seconds...")
                if attempt < max_retries - 1:
                    time.sleep(delay)
                    delay *= 2
                else:
                    raise LLMOverloadedError(f"Anthropic API call failed after {max_retries} attempts. Last error: {e}")
            else:
                raise e
        except Exception as e:
            print(f"ERROR: An unexpected error occurred during Anthropic API call: {e}")
            raise e

    return {"text": None, "input_tokens": 0, "output_tokens": 0, "flex": False}


def _dispatch_llm(prompt: str, model: str, temperature: float) -> Dict[str, Any]:
    """Route an LLM call to the right backend based on the model id.

    Backends (in resolution order):
        1. Gemini (model id starts with 'gemini') -> Vertex AI
        2. Claude with '@vertex' suffix          -> Anthropic on Vertex AI
        3. Llama / 'meta/llama-*'                -> Vertex AI MaaS (OpenAI proto)
        4. Anything else                         -> Anthropic direct API

    Returns dict {text, input_tokens, output_tokens, flex} so call sites
    can do token tracking uniformly.
    """
    if _is_gemini_model(model):
        return _call_gemini(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)
    if _is_medgemma_model(model):
        return _call_medgemma_vertex(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)
    if _is_qwen_model(model):
        return _call_qwen_vertex(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)
    if _is_anthropic_vertex_model(model):
        return _call_anthropic_vertex(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)
    if _is_llama_vertex_model(model):
        return _call_llama_vertex(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)
    if _is_deepseek_model(model):
        return _call_deepseek_maas(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)
    return _call_anthropic_direct(prompt, model, temperature, LLM_MAX_RETRIES, LLM_INITIAL_DELAY)


def call_llm(
    prompt: str,
    model: str = "claude-haiku-4-5",
    temperature: float = 0.0,
    agent_name: Optional[str] = None,
) -> str:
    """
    Calls the specified LLM with a given prompt.
    Includes retry logic with exponential backoff for overloaded cases.
    Supports both Anthropic (Claude) and Google (Gemini) models.
    
    Args:
        prompt: The prompt to send to the LLM
        model: Model to use (default: claude-haiku-4-5). Use "gemini-2.5-flash" for Gemini.
        temperature: Temperature setting (default: 0.0)
        agent_name: Optional name of the calling agent for token tracking
        
    Returns:
        The text response from the LLM
    """
    result = _dispatch_llm(prompt, model, temperature)

    if agent_name and result.get("input_tokens", 0) > 0:
        tracker = get_tracker()
        # We strip the '@vertex' suffix for accounting so a Claude model is
        # recorded under its canonical name regardless of the backend used.
        tracking_model = _normalize_anthropic_vertex_model(model) if _is_anthropic_vertex_model(model) else model
        tracker.record(
            agent_name,
            tracking_model,
            result["input_tokens"],
            result["output_tokens"],
            flex=result.get("flex", False),
        )

    return result["text"]


def call_llm_with_usage(
    prompt: str,
    model: str = "claude-haiku-4-5",
    temperature: float = 0.0,
    agent_name: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Calls the specified LLM and returns both text and token usage.
    Supports both Anthropic (Claude) and Google (Gemini) models.
    
    Args:
        prompt: The prompt to send to the LLM
        model: Model to use (default: claude-haiku-4-5). Use "gemini-2.5-flash" for Gemini.
        temperature: Temperature setting (default: 0.0)
        agent_name: Optional name of the calling agent for token tracking
        
    Returns:
        Dict with 'text', 'input_tokens', 'output_tokens'
    """
    result = _dispatch_llm(prompt, model, temperature)

    if agent_name and result.get("input_tokens", 0) > 0:
        tracker = get_tracker()
        tracking_model = _normalize_anthropic_vertex_model(model) if _is_anthropic_vertex_model(model) else model
        tracker.record(
            agent_name,
            tracking_model,
            result["input_tokens"],
            result["output_tokens"],
            flex=result.get("flex", False),
        )

    # call_llm_with_usage historically returned 3 keys; we keep the contract
    # unchanged but pass through 'flex' for callers that want it.
    return {
        "text": result.get("text"),
        "input_tokens": result.get("input_tokens", 0),
        "output_tokens": result.get("output_tokens", 0),
        "flex": result.get("flex", False),
    }

