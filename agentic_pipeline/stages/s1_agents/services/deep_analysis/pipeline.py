"""Deep analysis pipeline orchestrating secondary Mendelian agents."""
from __future__ import annotations

import os
import random
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from time import perf_counter
from typing import Dict, List, Optional, Callable, Tuple

from parsers.llm_parsers import (
    parse_deep_analysis_a1,
    parse_deep_analysis_a2,
    parse_deep_analysis_a3,
    parse_deep_analysis_a4,
    parse_onset_severity_response,
    parse_pathogenicity_response,
    parse_algorithmic_summary_response,
)
from services.llm_service import call_llm, call_llm_with_usage, LLMOverloadedError
from services.token_tracker import get_tracker
from utils.helpers import load_prompt

from .constants import (
    ASSOCIATION_LABELS,
    ASSOCIATION_PENALTIES,
    DEFAULT_ALGO_LEVEL,
    INHERITANCE_PENALTIES,
    INHERITANCE_SCORE_LABELS,
    MECHANISM_KEYWORDS,
    ONSET_PENALTIES,
    ONSET_SCORE_LABELS,
    PENETRANCE_KEYWORDS,
    INHERITANCE_KEYWORDS,
    MAX_DISEASES_FOR_QUERY,
    PENETRANCE_LEVEL_LABELS,
    SEVERITY_LEVEL_LABELS,
)
from .helpers import (
    articles_to_text,
    combine_articles,
    fetch_articles_for_query,
    format_association_summaries,
    format_disease_list,
    prepare_disease_names,
    build_keyword_clause,
    build_or_clause,
)


# ---------------------------------------------------------------------------
# Helper for prompt variants (knowledge, simple, proba)
# ---------------------------------------------------------------------------

# Prompts that have knowledge/simple variants
PROMPTS_WITH_ARTICLE_VARIANTS = (
    "deep_analysis_disease_agent",
    "deep_analysis_penetrance_agent",
    "deep_analysis_inheritance_agent",
    "deep_analysis_mechanism_agent",
)

# All prompts that have simple variants
PROMPTS_WITH_SIMPLE_VARIANTS = (
    "deep_analysis_disease_agent",
    "deep_analysis_penetrance_agent",
    "deep_analysis_inheritance_agent",
    "deep_analysis_mechanism_agent",
    "deep_analysis_onset_severity_agent",
    "deep_analysis_pathogenicity_agent",
    "deep_analysis_algorithmic_summary_agent",
)

# Prompts that have proba variants (probabilistic output)
PROMPTS_WITH_PROBA_VARIANTS = (
    "deep_analysis_penetrance_agent",
    "deep_analysis_inheritance_agent",
    "deep_analysis_onset_severity_agent",
)


def _get_prompt_name(base_name: str, knowledge_mode: bool = False, simple_mode: bool = False, proba_mode: bool = False) -> str:
    """Return the prompt path with subdirectory based on mode.
    
    Returns paths like 'simple/deep_analysis_disease_agent' instead of 
    'deep_analysis_disease_agent_simple'.
    
    Priority order: knowledge > proba > simple > normal

    Knowledge mode takes highest precedence so that, when it is enabled, the
    article-free knowledge/ variants are loaded even if proba is also on. The
    knowledge/ penetrance & inheritance prompts are kept in the same
    probabilistic output format as their proba/ counterparts, so the only
    difference between knowledge+proba and proba alone is abstracts vs.
    parametric knowledge. The response parser is selected independently via the
    proba_mode flag, so it stays compatible with these knowledge/ prompts.
    """
    # Mechanism agent: v2 by default in this repository.
    #
    # run_016 was originally scored with the v1 prompt, which forced a single
    # mechanism per disease and labelled multi-mechanism cases "Conflicting".
    # In March 2026 the 562 affected genes were re-annotated with v2, which
    # allows composite mechanisms in slash notation (e.g. "DN/LoF"), producing
    # 480 of them; the published tables reflect that state. A fresh run on the
    # v1 prompt would therefore produce no composite mechanism at all, and
    # --composite-mode strict would exclude nothing, silently changing DisPo.
    # Set PEPPER_MECHANISM_PROMPT_VERSION=v1 to reproduce the original scoring.
    if base_name == "deep_analysis_mechanism_agent":
        version = os.environ.get("PEPPER_MECHANISM_PROMPT_VERSION", "v2").lower()
        if version == "v2":
            return "normal/deep_analysis_mechanism_agent_v2"

    # Knowledge mode takes highest precedence for its article-free variants
    if knowledge_mode and base_name in PROMPTS_WITH_ARTICLE_VARIANTS:
        return f"knowledge/{base_name}"
    # Proba mode next for supported prompts
    if proba_mode and base_name in PROMPTS_WITH_PROBA_VARIANTS:
        return f"proba/{base_name}"
    # Simple mode
    if simple_mode and base_name in PROMPTS_WITH_SIMPLE_VARIANTS:
        return f"simple/{base_name}"
    return f"normal/{base_name}"


# ---------------------------------------------------------------------------
# A1 Validation and Retry Logic
# ---------------------------------------------------------------------------

# Patterns that indicate parsing errors (not valid disease names)
INVALID_DISEASE_PATTERNS = [
    (r'^PMID\s*:', "looks like a PMID reference, not a disease name"),
    (r'^\d{5,}$', "looks like a numeric ID, not a disease name"),
    (r'^article$', "looks like a table header"),
    (r'^evidence(\s+type)?$', "looks like a table header"),
    (r'^disease(\s+name)?$', "looks like a table header"),
    (r'^human\s+disease', "looks like a table header"),
]


def validate_a1_response(parsed_results: List[Dict]) -> Tuple[bool, Optional[str]]:
    """
    Validate parsed A1 results for format issues.
    Returns (is_valid, error_message).
    """
    if not parsed_results:
        return False, "No diseases were parsed from your response. Please output a valid Markdown table."
    
    # Check for PMIDs or headers mistakenly parsed as disease names
    for r in parsed_results:
        disease = r.get("disease", "").strip()
        disease_lower = disease.lower()
        
        for pattern, error_desc in INVALID_DISEASE_PATTERNS:
            if re.match(pattern, disease_lower, re.I):
                return False, f"'{disease}' {error_desc}. Please use actual disease names in the Disease column."
    
    # If we have multiple diseases but ALL have None scores, something is wrong
    # (Exception: single "None" entry for no-phenotype case is valid)
    if len(parsed_results) > 1:
        all_scores_none = all(r.get("association_score") is None for r in parsed_results)
        if all_scores_none:
            return False, "All association scores are missing. Please include numeric scores (1-4) in the AssociationScore column."
    
    return True, None


def call_a1_with_retry(
    base_prompt: str,
    model: str,
    temperature: float,
    gene_name: str,
    max_retries: int = 3,
) -> Tuple[List[Dict], int, Dict]:
    """
    Call A1 agent with retry logic on format errors.
    Returns (parsed_results, attempt_count).
    """
    prompt = base_prompt
    raw_attempts = []
    
    for attempt in range(1, max_retries + 1):
        # Add retry feedback after first attempt
        if attempt > 1:
            prompt = base_prompt + f"""

⚠️ IMPORTANT FORMAT CORRECTION (Attempt {attempt}/{max_retries}):
Your previous response had format issues: {error_msg}

Please respond with a proper Markdown table with these exact columns:
| Disease | AssociationScore | Justification | PMIDs |

- Disease: The actual disease name (not PMIDs or headers)
- AssociationScore: Integer 1-4 (or NA if no disease)
- Justification: Brief explanation
- PMIDs: Comma-separated list of cited PMIDs

If no disease is associated with {gene_name}, return: | None | NA | [3-sentence summary] | NA |
"""
        
        response = call_llm(prompt, model, temperature, agent_name="disease_agent")
        raw_attempts.append(
            {
                "attempt": attempt,
                "prompt": prompt,
                "response": response,
                "model": model,
                "temperature": temperature,
            }
        )
        parsed = parse_deep_analysis_a1(response)
        
        is_valid, error_msg = validate_a1_response(parsed)
        if is_valid:
            if attempt > 1:
                print(f"[DeepAnalysis] A1 succeeded on attempt {attempt}")
            return parsed, attempt, {"attempts": raw_attempts}
        
        print(f"[DeepAnalysis] A1 format error (attempt {attempt}): {error_msg}")
    
    # All retries failed
    print(f"[DeepAnalysis] A1 failed after {max_retries} attempts for {gene_name}")
    return [], max_retries, {"attempts": raw_attempts}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_deep_analysis(gene_name: str, cfg: Dict, base_articles: List[Dict]) -> Optional[Dict]:
    """
    Standard deep analysis with parallel LLM calls (for webapp).
    """
    try:
        return _run_deep_analysis(gene_name, cfg, base_articles)
    except LLMOverloadedError:
        raise
    except Exception as exc:  # pragma: no cover - defensive logging
        print(f"[DeepAnalysis] Unexpected failure for {gene_name}: {exc}")
        return None


def prepare_deep_analysis_pubmed(gene_name: str, cfg: Dict, base_articles: List[Dict]) -> Optional[Dict]:
    """
    Phase 1: Collect all PubMed data for a gene (sequential).
    Returns a dict with all articles and prompts ready for LLM calls.
    
    Used by scorer v3 to separate PubMed collection from LLM execution.
    In knowledge mode, skips PubMed searches and uses knowledge-based prompts.
    In simple mode, uses simplified prompts without prompt engineering.
    In proba mode, uses probabilistic prompts for penetrance, inheritance, onset_severity.
    """
    try:
        knowledge_mode = cfg.get("knowledge", False)
        simple_mode = cfg.get("simple", False)
        proba_mode = cfg.get("proba", False)
        
        # Case 1: No articles found (and not in knowledge mode) → return special result
        if not knowledge_mode and not base_articles:
            return {
                "gene_name": gene_name,
                "cfg": cfg,
                "no_articles": True,
                "algorithmic_summary": f"No relevant articles found for {gene_name}.",
            }
        
        # Disease and mechanism agents don't have proba variants
        prompt_a1 = load_prompt(_get_prompt_name("deep_analysis_disease_agent", knowledge_mode, simple_mode))
        prompt_a4 = load_prompt(_get_prompt_name("deep_analysis_mechanism_agent", knowledge_mode, simple_mode))
        # Penetrance and inheritance agents have proba variants
        prompt_a2 = load_prompt(_get_prompt_name("deep_analysis_penetrance_agent", knowledge_mode, simple_mode, proba_mode))
        prompt_a3 = load_prompt(_get_prompt_name("deep_analysis_inheritance_agent", knowledge_mode, simple_mode, proba_mode))

        if not prompt_a1:
            prompt_name = _get_prompt_name("deep_analysis_disease_agent", knowledge_mode, simple_mode)
            print(f"[DeepAnalysis] Missing prompt {prompt_name} for {gene_name}")
            return None

        model = cfg.get("model")
        temperature = cfg.get("temperature", 0.0)
        use_abstracts = cfg.get("use_abstracts", False)
        top_abstracts = cfg.get("top_abstracts")  # None = all (if use_abstracts), int = hybrid
        num_papers = cfg.get("num_papers") or 5
        max_results = max(num_papers, 1)

        base_articles_copy = list(base_articles or [])
        base_by_pmid = {art.get("pmid"): art for art in base_articles_copy if art.get("pmid")}
        
        # In knowledge mode, we don't send articles
        if knowledge_mode:
            articles_text = ""
            a1_prompt = prompt_a1.format(gene_name=gene_name)
        else:
            articles_text = articles_to_text(base_articles_copy, use_abstracts, top_abstracts)
            a1_prompt = prompt_a1.format(gene_name=gene_name, articles_text=articles_text)

        # We need to call A1 first to get disease names for other agents
        # This is unavoidable - A1 must run before we can prepare A2/A3/A4
        # Use retry logic for format errors
        a1_results, a1_attempts, a1_raw = call_a1_with_retry(a1_prompt, model, temperature, gene_name)
        a1_results, a1_attempts, a1_raw = call_a1_with_retry(a1_prompt, model, temperature, gene_name)
        
        # If A1 failed after all retries, return error result
        if not a1_results:
            print(f"[DeepAnalysis] A1 failed after {a1_attempts} attempts for {gene_name}.")
            return {
                "gene_name": gene_name,
                "cfg": cfg,
                "format_error": True,
                "error_message": f"Disease agent failed to return valid results after {a1_attempts} attempts.",
            }

        # Case 2: Disease agent found no phenotype (returned "None" with justification)
        first_disease = a1_results[0].get("disease", "").strip().lower() if a1_results else ""
        if first_disease == "none" and len(a1_results) == 1:
            justification = a1_results[0].get("justification", "No phenotype caused by mutations in this gene was identified.")
            return {
                "gene_name": gene_name,
                "cfg": cfg,
                "no_phenotype": True,
                "a1_results": a1_results,
                "justification": justification,
                "base_articles": list(base_articles or []),
            }

        disease_names = prepare_disease_names(a1_results)
        if not disease_names:
            print(f"[DeepAnalysis] No valid disease names parsed for {gene_name}.")
            return None

        cited_pmids = sorted({pmid for entry in a1_results for pmid in entry.get("pmids", [])})
        association_summaries = format_association_summaries(a1_results)
        disease_clause = build_or_clause(disease_names[:MAX_DISEASES_FOR_QUERY])

        # Collect PubMed articles for A2, A3, A4 (skip in knowledge mode)
        pubmed_data = {
            "a1": {"articles": base_articles_copy if not knowledge_mode else [], "query": "gene-wide Mendelian search" if not knowledge_mode else "knowledge-based"},
            "a2": {"articles": [], "query": None},
            "a3": {"articles": [], "query": None},
            "a4": {"articles": [], "query": None},
        }

        if not knowledge_mode:
            # A2: Penetrance
            if prompt_a2:
                query_a2 = f"{disease_clause} AND {build_keyword_clause(PENETRANCE_KEYWORDS)}"
                articles_a2 = fetch_articles_for_query(query_a2, max_results, use_abstracts)
                combined_a2 = combine_articles(articles_a2, base_by_pmid, a1_results)
                pubmed_data["a2"] = {"articles": list(combined_a2), "query": query_a2}

            # A3: Inheritance
            if prompt_a3:
                query_a3 = f'"{gene_name}" AND {disease_clause} AND {build_keyword_clause(INHERITANCE_KEYWORDS)}'
                articles_a3 = fetch_articles_for_query(query_a3, max_results, use_abstracts)
                combined_a3 = combine_articles(articles_a3, base_by_pmid, a1_results)
                pubmed_data["a3"] = {"articles": list(combined_a3), "query": query_a3}

            # A4: Mechanism
            if prompt_a4:
                query_a4 = f'"{gene_name}" AND {disease_clause} AND {build_keyword_clause(MECHANISM_KEYWORDS)}'
                articles_a4 = fetch_articles_for_query(query_a4, max_results, use_abstracts)
                combined_a4 = combine_articles(articles_a4, base_by_pmid, a1_results)
                pubmed_data["a4"] = {"articles": list(combined_a4), "query": query_a4}
        else:
            # In knowledge mode, set queries to indicate knowledge-based
            pubmed_data["a2"]["query"] = "knowledge-based"
            pubmed_data["a3"]["query"] = "knowledge-based"
            pubmed_data["a4"]["query"] = "knowledge-based"

        return {
            "gene_name": gene_name,
            "cfg": cfg,
            "a1_results": a1_results,
            "disease_names": disease_names,
            "association_summaries": association_summaries,
            "pubmed_data": pubmed_data,
            "prompts": {
                "a1": prompt_a1,
                "a2": prompt_a2,
                "a3": prompt_a3,
                "a4": prompt_a4,
            },
            "base_articles": base_articles_copy,
            "use_abstracts": use_abstracts,
            "top_abstracts": top_abstracts,
            "knowledge_mode": knowledge_mode,
        }

    except LLMOverloadedError:
        raise
    except Exception as exc:
        print(f"[DeepAnalysis] PubMed preparation failed for {gene_name}: {exc}")
        return None


def run_deep_analysis_sequential(prepared_data: Dict) -> Optional[Dict]:
    """
    Phase 2: Run LLM calls sequentially using pre-collected PubMed data.
    
    Used by scorer v3 for batch processing with gene-level parallelism.
    All LLM calls for a single gene are sequential (no internal parallelism).
    """
    if not prepared_data:
        return None

    try:
        gene_name = prepared_data["gene_name"]
        cfg = prepared_data["cfg"]
        
        # Handle special case: no articles found
        if prepared_data.get("no_articles"):
            return {
                "diseases": [],
                "algorithmic_level": None,
                "algorithmic_summary": prepared_data.get("algorithmic_summary", f"No relevant articles found for {gene_name}."),
                "no_articles": True,
                "raw": {"a1": [], "a2": [], "a3": [], "a4": []},
                "inputs": {
                    "a1": {"query": "gene-wide Mendelian search", "article_count": 0, "articles": [], "cited_pmids": []},
                    "a2": {"query": None, "article_count": 0, "articles": []},
                    "a3": {"query": None, "article_count": 0, "articles": []},
                    "a4": {"query": None, "article_count": 0, "articles": []},
                },
                "timing": {"total": 0},
            }
        
        # Handle special case: format error (A1 failed after retries)
        if prepared_data.get("format_error"):
            error_msg = prepared_data.get("error_message", "Disease agent failed to return valid results.")
            return {
                "diseases": [],
                "algorithmic_level": None,
                "algorithmic_summary": error_msg,
                "format_error": True,
                "raw": {"a1": [], "a2": [], "a3": [], "a4": []},
                "inputs": {
                    "a1": {"query": "gene-wide Mendelian search", "article_count": 0, "articles": [], "cited_pmids": []},
                    "a2": {"query": None, "article_count": 0, "articles": []},
                    "a3": {"query": None, "article_count": 0, "articles": []},
                    "a4": {"query": None, "article_count": 0, "articles": []},
                },
                "timing": {"total": 0},
            }
        
        # Handle special case: no phenotype found (disease agent returned "None")
        if prepared_data.get("no_phenotype"):
            a1_results = prepared_data.get("a1_results", [])
            justification = prepared_data.get("justification", "No phenotype caused by mutations in this gene was identified.")
            base_articles = prepared_data.get("base_articles", [])
            return {
                "diseases": [{
                    "name": "None",
                    "association_score": None,
                    "association_justification": justification,
                    "association_pmids": [],
                    "algorithmic_level": DEFAULT_ALGO_LEVEL,  # 7 = no phenotype found
                }],
                "algorithmic_level": DEFAULT_ALGO_LEVEL,  # 7 = no phenotype found
                "algorithmic_summary": justification,
                "no_phenotype": True,
                "raw": {"a1": a1_results, "a2": [], "a3": [], "a4": []},
                "inputs": {
                    "a1": {
                        "query": "gene-wide Mendelian search",
                        "article_count": len(base_articles),
                        "articles": base_articles,
                        "cited_pmids": [],
                    },
                    "a2": {"query": None, "article_count": 0, "articles": []},
                    "a3": {"query": None, "article_count": 0, "articles": []},
                    "a4": {"query": None, "article_count": 0, "articles": []},
                },
                "timing": {"total": 0},
            }
        
        a1_results = prepared_data["a1_results"]
        disease_names = prepared_data["disease_names"]
        association_summaries = prepared_data["association_summaries"]
        pubmed_data = prepared_data["pubmed_data"]
        prompts = prepared_data["prompts"]
        use_abstracts = prepared_data["use_abstracts"]
        top_abstracts = prepared_data.get("top_abstracts")

        model = cfg.get("model")
        temperature = cfg.get("temperature", 0.0)

        timing: Dict = {
            "agents": {},
            "agent_totals": {},
            "pubmed_queries": {},
            "classifiers": {},
            "other": {},
        }
        total_start = perf_counter()

        inputs = {
            "a1": {
                "query": pubmed_data["a1"]["query"],
                "article_count": len(pubmed_data["a1"]["articles"]),
                "articles": pubmed_data["a1"]["articles"],
                "cited_pmids": sorted({pmid for entry in a1_results for pmid in entry.get("pmids", [])}),
            },
            "a2": {"query": pubmed_data["a2"]["query"], "article_count": len(pubmed_data["a2"]["articles"]), "articles": pubmed_data["a2"]["articles"]},
            "a3": {"query": pubmed_data["a3"]["query"], "article_count": len(pubmed_data["a3"]["articles"]), "articles": pubmed_data["a3"]["articles"]},
            "a4": {"query": pubmed_data["a4"]["query"], "article_count": len(pubmed_data["a4"]["articles"]), "articles": pubmed_data["a4"]["articles"]},
        }

        knowledge_mode = prepared_data.get("knowledge_mode", False)
        simple_mode = cfg.get("simple", False)
        proba_mode = cfg.get("proba", False)

        # A2: Penetrance (sequential)
        a2_results = []
        a2_distributions = {}
        if prompts["a2"] and (pubmed_data["a2"]["articles"] or knowledge_mode):
            disease_list = format_disease_list(disease_names)
            if knowledge_mode:
                a2_prompt = prompts["a2"].format(gene_name=gene_name, disease_list=disease_list)
            else:
                articles_text_a2 = articles_to_text(pubmed_data["a2"]["articles"], use_abstracts, top_abstracts)
                a2_prompt = prompts["a2"].format(gene_name=gene_name, disease_list=disease_list, articles_text=articles_text_a2)
            
            a2_start = perf_counter()
            a2_response = call_llm(a2_prompt, model, temperature, agent_name="penetrance_agent")
            timing["agents"]["a2"] = perf_counter() - a2_start
            if proba_mode:
                from parsers.llm_parsers import parse_penetrance_proba_response
                a2_distributions = parse_penetrance_proba_response(a2_response) if a2_response else {}
                # Convert distributions to expected scores for compatibility
                a2_results = _distributions_to_scores_penetrance(a2_distributions)
            else:
                a2_results = parse_deep_analysis_a2(a2_response) if a2_response else []

        # A3: Inheritance (sequential)
        a3_results = []
        a3_distributions = {}
        if prompts["a3"] and (pubmed_data["a3"]["articles"] or knowledge_mode):
            disease_list = format_disease_list(disease_names)
            if knowledge_mode:
                a3_prompt = prompts["a3"].format(gene_name=gene_name, disease_list=disease_list)
            else:
                articles_text_a3 = articles_to_text(pubmed_data["a3"]["articles"], use_abstracts, top_abstracts)
                a3_prompt = prompts["a3"].format(gene_name=gene_name, disease_list=disease_list, articles_text=articles_text_a3)
            
            a3_start = perf_counter()
            a3_response = call_llm(a3_prompt, model, temperature, agent_name="inheritance_agent")
            timing["agents"]["a3"] = perf_counter() - a3_start
            if proba_mode:
                from parsers.llm_parsers import parse_inheritance_proba_response
                a3_distributions = parse_inheritance_proba_response(a3_response) if a3_response else {}
                # Convert distributions to expected scores for compatibility
                a3_results = _distributions_to_scores_inheritance(a3_distributions)
            else:
                a3_results = parse_deep_analysis_a3(a3_response) if a3_response else []

        # A4: Mechanism (sequential)
        a4_results = []
        if prompts["a4"] and (pubmed_data["a4"]["articles"] or knowledge_mode):
            disease_list = format_disease_list(disease_names)
            if knowledge_mode:
                a4_prompt = prompts["a4"].format(gene_name=gene_name, disease_list=disease_list)
            else:
                articles_text_a4 = articles_to_text(pubmed_data["a4"]["articles"], use_abstracts, top_abstracts)
                a4_prompt = prompts["a4"].format(gene_name=gene_name, disease_list=disease_list, articles_text=articles_text_a4)
            
            a4_start = perf_counter()
            a4_response = call_llm(a4_prompt, model, temperature, agent_name="mechanism_agent")
            timing["agents"]["a4"] = perf_counter() - a4_start
            a4_results = parse_deep_analysis_a4(a4_response) if a4_response else []

        # Pathogenicity classifier (sequential)
        pathogenicity_map = {}
        if disease_names and association_summaries:
            path_prompt_template = load_prompt(_get_prompt_name("deep_analysis_pathogenicity_agent", simple_mode=simple_mode))
            if path_prompt_template:
                diseases_list = format_disease_list(disease_names)
                path_prompt = path_prompt_template.format(
                    gene_name=gene_name,
                    disease_list=diseases_list,
                    association_summaries=association_summaries,
                )
                path_start = perf_counter()
                path_response = call_llm(path_prompt, model, temperature, agent_name="pathogenicity_agent")
                timing["classifiers"]["pathogenicity"] = perf_counter() - path_start
                pathogenicity_map = parse_pathogenicity_response(path_response) if path_response else {}

        # Onset/severity classifier (sequential)
        onset_severity_map = {}
        onset_severity_distributions = {}
        if disease_names:
            onset_prompt_template = load_prompt(_get_prompt_name("deep_analysis_onset_severity_agent", simple_mode=simple_mode, proba_mode=proba_mode))
            if onset_prompt_template:
                diseases_list = format_disease_list(disease_names)
                onset_prompt = onset_prompt_template.format(disease_list=diseases_list, gene_name=gene_name)
                onset_start = perf_counter()
                onset_response = call_llm(onset_prompt, model, temperature, agent_name="onset_severity_agent")
                timing["classifiers"]["onset_severity"] = perf_counter() - onset_start
                if proba_mode:
                    from parsers.llm_parsers import parse_onset_severity_proba_response
                    onset_severity_distributions = parse_onset_severity_proba_response(onset_response) if onset_response else {}
                    # Convert distributions to expected scores for compatibility
                    onset_severity_map = _distributions_to_scores_onset_severity(onset_severity_distributions)
                else:
                    onset_severity_map = parse_onset_severity_response(onset_response) if onset_response else {}

        inputs["pathogenicity"] = {"classification": pathogenicity_map}

        # Aggregate results
        aggregate_start = perf_counter()
        aggregate, min_algo_level = _aggregate_deep_analysis(
            disease_names, a1_results, a2_results, a3_results, a4_results, onset_severity_map
        )
        timing["other"]["aggregate"] = perf_counter() - aggregate_start
        _apply_pathogenicity_map(aggregate, pathogenicity_map)

        # Algorithmic summary (sequential)
        algo_summary_start = perf_counter()
        algorithmic_summary = _run_algorithmic_summary_agent(gene_name, cfg, aggregate)
        timing["classifiers"]["algorithmic_summary"] = perf_counter() - algo_summary_start

        timing["total"] = perf_counter() - total_start

        # Add distributions if in proba mode
        proba_metadata = None
        if proba_mode:
            # Apply distributions to disease records
            _apply_distributions_to_diseases(aggregate, a2_distributions, a3_distributions, onset_severity_distributions)
            # Compute Monte Carlo distribution for algorithmic level
            _compute_monte_carlo_levels(aggregate)
            # Use expected_level from Monte Carlo as the algorithmic_level (decimal)
            proba_metadata = _get_min_expected_level_with_metadata(aggregate)
            if proba_metadata["expected_level"] is not None:
                min_algo_level = proba_metadata["expected_level"]
        
        result = {
            "diseases": aggregate,
            "algorithmic_level": min_algo_level,
            "algorithmic_summary": algorithmic_summary,
            "raw": {
                "a1": a1_results,
                "a2": a2_results,
                "a3": a3_results,
                "a4": a4_results,
            },
            "inputs": inputs,
            "pathogenicity_agent": pathogenicity_map,
            "timing": timing,
        }
        
        # Add distributions metadata if in proba mode
        if proba_mode:
            result["proba_mode"] = True
            result["distributions"] = {
                "penetrance": a2_distributions,
                "inheritance": a3_distributions,
                "onset_severity": onset_severity_distributions,
            }
            # Add Monte Carlo metadata at deep_analysis level
            if proba_metadata:
                result["expected_level"] = proba_metadata["expected_level"]
                result["level_variance"] = proba_metadata["level_variance"]
                result["level_distribution"] = proba_metadata["level_distribution"]
                result["kappa"] = proba_metadata["kappa"]
        
        return result

    except LLMOverloadedError:
        raise
    except Exception as exc:
        print(f"[DeepAnalysis] Sequential execution failed: {exc}")
        return None


# ---------------------------------------------------------------------------
# Pipeline orchestration (original with parallel LLM calls for webapp)
# ---------------------------------------------------------------------------

def _run_deep_analysis(gene_name: str, cfg: Dict, base_articles: List[Dict]) -> Optional[Dict]:
    timing: Dict = {
        "agents": {},
        "agent_totals": {},
        "pubmed_queries": {},
        "classifiers": {},
        "other": {},
    }
    total_start = perf_counter()
    def finalize_timing():
        if "total" not in timing:
            timing["total"] = perf_counter() - total_start

    knowledge_mode = cfg.get("knowledge", False)
    simple_mode = cfg.get("simple", False)
    proba_mode = cfg.get("proba", False)

    # Case 1: No articles found (and not in knowledge mode) → return special message
    if not knowledge_mode and not base_articles:
        finalize_timing()
        return {
            "diseases": [],
            "algorithmic_level": None,
            "algorithmic_summary": f"No relevant articles found for {gene_name}.",
            "no_articles": True,
            "raw": {"a1": [], "a2": [], "a3": [], "a4": []},
            "inputs": {
                "a1": {"query": "gene-wide Mendelian search", "article_count": 0, "articles": [], "cited_pmids": []},
                "a2": {"query": None, "article_count": 0, "articles": []},
                "a3": {"query": None, "article_count": 0, "articles": []},
                "a4": {"query": None, "article_count": 0, "articles": []},
            },
            "timing": timing,
        }

    # Disease and mechanism agents don't have proba variants
    prompt_a1 = load_prompt(_get_prompt_name("deep_analysis_disease_agent", knowledge_mode, simple_mode))
    prompt_a4 = load_prompt(_get_prompt_name("deep_analysis_mechanism_agent", knowledge_mode, simple_mode))
    # Penetrance and inheritance agents have proba variants
    prompt_a2 = load_prompt(_get_prompt_name("deep_analysis_penetrance_agent", knowledge_mode, simple_mode, proba_mode))
    prompt_a3 = load_prompt(_get_prompt_name("deep_analysis_inheritance_agent", knowledge_mode, simple_mode, proba_mode))

    if not prompt_a1:
        prompt_name = _get_prompt_name("deep_analysis_disease_agent", knowledge_mode, simple_mode)
        print(f"[DeepAnalysis] Missing prompt {prompt_name}")
        finalize_timing()
        return None

    model = cfg.get("model")
    temperature = cfg.get("temperature", 0.0)
    use_abstracts = cfg.get("use_abstracts", False)
    top_abstracts = cfg.get("top_abstracts")  # None = all (if use_abstracts), int = hybrid
    num_papers = cfg.get("num_papers") or 5
    max_results = max(num_papers, 1)

    base_articles_copy = list(base_articles or [])
    base_by_pmid = {art.get("pmid"): art for art in base_articles_copy if art.get("pmid")}
    
    if knowledge_mode:
        articles_text = ""
    else:
        articles_text = articles_to_text(base_articles_copy, use_abstracts, top_abstracts)

    inputs = {
        "a1": {
            "query": "knowledge-based" if knowledge_mode else "gene-wide Mendelian search",
            "article_count": 0 if knowledge_mode else len(base_articles_copy),
            "articles": [] if knowledge_mode else base_articles_copy,
            "cited_pmids": [],
        },
        "a2": {"query": None, "article_count": 0, "articles": []},
        "a3": {"query": None, "article_count": 0, "articles": []},
        "a4": {"query": None, "article_count": 0, "articles": []},
    }

    if knowledge_mode:
        a1_prompt = prompt_a1.format(gene_name=gene_name)
    else:
        a1_prompt = prompt_a1.format(gene_name=gene_name, articles_text=articles_text)
    
    # Call A1 with retry logic for format errors
    a1_start = perf_counter()
    a1_results, a1_attempts, a1_raw = call_a1_with_retry(a1_prompt, model, temperature, gene_name)
    timing["agents"]["a1"] = perf_counter() - a1_start
    timing["agent_totals"]["a1"] = timing["agents"]["a1"]
    
    # If A1 failed after all retries, return error result
    if not a1_results:
        print(f"[DeepAnalysis] A1 failed after {a1_attempts} attempts; returning error.")
        finalize_timing()
        return {
            "diseases": [],
            "algorithmic_level": None,
            "algorithmic_summary": f"Disease agent failed to return valid results after {a1_attempts} attempts.",
            "format_error": True,
            "raw": {"a1": [], "a2": [], "a3": [], "a4": []},
            "inputs": {
                "a1": {
                    "query": "knowledge-based" if knowledge_mode else "gene-wide Mendelian search",
                    "article_count": 0 if knowledge_mode else len(base_articles_copy),
                    "articles": [] if knowledge_mode else base_articles_copy,
                    "cited_pmids": [],
                },
                "a2": {"query": None, "article_count": 0, "articles": []},
                "a3": {"query": None, "article_count": 0, "articles": []},
                "a4": {"query": None, "article_count": 0, "articles": []},
            },
            "timing": timing,
        }

    # Case 2: Disease agent found no phenotype (returned "None" with justification)
    first_disease = a1_results[0].get("disease", "").strip().lower() if a1_results else ""
    if first_disease == "none" and len(a1_results) == 1:
        justification = a1_results[0].get("justification", "No phenotype caused by mutations in this gene was identified.")
        finalize_timing()
        return {
            "diseases": [{
                "name": "None",
                "association_score": None,
                "association_justification": justification,
                "association_pmids": [],
                "algorithmic_level": DEFAULT_ALGO_LEVEL,  # 7 = no phenotype found
            }],
            "algorithmic_level": DEFAULT_ALGO_LEVEL,  # 7 = no phenotype found
            "algorithmic_summary": justification,
            "no_phenotype": True,
            "raw": {"a1": a1_results, "a2": [], "a3": [], "a4": []},
            "inputs": {
                "a1": {
                    "query": "knowledge-based" if knowledge_mode else "gene-wide Mendelian search",
                    "article_count": 0 if knowledge_mode else len(base_articles_copy),
                    "articles": [] if knowledge_mode else base_articles_copy,
                    "cited_pmids": [],
                },
                "a2": {"query": None, "article_count": 0, "articles": []},
                "a3": {"query": None, "article_count": 0, "articles": []},
                "a4": {"query": None, "article_count": 0, "articles": []},
            },
            "timing": timing,
        }

    disease_names = prepare_disease_names(a1_results)
    if not disease_names:
        print("[DeepAnalysis] No valid disease names parsed.")
        finalize_timing()
        return None

    cited_pmids = sorted({pmid for entry in a1_results for pmid in entry.get("pmids", [])})
    inputs["a1"]["cited_pmids"] = cited_pmids

    association_summaries = format_association_summaries(a1_results)

    disease_clause = build_or_clause(disease_names[:MAX_DISEASES_FOR_QUERY])

    agent_tasks = []
    a2_task = _prepare_agent_task(
        prompt=prompt_a2,
        query=f"{disease_clause} AND {build_keyword_clause(PENETRANCE_KEYWORDS)}"
        if prompt_a2 and not knowledge_mode
        else None,
        disease_names=disease_names,
        gene_name=gene_name,
        model=model,
        temperature=temperature,
        base_by_pmid=base_by_pmid,
        base_results=a1_results,
        max_results=max_results,
        use_abstracts=use_abstracts,
        top_abstracts=top_abstracts,
        parsing_fn_name="a2_proba" if proba_mode else "a2",
        inputs_section=inputs["a2"],
        timing=timing,
        agent_key="a2",
        knowledge_mode=knowledge_mode,
    )
    if a2_task:
        agent_tasks.append(a2_task)

    a3_task = _prepare_agent_task(
        prompt=prompt_a3,
        query=f'"{gene_name}" AND {disease_clause} AND {build_keyword_clause(INHERITANCE_KEYWORDS)}'
        if prompt_a3 and not knowledge_mode
        else None,
        disease_names=disease_names,
        gene_name=gene_name,
        model=model,
        temperature=temperature,
        base_by_pmid=base_by_pmid,
        base_results=a1_results,
        max_results=max_results,
        use_abstracts=use_abstracts,
        top_abstracts=top_abstracts,
        parsing_fn_name="a3_proba" if proba_mode else "a3",
        inputs_section=inputs["a3"],
        timing=timing,
        agent_key="a3",
        knowledge_mode=knowledge_mode,
    )
    if a3_task:
        agent_tasks.append(a3_task)

    a4_task = _prepare_agent_task(
        prompt=prompt_a4,
        query=f'"{gene_name}" AND {disease_clause} AND {build_keyword_clause(MECHANISM_KEYWORDS)}'
        if prompt_a4 and not knowledge_mode
        else None,
        disease_names=disease_names,
        gene_name=gene_name,
        model=model,
        temperature=temperature,
        base_by_pmid=base_by_pmid,
        base_results=a1_results,
        max_results=max_results,
        use_abstracts=use_abstracts,
        top_abstracts=top_abstracts,
        parsing_fn_name="a4",
        inputs_section=inputs["a4"],
        timing=timing,
        agent_key="a4",
        knowledge_mode=knowledge_mode,
    )
    if a4_task:
        agent_tasks.append(a4_task)

    classifier_tasks = _prepare_classifier_tasks(
        gene_name=gene_name,
        disease_names=disease_names,
        association_summaries=association_summaries,
        cfg=cfg,
    )

    parallel_results = _run_parallel_agents(agent_tasks + classifier_tasks, timing)
    raw_parallel = parallel_results.pop("_raw_llm_parallel", {}) if isinstance(parallel_results, dict) else {}
    # Collect raw prompts/responses from parallel results (attached by _execute_agent_task)
    raw_llm = {
        "a1": a1_raw,
        "parallel": raw_parallel,
    }
    
    # Handle proba mode: convert distributions to scores for compatibility
    a2_distributions = {}
    a3_distributions = {}
    onset_severity_distributions = {}
    
    if proba_mode:
        # In proba mode, results are distributions (dicts)
        a2_distributions = parallel_results.get("a2") or {}
        a3_distributions = parallel_results.get("a3") or {}
        onset_severity_distributions = parallel_results.get("onset_severity") or {}
        
        # Convert distributions to score results for pipeline compatibility
        a2_results = _distributions_to_scores_penetrance(a2_distributions)
        a3_results = _distributions_to_scores_inheritance(a3_distributions)
        onset_severity_map = _distributions_to_scores_onset_severity(onset_severity_distributions)
    else:
        a2_results = parallel_results.get("a2", [])
        a3_results = parallel_results.get("a3", [])
        onset_severity_map = parallel_results.get("onset_severity") or {}
    
    a4_results = parallel_results.get("a4", [])
    pathogenicity_map = parallel_results.get("pathogenicity") or {}
    inputs["pathogenicity"] = {"classification": pathogenicity_map}

    aggregate_start = perf_counter()
    aggregate, min_algo_level = _aggregate_deep_analysis(
        disease_names,
        a1_results,
        a2_results,
        a3_results,
        a4_results,
        onset_severity_map,
    )
    timing["other"]["aggregate"] = perf_counter() - aggregate_start
    _apply_pathogenicity_map(aggregate, pathogenicity_map)
    
    # Apply distributions to diseases in proba mode
    proba_metadata = None
    if proba_mode:
        _apply_distributions_to_diseases(aggregate, a2_distributions, a3_distributions, onset_severity_distributions)
        # Compute Monte Carlo distribution for algorithmic level
        _compute_monte_carlo_levels(aggregate)
        # Use expected_level from Monte Carlo as the algorithmic_level (decimal)
        proba_metadata = _get_min_expected_level_with_metadata(aggregate)
        if proba_metadata["expected_level"] is not None:
            min_algo_level = proba_metadata["expected_level"]

    algo_summary_start = perf_counter()
    algorithmic_summary, algo_summary_raw = _run_algorithmic_summary_agent_with_raw(gene_name, cfg, aggregate)
    timing["classifiers"]["algorithmic_summary"] = perf_counter() - algo_summary_start

    finalize_timing()

    raw_llm["algorithmic_summary"] = algo_summary_raw

    result = {
        "diseases": aggregate,
        "algorithmic_level": min_algo_level,
        "algorithmic_summary": algorithmic_summary,
        "_raw_llm": raw_llm,
        "raw": {
            "a1": a1_results,
            "a2": a2_results,
            "a3": a3_results,
            "a4": a4_results,
        },
        "inputs": inputs,
        "pathogenicity_agent": pathogenicity_map,
        "timing": timing,
    }
    
    # Add distributions if in proba mode
    if proba_mode:
        result["proba_mode"] = True
        result["distributions"] = {
            "penetrance": a2_distributions,
            "inheritance": a3_distributions,
            "onset_severity": onset_severity_distributions,
        }
        # Add Monte Carlo metadata at deep_analysis level
        if proba_metadata:
            result["expected_level"] = proba_metadata["expected_level"]
            result["level_variance"] = proba_metadata["level_variance"]
            result["level_distribution"] = proba_metadata["level_distribution"]
            result["kappa"] = proba_metadata["kappa"]
    
    return result


def _prepare_agent_task(
    *,
    prompt,
    query,
    disease_names,
    gene_name,
    model,
    temperature,
    base_by_pmid,
    base_results,
    max_results,
    use_abstracts,
    top_abstracts=None,
    parsing_fn_name,
    inputs_section,
    timing=None,
    agent_key=None,
    knowledge_mode=False,
):
    if not prompt:
        return None
    
    # In knowledge mode, we don't need query or PubMed search
    if knowledge_mode:
        inputs_section.update(
            {
                "query": "knowledge-based",
                "article_count": 0,
                "articles": [],
            }
        )
        
        disease_list = format_disease_list(disease_names)
        agent_prompt = prompt.format(
            gene_name=gene_name,
            disease_list=disease_list,
        )
        
        return {
            "agent_key": agent_key,
            "prompt": agent_prompt,
            "model": model,
            "temperature": temperature,
            "parsing_fn_name": parsing_fn_name,
            "fetch_duration": 0.0,
            "timing_bucket": "agents",
        }
    
    # Normal mode with PubMed search
    if not query:
        return None

    fetch_start = perf_counter()
    articles = fetch_articles_for_query(query, max_results, use_abstracts)
    fetch_duration = perf_counter() - fetch_start
    if timing is not None and agent_key:
        timing["pubmed_queries"][agent_key] = fetch_duration

    combined = combine_articles(articles, base_by_pmid, base_results)

    inputs_section.update(
        {
            "query": query,
            "article_count": len(combined),
            "articles": list(combined),
        }
    )

    articles_text = articles_to_text(combined, use_abstracts, top_abstracts)
    disease_list = format_disease_list(disease_names)

    agent_prompt = prompt.format(
        gene_name=gene_name,
        disease_list=disease_list,
        articles_text=articles_text,
    )

    return {
        "agent_key": agent_key,
        "prompt": agent_prompt,
        "model": model,
        "temperature": temperature,
        "parsing_fn_name": parsing_fn_name,
        "fetch_duration": fetch_duration,
        "timing_bucket": "agents",
    }


def _run_parallel_agents(tasks, timing):
    if not tasks:
        return {}

    results = {}
    raw_payloads = {}
    tracker = get_tracker()
    workers = max(1, len(tasks))
    # Agent tasks are pure I/O-bound LLM calls (_execute_agent_task), so threads
    # are the right tool: they share the worker's memory instead of spawning
    # extra processes that each re-import this module and reload the big
    # reference files (GenCC/GoF/gene scores) -> that nested process fan-out was
    # blowing up RAM under high --n_core. Override with DEEP_ANALYSIS_EXECUTOR=process.
    _executor_kind = os.environ.get("DEEP_ANALYSIS_EXECUTOR", "thread").lower()
    Executor = ProcessPoolExecutor if _executor_kind == "process" else ThreadPoolExecutor
    with Executor(max_workers=workers) as executor:
        future_map = {executor.submit(_execute_agent_task, task): task for task in tasks}
        for future in as_completed(future_map):
            task = future_map[future]
            agent_key = task.get("agent_key")
            try:
                output = future.result()
            except Exception as exc:  # pragma: no cover
                print(f"[DeepAnalysis] Parallel agent failure ({agent_key}): {exc}")
                output = {"result": [], "llm_duration": 0.0, "input_tokens": 0, "output_tokens": 0}

            results[agent_key] = output.get("result") or []
            # Capture raw artifacts when available (prompt + response + tokens)
            if agent_key:
                raw_payloads[agent_key] = {
                    "agent_name": output.get("agent_name", agent_key),
                    "model": output.get("model", task.get("model", "unknown")),
                    "temperature": task.get("temperature"),
                    "prompt": output.get("prompt", task.get("prompt")),
                    "response": output.get("response"),
                    "input_tokens": output.get("input_tokens", 0),
                    "output_tokens": output.get("output_tokens", 0),
                    "llm_duration": output.get("llm_duration", 0.0),
                }
            
            # Record tokens in the main process tracker
            agent_name = output.get("agent_name", agent_key)
            model = output.get("model", task.get("model", "unknown"))
            input_tokens = output.get("input_tokens", 0)
            output_tokens = output.get("output_tokens", 0)
            if input_tokens > 0 or output_tokens > 0:
                tracker.record(agent_name, model, input_tokens, output_tokens)
            
            if timing is not None and agent_key:
                llm_duration = output.get("llm_duration", 0.0)
                fetch_duration = task.get("fetch_duration", 0.0)
                bucket = task.get("timing_bucket", "agents")
                if bucket == "agents":
                    timing["agents"][agent_key] = llm_duration
                    timing["agent_totals"][agent_key] = fetch_duration + llm_duration
                elif bucket == "classifiers":
                    timing["classifiers"][agent_key] = llm_duration

    # Surface raw artifacts under a reserved key (used internally for persistence only).
    results["_raw_llm_parallel"] = raw_payloads
    return results


PARSING_FUNCTIONS: Dict[str, Callable] = {
    "a2": parse_deep_analysis_a2,
    "a3": parse_deep_analysis_a3,
    "a4": parse_deep_analysis_a4,
    "pathogenicity": parse_pathogenicity_response,
    "onset_severity": parse_onset_severity_response,
}

# Lazy import of proba parsers to avoid circular imports
def _get_proba_parsing_functions():
    from parsers.llm_parsers import (
        parse_penetrance_proba_response,
        parse_inheritance_proba_response,
        parse_onset_severity_proba_response,
    )
    return {
        "a2_proba": parse_penetrance_proba_response,
        "a3_proba": parse_inheritance_proba_response,
        "onset_severity_proba": parse_onset_severity_proba_response,
    }


AGENT_KEY_TO_NAME = {
    "a2": "penetrance_agent",
    "a3": "inheritance_agent",
    "a4": "mechanism_agent",
    "pathogenicity": "pathogenicity_agent",
    "onset_severity": "onset_severity_agent",
}


def _execute_agent_task(task):
    prompt = task["prompt"]
    model = task["model"]
    temperature = task["temperature"]
    parsing_fn_name = task["parsing_fn_name"]
    agent_key = task.get("agent_key")
    
    # Check if this is a proba parser
    parsing_fn = PARSING_FUNCTIONS.get(parsing_fn_name)
    if parsing_fn is None and parsing_fn_name.endswith("_proba"):
        proba_parsers = _get_proba_parsing_functions()
        parsing_fn = proba_parsers.get(parsing_fn_name)
    
    if parsing_fn is None:
        return {"agent_key": agent_key, "result": [], "llm_duration": 0.0, "input_tokens": 0, "output_tokens": 0}

    # Map agent_key to a readable agent name for token tracking
    agent_name = AGENT_KEY_TO_NAME.get(agent_key, agent_key)

    llm_start = perf_counter()
    input_tokens = 0
    output_tokens = 0
    response = None
    try:
        # Use call_llm_with_usage to get tokens (needed for multiprocessing)
        # Don't pass agent_name here - we'll record in parent process
        llm_result = call_llm_with_usage(prompt, model, temperature)
        response = llm_result["text"]
        input_tokens = llm_result["input_tokens"]
        output_tokens = llm_result["output_tokens"]
        parsed = parsing_fn(response)
    except Exception as exc:  # pragma: no cover - defensive logging
        print(f"[DeepAnalysis] Agent execution failed: {exc}")
        parsed = []
    llm_duration = perf_counter() - llm_start
    return {
        "agent_key": agent_key,
        "agent_name": agent_name,
        "model": model,
        "result": parsed,
        "prompt": prompt,
        "response": response,
        "llm_duration": llm_duration,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
    }


def _run_algorithmic_summary_agent_with_raw(gene_name: str, cfg: Dict, records: List[Dict]) -> Tuple[Optional[str], Dict]:
    """
    Wrapper that returns both parsed algorithmic summary and raw artifacts.
    """
    prompt_template = load_prompt(_get_prompt_name("deep_analysis_algorithmic_summary_agent", simple_mode=cfg.get("simple", False)))
    if not prompt_template:
        return None, {"prompt": None, "response": None, "model": cfg.get("model"), "temperature": cfg.get("temperature", 0.0)}

    selected_records = list(records or [])[:MAX_DISEASES_FOR_QUERY]
    disease_blocks = []
    for rec in selected_records:
        disease_blocks.append(_format_algorithmic_summary_block(rec))

    prompt = prompt_template.format(
        gene_name=gene_name,
        disease_blocks="\n\n".join(block for block in disease_blocks if block.strip()) or "No disease data available."
    )

    try:
        response = call_llm(
            prompt,
            cfg["model"],
            cfg.get("temperature", 0.0),
            agent_name="algorithmic_summary_agent",
        )
        parsed = parse_algorithmic_summary_response(response)
        return parsed, {"prompt": prompt, "response": response, "model": cfg.get("model"), "temperature": cfg.get("temperature", 0.0)}
    except Exception as exc:  # pragma: no cover
        print(f"[DeepAnalysis] Algorithmic summary agent failed: {exc}")
        return None, {"prompt": prompt, "response": None, "model": cfg.get("model"), "temperature": cfg.get("temperature", 0.0)}


# ---------------------------------------------------------------------------
# Supporting classifiers
# ---------------------------------------------------------------------------

def _prepare_classifier_tasks(gene_name, disease_names, association_summaries, cfg):
    tasks = []
    model = cfg["model"]
    temperature = cfg.get("temperature", 0.0)
    simple_mode = cfg.get("simple", False)
    proba_mode = cfg.get("proba", False)

    if disease_names and association_summaries:
        path_prompt_template = load_prompt(_get_prompt_name("deep_analysis_pathogenicity_agent", simple_mode=simple_mode))
        if path_prompt_template:
            diseases_list = format_disease_list(disease_names)
            prompt = path_prompt_template.format(
                gene_name=gene_name,
                disease_list=diseases_list,
                association_summaries=association_summaries,
            )
            tasks.append(
                {
                    "agent_key": "pathogenicity",
                    "prompt": prompt,
                    "model": model,
                    "temperature": temperature,
                    "parsing_fn_name": "pathogenicity",
                    "timing_bucket": "classifiers",
                    "fetch_duration": 0.0,
                }
            )

    if disease_names:
        # Onset/severity agent has proba variant
        onset_prompt_template = load_prompt(_get_prompt_name("deep_analysis_onset_severity_agent", simple_mode=simple_mode, proba_mode=proba_mode))
        if onset_prompt_template:
            diseases_list = format_disease_list(disease_names)
            prompt = onset_prompt_template.format(disease_list=diseases_list, gene_name=gene_name)
            tasks.append(
                {
                    "agent_key": "onset_severity",
                    "prompt": prompt,
                    "model": model,
                    "temperature": temperature,
                    "parsing_fn_name": "onset_severity_proba" if proba_mode else "onset_severity",
                    "timing_bucket": "classifiers",
                    "fetch_duration": 0.0,
                }
            )

    return tasks


# ---------------------------------------------------------------------------
# Aggregation & scoring
# ---------------------------------------------------------------------------

def _aggregate_deep_analysis(disease_names, a1, a2, a3, a4, onset_severity_map=None):
    records = {}

    def ensure(name):
        if name not in records:
            records[name] = {
                "name": name,
                "association_score": None,
                "association_justification": None,
                "association_pmids": [],
                "association_is_protective": False,
                "association_is_neutral": False,
                "penetrance_level": None,
                "penetrance_justification": None,
                "penetrance_pmids": [],
                "inheritance": None,
                "inheritance_score": None,
                "inheritance_justification": None,
                "inheritance_pmids": [],
                "mechanism": None,
                "mechanism_confidence": None,
                "mechanism_justification": None,
                "mechanism_pmids": [],
                "disease_onset": None,
                "disease_onset_score": None,
                "severity_score": None,
                "severity_justification": None,
            }
        return records[name]

    for entry in a1:
        name = entry.get("disease")
        if not name:
            continue
        rec = ensure(name)
        rec["association_score"] = entry.get("association_score")
        rec["association_justification"] = entry.get("justification")
        rec["association_pmids"] = entry.get("pmids", [])

    for entry in a2:
        name = entry.get("disease")
        if not name:
            continue
        rec = ensure(name)
        rec["penetrance_level"] = entry.get("penetrance_level")
        rec["penetrance_justification"] = entry.get("justification")
        rec["penetrance_pmids"] = entry.get("pmids", [])

    for entry in a3:
        name = entry.get("disease")
        if not name:
            continue
        rec = ensure(name)
        rec["inheritance_score"] = entry.get("inheritance_score")
        rec["inheritance"] = entry.get("inheritance") or _label_from_inheritance_score(
            entry.get("inheritance_score")
        )
        rec["inheritance_justification"] = entry.get("justification")
        rec["inheritance_pmids"] = entry.get("pmids", [])

    for entry in a4:
        name = entry.get("disease")
        if not name:
            continue
        rec = ensure(name)
        rec["mechanism"] = entry.get("mechanism")
        rec["mechanism_confidence"] = entry.get("confidence_score")
        rec["mechanism_justification"] = entry.get("justification")
        rec["mechanism_pmids"] = entry.get("pmids", [])

    for entry_list in (a2, a3, a4):
        for entry in entry_list:
            name = entry.get("disease")
            if name:
                ensure(name)

    if onset_severity_map:
        for name, payload in onset_severity_map.items():
            if not name:
                continue
            rec = ensure(name)
            onset_score = payload.get("onset_score")
            rec["disease_onset_score"] = onset_score
            rec["disease_onset"] = _label_from_onset_score(onset_score) or payload.get("onset")
            rec["severity_score"] = payload.get("severity_score")
            rec["severity_justification"] = payload.get("justification")

    ordered_records = []
    seen = set()
    for name in disease_names + [n for n in records.keys() if n not in disease_names]:
        if name not in seen:
            seen.add(name)
            ordered_records.append(records[name])

    ordered_records.sort(key=lambda item: (item.get("association_score") is None, item.get("association_score") or 99))

    algo_min = None
    for rec in ordered_records:
        if rec.get("association_is_protective") or rec.get("association_is_neutral"):
            algo_level = DEFAULT_ALGO_LEVEL
        else:
            algo_level = _compute_algorithmic_level(rec)
        rec["algorithmic_level"] = algo_level
        if algo_min is None or algo_level < algo_min:
            algo_min = algo_level

    # Reorder diseases by evidence first, then algorithmic level
    ordered_records.sort(
        key=lambda rec: (
            rec.get("association_score") is None,
            rec.get("association_score") if rec.get("association_score") is not None else float("inf"),
            rec.get("algorithmic_level") is None,
            rec.get("algorithmic_level") if rec.get("algorithmic_level") is not None else float("inf"),
        )
    )

    return ordered_records, algo_min


# ---------------------------------------------------------------------------
# Score calculators & helpers
# ---------------------------------------------------------------------------

def _apply_pathogenicity_map(records, pathogenicity_map):
    if not pathogenicity_map:
        return
    for rec in records:
        name = rec.get("name")
        classification = pathogenicity_map.get(name) if name else None
        rec["association_is_protective"] = classification == "protective"
        rec["association_is_neutral"] = classification == "neutral"


# ---------------------------------------------------------------------------
# Probabilistic mode helpers
# ---------------------------------------------------------------------------

def _distribution_to_expected(distribution: dict, values: list) -> float:
    """
    Convert a probability distribution to expected value (mean).
    
    Args:
        distribution: Dict mapping category names to probabilities (0-100)
        values: List of numeric values corresponding to each category
    
    Returns:
        Expected value (weighted mean)
    """
    if not distribution:
        return None
    
    total_prob = sum(distribution.values())
    if total_prob == 0:
        return None
    
    # Normalize probabilities to sum to 1
    probs = [p / total_prob for p in distribution.values()]
    expected = sum(v * p for v, p in zip(values, probs))
    return expected


def _distribution_variance(distribution: dict, values: list) -> float:
    """
    Compute variance of a probability distribution.
    
    Args:
        distribution: Dict mapping category names to probabilities (0-100)
        values: List of numeric values corresponding to each category
    
    Returns:
        Variance of the distribution
    """
    if not distribution:
        return None
    
    total_prob = sum(distribution.values())
    if total_prob == 0:
        return None
    
    probs = [p / total_prob for p in distribution.values()]
    expected = sum(v * p for v, p in zip(values, probs))
    variance = sum(p * (v - expected) ** 2 for v, p in zip(values, probs))
    return variance


def _distributions_to_scores_penetrance(distributions: dict) -> list:
    """
    Convert penetrance distributions to score results compatible with existing pipeline.
    
    Args:
        distributions: Dict mapping disease names to distribution dicts
            New format: {"Huntington": {"distribution": {...}, "justification": "...", "pmids": [...]}}
            Legacy format: {"Huntington": {"mendelian": 85, ...}}
    
    Returns:
        List of dicts with disease name and penetrance_level (expected value rounded)
    """
    results = []
    values = [1, 2, 3, 4]  # Score values for mendelian, high, moderate, complex
    
    for disease_name, data in distributions.items():
        # Handle both new format (with "distribution" key) and legacy format
        if "distribution" in data:
            dist = data["distribution"]
            justification = data.get("justification")
            pmids = data.get("pmids", [])
        else:
            dist = data
            justification = None
            pmids = []
        
        expected = _distribution_to_expected(dist, values)
        score = round(expected) if expected is not None else None
        results.append({
            "disease": disease_name,
            "penetrance_level": score,
            "penetrance_distribution": dist,
            "penetrance_expected": expected,
            "penetrance_variance": _distribution_variance(dist, values),
            "justification": justification,
            "pmids": pmids,
        })
    
    return results


def _distributions_to_scores_inheritance(distributions: dict) -> list:
    """
    Convert inheritance distributions to score results compatible with existing pipeline.
    
    Args:
        distributions: Dict mapping disease names to distribution dicts
            New format: {"Huntington": {"distribution": {...}, "justification": "...", "pmids": [...]}}
            Legacy format: {"Huntington": {"dominant": 90, ...}}
    
    Returns:
        List of dicts with disease name and inheritance_score (expected value rounded)
    """
    results = []
    values = [1, 2, 3, 4, 5, 6]  # Score values for dominant through recessive
    
    for disease_name, data in distributions.items():
        # Handle both new format (with "distribution" key) and legacy format
        if "distribution" in data:
            dist = data["distribution"]
            justification = data.get("justification")
            pmids = data.get("pmids", [])
        else:
            dist = data
            justification = None
            pmids = []
        
        expected = _distribution_to_expected(dist, values)
        score = round(expected) if expected is not None else None
        results.append({
            "disease": disease_name,
            "inheritance_score": score,
            "inheritance_distribution": dist,
            "inheritance_expected": expected,
            "inheritance_variance": _distribution_variance(dist, values),
            "justification": justification,
            "pmids": pmids,
        })
    
    return results


def _distributions_to_scores_onset_severity(distributions: dict) -> dict:
    """
    Convert onset/severity distributions to score map compatible with existing pipeline.
    
    Args:
        distributions: Dict mapping disease names to dicts with onset_distribution and severity_distribution
            e.g., {"Huntington": {"onset": {...}, "severity": {...}, "justification": "..."}}
    
    Returns:
        Dict mapping disease names to (onset_score, severity_score) tuples
    """
    result_map = {}
    onset_values = [1, 2, 3, 4, 5, 6, 7]  # Prenatal through late-onset
    severity_values = [1, 2, 3, 4, 5]  # Lethal through very mild
    
    for disease_name, dist_data in distributions.items():
        onset_dist = dist_data.get("onset", {})
        severity_dist = dist_data.get("severity", {})
        justification = dist_data.get("justification")
        
        onset_expected = _distribution_to_expected(onset_dist, onset_values)
        severity_expected = _distribution_to_expected(severity_dist, severity_values)
        
        onset_score = round(onset_expected) if onset_expected is not None else None
        severity_score = round(severity_expected) if severity_expected is not None else None
        
        result_map[disease_name] = {
            "onset_score": onset_score,
            "severity_score": severity_score,
            "onset_distribution": onset_dist,
            "severity_distribution": severity_dist,
            "onset_expected": onset_expected,
            "severity_expected": severity_expected,
            "onset_variance": _distribution_variance(onset_dist, onset_values),
            "severity_variance": _distribution_variance(severity_dist, severity_values),
            "justification": justification,
        }
    
    return result_map


def _apply_distributions_to_diseases(records, penetrance_dists, inheritance_dists, onset_severity_dists):
    """
    Apply probability distributions to disease records.
    
    Args:
        records: List of disease record dicts
        penetrance_dists: Dict from disease name to penetrance data (new format: {"distribution": {...}, "justification": ..., "pmids": [...]})
        inheritance_dists: Dict from disease name to inheritance data (new format: {"distribution": {...}, "justification": ..., "pmids": [...]})
        onset_severity_dists: Dict from disease name to onset/severity data ({"onset": {...}, "severity": {...}, "justification": ...})
    """
    for rec in records:
        name = rec.get("name")
        if not name:
            continue
        
        # Apply penetrance distribution
        if name in penetrance_dists:
            pen_data = penetrance_dists[name]
            # Handle both new format (with "distribution" key) and legacy format
            if "distribution" in pen_data:
                pen_dist = pen_data["distribution"]
                rec["penetrance_justification"] = pen_data.get("justification")
                rec["penetrance_pmids"] = pen_data.get("pmids", [])
            else:
                pen_dist = pen_data
            rec["penetrance_distribution"] = pen_dist
            values = [1, 2, 3, 4]
            rec["penetrance_expected"] = _distribution_to_expected(pen_dist, values)
            rec["penetrance_variance"] = _distribution_variance(pen_dist, values)
        
        # Apply inheritance distribution
        if name in inheritance_dists:
            inh_data = inheritance_dists[name]
            # Handle both new format (with "distribution" key) and legacy format
            if "distribution" in inh_data:
                inh_dist = inh_data["distribution"]
                rec["inheritance_justification"] = inh_data.get("justification")
                rec["inheritance_pmids"] = inh_data.get("pmids", [])
            else:
                inh_dist = inh_data
            rec["inheritance_distribution"] = inh_dist
            values = [1, 2, 3, 4, 5, 6]
            rec["inheritance_expected"] = _distribution_to_expected(inh_dist, values)
            rec["inheritance_variance"] = _distribution_variance(inh_dist, values)
        
        # Apply onset/severity distributions
        if name in onset_severity_dists:
            dist_data = onset_severity_dists[name]
            # Get distributions (use both key formats for compatibility)
            onset_dist = dist_data.get("onset") or dist_data.get("onset_distribution", {})
            severity_dist = dist_data.get("severity") or dist_data.get("severity_distribution", {})
            rec["onset_distribution"] = onset_dist
            rec["severity_distribution"] = severity_dist
            # Apply justification for onset/severity
            rec["severity_justification"] = dist_data.get("justification")
            # Compute expected values and variance if not already present
            onset_values = [1, 2, 3, 4, 5, 6, 7]
            severity_values = [1, 2, 3, 4, 5]
            rec["onset_expected"] = dist_data.get("onset_expected") or _distribution_to_expected(onset_dist, onset_values)
            rec["severity_expected"] = dist_data.get("severity_expected") or _distribution_to_expected(severity_dist, severity_values)
            rec["onset_variance"] = dist_data.get("onset_variance") or _distribution_variance(onset_dist, onset_values)
            rec["severity_variance"] = dist_data.get("severity_variance") or _distribution_variance(severity_dist, severity_values)


# ---------------------------------------------------------------------------
# Monte Carlo simulation for probabilistic algorithmic level
# ---------------------------------------------------------------------------

MONTE_CARLO_SEED = 42
MONTE_CARLO_SAMPLES = 10000


def _sample_from_distribution(distribution: dict, values: list) -> int:
    """
    Sample a single value from a probability distribution.
    
    Args:
        distribution: Dict mapping category names to probabilities (0-100)
        values: List of numeric values corresponding to each category
    
    Returns:
        Sampled value from the distribution
    """
    if not distribution:
        return None
    
    probs = list(distribution.values())
    total = sum(probs)
    if total == 0:
        return None
    
    # Normalize probabilities
    probs = [p / total for p in probs]
    
    # Sample using random.choices
    sampled_value = random.choices(values, weights=probs, k=1)[0]
    return sampled_value


def compute_algorithmic_level_distribution(record: dict, n_samples: int = MONTE_CARLO_SAMPLES, seed: int = MONTE_CARLO_SEED) -> dict:
    """
    Compute the distribution of algorithmic levels using Monte Carlo simulation.
    
    Samples from the probability distributions of penetrance, inheritance, onset,
    and severity to estimate the distribution of the final algorithmic level.
    
    Args:
        record: Disease record with distributions
        n_samples: Number of Monte Carlo samples
        seed: Random seed for reproducibility
    
    Returns:
        Dict with:
            - "level_distribution": Dict mapping levels 1-7 to probabilities (0-100)
            - "expected_level": Expected (mean) level
            - "level_variance": Variance of the level distribution
            - "samples": Raw sample counts for each level (for recalculation)
    """
    random.seed(seed)
    
    # Fixed values from the record (not probabilistic)
    a1 = record.get("association_score")
    
    # Get distributions
    pen_dist = record.get("penetrance_distribution", {})
    inh_dist = record.get("inheritance_distribution", {})
    onset_dist = record.get("onset_distribution", {})
    severity_dist = record.get("severity_distribution", {})
    
    # Values for each distribution
    pen_values = [1, 2, 3, 4]  # mendelian, high, moderate, complex
    inh_values = [1, 2, 3, 4, 5, 6]  # dominant through recessive
    onset_values = [1, 2, 3, 4, 5, 6, 7]  # prenatal through late
    severity_values = [1, 2, 3, 4, 5]  # lethal through very mild
    
    # If no distributions available, return single deterministic level
    if not any([pen_dist, inh_dist, onset_dist, severity_dist]):
        level = _compute_algorithmic_level(record)
        return {
            "level_distribution": {str(level): 100},
            "expected_level": float(level),
            "level_variance": 0.0,
            "samples": {str(level): n_samples},
        }
    
    # Monte Carlo sampling
    level_counts = {str(i): 0 for i in range(1, 8)}
    
    for _ in range(n_samples):
        # Sample each variable
        a2 = _sample_from_distribution(pen_dist, pen_values) if pen_dist else record.get("penetrance_level")
        a3 = _sample_from_distribution(inh_dist, inh_values) if inh_dist else record.get("inheritance_score")
        onset = _sample_from_distribution(onset_dist, onset_values) if onset_dist else record.get("disease_onset_score")
        severity = _sample_from_distribution(severity_dist, severity_values) if severity_dist else record.get("severity_score")
        
        # Compute level for this sample
        level = _compute_algorithmic_level_from_scores(a1, a2, a3, onset, severity)
        level_counts[str(level)] += 1
    
    # Convert counts to probabilities
    level_distribution = {k: round(v * 100 / n_samples, 2) for k, v in level_counts.items()}
    
    # Compute expected value and variance
    expected = sum(int(k) * v / 100 for k, v in level_distribution.items())
    variance = sum((int(k) - expected) ** 2 * v / 100 for k, v in level_distribution.items())
    
    return {
        "level_distribution": level_distribution,
        "expected_level": round(expected, 4),
        "level_variance": round(variance, 4),
        "samples": level_counts,
    }


def _compute_algorithmic_level_from_scores(a1, a2, a3, onset, severity):
    """
    Compute algorithmic level from raw scores (used by Monte Carlo).
    
    This is a standalone version of _compute_algorithmic_level that takes
    scores directly rather than extracting from a record dict.
    """
    def within(value, allowed):
        if value is None:
            return False
        return allowed(value)
    
    # V4: 7 levels, with level 1 very strict (a2=1)
    checks = [
        (
            1,  # Very strict - textbook Mendelian with high penetrance
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v == 1)
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 3)
            and within(severity, lambda v: v == 1),
        ),
        (
            2,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v in (1, 2))
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 4)
            and within(severity, lambda v: v == 1),
        ),
        (
            3,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v in (1, 2))
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 5)
            and within(severity, lambda v: v in (1, 2)),
        ),
        (
            4,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: 1 <= v <= 3)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 6)
            and within(severity, lambda v: v in (1, 2)),
        ),
        (
            5,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: 1 <= v <= 3)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 7)
            and within(severity, lambda v: 1 <= v <= 3),
        ),
        (
            6,
            lambda: within(a1, lambda v: v in (1, 2))
            and within(a2, lambda v: 1 <= v <= 5)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 9)
            and within(severity, lambda v: 1 <= v <= 4),
        ),
    ]
    
    for level, predicate in checks:
        if predicate():
            return level
    return DEFAULT_ALGO_LEVEL


def compute_kappa_from_variance(variance: float, min_kappa: float = 2.0, max_kappa: float = 50.0) -> float:
    """
    Compute kappa parameter for Bayesian prior from the variance of the level distribution.
    
    Higher variance (more uncertainty) → lower kappa (weaker prior)
    Lower variance (more confidence) → higher kappa (stronger prior)
    
    The mapping uses an inverse relationship scaled to the expected variance range:
    - Variance near 0 → kappa near max_kappa (very confident LLM)
    - Variance near 3-4 (max possible ~9) → kappa near min_kappa (very uncertain LLM)
    
    Args:
        variance: Variance of the algorithmic level distribution
        min_kappa: Minimum kappa value (for very uncertain cases)
        max_kappa: Maximum kappa value (for very confident cases)
    
    Returns:
        kappa: Parameter for the Bayesian prior
    """
    # Maximum theoretical variance for discrete uniform on 1-7 is about 4
    # In practice, distributions are usually more concentrated
    max_variance = 4.0
    
    # Clamp variance to expected range
    variance = max(0.0, min(variance, max_variance))
    
    # Linear interpolation: high variance → low kappa, low variance → high kappa
    # kappa = max_kappa - (max_kappa - min_kappa) * (variance / max_variance)
    kappa = max_kappa - (max_kappa - min_kappa) * (variance / max_variance)
    
    return round(kappa, 2)


def _compute_monte_carlo_levels(records: list) -> None:
    """
    Compute Monte Carlo algorithmic level distributions for all disease records.
    
    Modifies records in place, adding:
        - level_distribution: Dict mapping levels 1-7 to probabilities
        - expected_level: Expected (mean) level from distribution
        - level_variance: Variance of the level distribution
        - level_samples: Raw sample counts for each level
        - kappa: Computed kappa from variance for Bayesian prior
    
    Args:
        records: List of disease record dicts with distributions applied
    """
    for rec in records:
        # Skip if protective or neutral
        if rec.get("association_is_protective") or rec.get("association_is_neutral"):
            continue
        
        # Check if record has any distributions
        has_distributions = any([
            rec.get("penetrance_distribution"),
            rec.get("inheritance_distribution"),
            rec.get("onset_distribution"),
            rec.get("severity_distribution"),
        ])
        
        if has_distributions:
            mc_result = compute_algorithmic_level_distribution(rec)
            rec["level_distribution"] = mc_result["level_distribution"]
            rec["expected_level"] = mc_result["expected_level"]
            rec["level_variance"] = mc_result["level_variance"]
            rec["level_samples"] = mc_result["samples"]
            rec["kappa"] = compute_kappa_from_variance(mc_result["level_variance"])


def _get_min_expected_level_with_metadata(records: list) -> dict:
    """
    Get the minimum expected_level from all disease records along with its metadata.
    
    In proba mode, this returns the decimal expected level from Monte Carlo,
    which replaces the integer algorithmic_level, along with variance and distribution.
    
    Args:
        records: List of disease record dicts
    
    Returns:
        Dict with keys: expected_level, level_variance, level_distribution, kappa
        All values are None if no expected_level found
    """
    min_level = None
    min_record = None
    
    for rec in records:
        # Skip if protective or neutral
        if rec.get("association_is_protective") or rec.get("association_is_neutral"):
            continue
        
        expected = rec.get("expected_level")
        if expected is not None:
            if min_level is None or expected < min_level:
                min_level = expected
                min_record = rec
    
    if min_record is not None:
        return {
            "expected_level": min_record.get("expected_level"),
            "level_variance": min_record.get("level_variance"),
            "level_distribution": min_record.get("level_distribution"),
            "kappa": min_record.get("kappa"),
        }
    
    return {
        "expected_level": None,
        "level_variance": None,
        "level_distribution": None,
        "kappa": None,
    }


def _compute_algorithmic_level(record):
    a1 = record.get("association_score")
    a2 = record.get("penetrance_level")
    a3 = record.get("inheritance_score")
    onset = record.get("disease_onset_score")
    severity = record.get("severity_score")

    def within(value, allowed):
        if value is None:
            return False
        return allowed(value)

    # V4: 7 levels, with level 1 very strict (a2=1)
    # | Level  | a1      | a2      | a3        | onset   | severity |
    # |--------|---------|---------|-----------|---------|----------|
    # | 1      | =1      | =1      | {1,2,4}   | [1-3]   | =1       |
    # | 2      | =1      | {1,2}   | {1,2,4}   | [1-4]   | =1       |
    # | 3      | =1      | {1,2}   | {1,2,4}   | [1-5]   | {1,2}    |
    # | 4      | =1      | [1-3]   | [1-8]     | [1-6]   | {1,2}    |
    # | 5      | =1      | [1-3]   | [1-8]     | [1-7]   | [1-3]    |
    # | 6      | {1,2}   | [1-5]   | [1-8]     | [1-9]   | [1-4]    |
    # | 7      | default |         |           |         |          |
    checks = [
        (
            1,  # Very strict - textbook Mendelian with high penetrance
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v == 1)
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 3)
            and within(severity, lambda v: v == 1),
        ),
        (
            2,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v in (1, 2))
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 4)
            and within(severity, lambda v: v == 1),
        ),
        (
            3,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v in (1, 2))
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 5)
            and within(severity, lambda v: v in (1, 2)),
        ),
        (
            4,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: 1 <= v <= 3)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 6)
            and within(severity, lambda v: v in (1, 2)),
        ),
        (
            5,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: 1 <= v <= 3)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 7)
            and within(severity, lambda v: 1 <= v <= 3),
        ),
        (
            6,
            lambda: within(a1, lambda v: v in (1, 2))
            and within(a2, lambda v: 1 <= v <= 5)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 9)
            and within(severity, lambda v: 1 <= v <= 4),
        ),
    ]

    for level, predicate in checks:
        if predicate():
            return level
    return DEFAULT_ALGO_LEVEL




def _run_algorithmic_summary_agent(gene_name, cfg, diseases):
    if not diseases:
        return None

    simple_mode = cfg.get("simple", False)
    prompt_template = load_prompt(_get_prompt_name("deep_analysis_algorithmic_summary_agent", simple_mode=simple_mode))
    if not prompt_template:
        return None

    pathogenic_records = [rec for rec in diseases if not rec.get("association_is_protective") and not rec.get("association_is_neutral")]
    definitive_records = [
        rec for rec in pathogenic_records
        if rec.get("association_score") == 1
    ]
    if definitive_records:
        selected_records = definitive_records
    else:
        selected_records = (pathogenic_records or diseases)[:3]

    disease_blocks = []
    for rec in selected_records:
        disease_blocks.append(_format_algorithmic_summary_block(rec))

    prompt = prompt_template.format(
        gene_name=gene_name,
        disease_blocks="\n\n".join(block for block in disease_blocks if block.strip()) or "No disease data available."
    )

    try:
        response = call_llm(
            prompt,
            cfg["model"],
            cfg.get("temperature", 0.0),
            agent_name="algorithmic_summary_agent",
        )
        return parse_algorithmic_summary_response(response)
    except Exception as exc:  # pragma: no cover
        print(f"[DeepAnalysis] Algorithmic summary agent failed: {exc}")
        return None


def _label_from_inheritance_score(score):
    if score is None:
        return None
    return INHERITANCE_SCORE_LABELS.get(score, "unknown-inheritance")


def _label_from_onset_score(score):
    if score is None:
        return None
    return ONSET_SCORE_LABELS.get(score, "Unknown onset")


def _format_algorithmic_summary_block(rec):
    name = rec.get("name") or "Unknown disease"
    association_label = ASSOCIATION_LABELS.get(rec.get("association_score"), "Unknown evidence")
    assoc_notes = rec.get("association_justification") or "No association notes provided."

    penetrance_text = PENETRANCE_LEVEL_LABELS.get(rec.get("penetrance_level"), "unknown penetrance pattern")
    if rec.get("penetrance_justification"):
        penetrance_text += f" ({rec['penetrance_justification']})"

    inheritance_label = rec.get("inheritance") or _label_from_inheritance_score(rec.get("inheritance_score")) or "unknown inheritance pattern"
    inheritance_label = inheritance_label.replace("_", " ").replace("-", " ")
    inheritance_label = inheritance_label.strip()
    if inheritance_label:
        inheritance_label = inheritance_label[0].upper() + inheritance_label[1:]

    inheritance_text = inheritance_label
    if rec.get("inheritance_justification"):
        inheritance_text += f" ({rec['inheritance_justification']})"

    onset_label = rec.get("disease_onset") or _label_from_onset_score(rec.get("disease_onset_score")) or "unknown onset"
    severity_text = SEVERITY_LEVEL_LABELS.get(rec.get("severity_score"), "severity not reported")

    mechanism_text = rec.get("mechanism") or "Mechanism uncertain"
    if rec.get("mechanism_justification"):
        mechanism_text += f" ({rec['mechanism_justification']})"

    pathogenicity_flag = "protective association" if rec.get("association_is_protective") else (
        "neutral association" if rec.get("association_is_neutral") else "pathogenic association"
    )

    lines = [
        f"Disease: {name}",
        f"Evidence: {association_label}. {assoc_notes}",
        f"Clinical profile: {severity_text}; onset typically {onset_label.lower() if isinstance(onset_label, str) else onset_label}; penetrance shows {penetrance_text}.",
        f"Inheritance: {inheritance_text}.",
        f"Mechanism: {mechanism_text}.",
        f"Pathogenicity: {pathogenicity_flag}.",
    ]
    return "\n".join(line.strip() for line in lines if line.strip())
