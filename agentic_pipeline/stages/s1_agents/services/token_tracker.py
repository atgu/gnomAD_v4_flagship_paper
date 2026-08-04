"""Service for tracking LLM token usage and costs."""
from datetime import datetime
from typing import Dict, Optional
import os
import threading


# Pricing per million tokens (MTok). Standard PayGo prices.
# Flex PayGo applies a 0.5 multiplier; see FLEX_DISCOUNT below.
PRICING = {
    # Anthropic Claude models
    "claude-sonnet-4-5-20250514": {"input": 3.0, "output": 15.0},
    "claude-sonnet-4-5": {"input": 3.0, "output": 15.0},
    "sonnet-4-5": {"input": 3.0, "output": 15.0},
    "claude-haiku-4-5-20241022": {"input": 1.0, "output": 5.0},
    "claude-haiku-4-5": {"input": 1.0, "output": 5.0},
    "haiku-4-5": {"input": 1.0, "output": 5.0},
    # Google Gemini models (Vertex AI Standard PayGo pricing, global endpoint)
    "gemini-3.5-flash": {"input": 1.50, "output": 9.00},   # GA 2026-05-19, served on location='global'
    "gemini-3.1-flash-lite": {"input": 0.25, "output": 1.50},  # GA 2026-05-07
    "gemini-2.5-flash": {"input": 0.15, "output": 0.60},  # < 200k context
    "gemini-2.5-pro": {"input": 1.25, "output": 10.0},    # < 200k context
    "gemini-2.0-flash": {"input": 0.10, "output": 0.40},
    "gemini-1.5-flash": {"input": 0.075, "output": 0.30},
    "gemini-1.5-pro": {"input": 1.25, "output": 5.0},
    # Meta Llama models served as MaaS on Vertex AI (us-east5). Prices from
    # https://cloud.google.com/vertex-ai/generative-ai/pricing (2026-05).
    # We register both the short aliases used in CLI and the canonical
    # MaaS model ids returned by Vertex so get_pricing() works for either.
    "llama-4-scout": {"input": 0.25, "output": 0.70},
    "meta/llama-4-scout-17b-16e-instruct-maas": {"input": 0.25, "output": 0.70},
    "llama-4-maverick": {"input": 0.35, "output": 1.15},
    "meta/llama-4-maverick-17b-128e-instruct-maas": {"input": 0.35, "output": 1.15},
    "llama-3.3-70b": {"input": 0.72, "output": 0.72},
    "meta/llama-3.3-70b-instruct-maas": {"input": 0.72, "output": 0.72},
    # Google MedGemma 27B (self-hosted on a Vertex dedicated endpoint, billed by
    # GPU-hour, NOT per token). The value below is an *effective* per-Mtok rate
    # calibrated to the GPU-time cost of scoring genes, so the cost figures (which
    # extrapolate on a fixed ~32.6k tokens/gene budget) reflect the true cost to
    # score N genes. Measured on the 4-replica run: 959 genes in 98.6 min on
    # 4x g4-standard-48 (NVIDIA RTX PRO 6000 @ $4.4999/h each = $18/h) = $29.6
    # => $0.0308/gene. At the figure's 32613 tok/gene budget that is
    # $0.0308 / 0.032613 Mtok = $0.9455/Mtok (input==output, GPU time is
    # mechanism-agnostic). Deploy overhead (~25 min spin-up) is negligible at scale.
    "medgemma-27b-text-it": {"input": 0.9455, "output": 0.9455},
    "medgemma-27b": {"input": 0.9455, "output": 0.9455},
    "medgemma": {"input": 0.9455, "output": 0.9455},
    # Qwen3.6-27B: self-hosted on Vertex (g4-standard-48 + RTX PRO 6000 x1 per
    # replica @ $4.4999/GPU-h), billed by GPU-hour. Effective per-Mtok rate
    # calibrated to the REAL GPU time of the 1010-gene run (run_027): scoring
    # took 1.690 h on 4 replicas = 6.76 GPU-h = $30.42 for 32.66 Mtok
    # (input==output, GPU time is mechanism-agnostic) -> $0.9313/Mtok. Deploy
    # provisioning (~49 min) is excluded as negligible at scale.
    "qwen3.6-27b": {"input": 0.9313, "output": 0.9313},
    "qwen3-6-27b": {"input": 0.9313, "output": 0.9313},
    "qwen": {"input": 0.9313, "output": 0.9313},

    # DeepSeek V3.2 on Vertex AI Model-as-a-Service (fully managed, pay-per-token,
    # served on location='global'). Published Vertex MaaS rates.
    "deepseek-v3.2-maas": {"input": 0.56, "output": 1.68},
    "deepseek-ai/deepseek-v3.2-maas": {"input": 0.56, "output": 1.68},
    "deepseek-v3.2": {"input": 0.56, "output": 1.68},
    "deepseek": {"input": 0.56, "output": 1.68},
}

# Flex PayGo (Vertex AI) charges 50% of Standard PayGo prices on the same
# global endpoint. Trade-off: higher latency, best-effort SLA. See
# https://cloud.google.com/vertex-ai/generative-ai/docs/flex-paygo
FLEX_DISCOUNT = 0.5


def get_pricing(model: str) -> Optional[Dict[str, float]]:
    """Get pricing for a model. Returns None if model not in pricing list."""
    model_lower = model.lower()
    
    # Direct match
    if model_lower in PRICING:
        return PRICING[model_lower]
    
    # Partial match (longest key wins, to avoid e.g. "gemini-1.5-flash" matching
    # against "gemini-2.5-flash" which used to happen with the previous
    # iteration order). Sorting by descending length matches the most specific
    # SKU first.
    for key in sorted(PRICING.keys(), key=len, reverse=True):
        if key in model_lower:
            return PRICING[key]
    
    # Fallback based on model family - Anthropic. Strip a potential '@vertex'
    # suffix so 'claude-haiku-4-5@vertex' falls back to haiku pricing too.
    stripped = model_lower.replace("@vertex", "")
    if "sonnet" in stripped:
        return PRICING["sonnet-4-5"]
    if "haiku" in stripped:
        return PRICING["haiku-4-5"]

    # Fallback based on model family - Meta Llama on Vertex
    if "llama-4-scout" in model_lower or "llama-4-17b-16e" in model_lower:
        return PRICING["llama-4-scout"]
    if "llama-4-maverick" in model_lower or "llama-4-17b-128e" in model_lower:
        return PRICING["llama-4-maverick"]
    if "llama-3.3-70b" in model_lower or "llama-3-3-70b" in model_lower:
        return PRICING["llama-3.3-70b"]

    # Fallback based on model family - DeepSeek on Vertex MaaS
    if "deepseek" in model_lower:
        return PRICING["deepseek-v3.2-maas"]
    
    # Fallback based on model family - Google Gemini
    if "gemini-3.5-flash" in model_lower:
        return PRICING["gemini-3.5-flash"]
    if "gemini-3.1-flash-lite" in model_lower:
        return PRICING["gemini-3.1-flash-lite"]
    if "gemini-2.5-flash" in model_lower:
        return PRICING["gemini-2.5-flash"]
    if "gemini-2.5-pro" in model_lower:
        return PRICING["gemini-2.5-pro"]
    if "gemini-2.0" in model_lower:
        return PRICING["gemini-2.0-flash"]
    if "gemini-1.5-flash" in model_lower:
        return PRICING["gemini-1.5-flash"]
    if "gemini-1.5-pro" in model_lower:
        return PRICING["gemini-1.5-pro"]
    if "gemini" in model_lower:
        return PRICING["gemini-2.5-flash"]  # Default Gemini pricing
    
    return None


class TokenTracker:
    """
    Thread-safe tracker for LLM token usage across agents.
    
    Usage:
        tracker = TokenTracker()
        tracker.record("disease_agent", "claude-haiku-4-5", 1000, 500)
        tracker.record("penetrance_agent", "claude-haiku-4-5", 800, 400)
        tracker.write_report("/path/to/run/token_usage.txt")
    """
    
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self):
        self._data: Dict[str, Dict] = {}
        self._model: Optional[str] = None
        self._run_start: Optional[datetime] = None
        self._num_papers: Optional[int] = None
        self._use_abstracts: Optional[bool] = None
        self._data_lock = threading.Lock()
    
    @classmethod
    def get_instance(cls) -> "TokenTracker":
        """Get or create the singleton instance."""
        with cls._lock:
            if cls._instance is None:
                cls._instance = cls()
            return cls._instance
    
    @classmethod
    def reset_instance(cls):
        """Reset the singleton instance (for testing or new runs)."""
        with cls._lock:
            cls._instance = None
    
    def reset(self):
        """Reset all tracking data for a new run."""
        with self._data_lock:
            self._data = {}
            self._model = None
            self._run_start = datetime.now()
            self._num_papers = None
            self._use_abstracts = None
    
    def set_run_config(self, model: str, num_papers: int, use_abstracts: bool):
        """
        Set the run configuration for the report header.
        
        Args:
            model: Model name (e.g., "claude-haiku-4-5")
            num_papers: Number of papers used
            use_abstracts: Whether abstracts are used
        """
        with self._data_lock:
            self._model = model
            self._num_papers = num_papers
            self._use_abstracts = use_abstracts
    
    def record(
        self,
        agent_name: str,
        model: str,
        input_tokens: int,
        output_tokens: int,
        flex: bool = False,
    ):
        """
        Record token usage for an agent.
        
        Args:
            agent_name: Name of the agent (e.g., "disease_agent", "penetrance_agent")
            model: Model used (e.g., "claude-haiku-4-5")
            input_tokens: Number of input tokens
            output_tokens: Number of output tokens
            flex: If True, tokens were billed at Vertex AI Flex PayGo rates
                  (-50% vs Standard PayGo). The source of truth for this flag
                  is Vertex's response field `traffic_type == "ON_DEMAND_FLEX"`.
        """
        with self._data_lock:
            if self._run_start is None:
                self._run_start = datetime.now()
            
            if self._model is None:
                self._model = model
            
            if agent_name not in self._data:
                self._data[agent_name] = {
                    "model": model,
                    "input_tokens": 0,
                    "output_tokens": 0,
                    "call_count": 0,
                    "flex_input_tokens": 0,
                    "flex_output_tokens": 0,
                    "flex_call_count": 0,
                }
            
            self._data[agent_name]["input_tokens"] += input_tokens
            self._data[agent_name]["output_tokens"] += output_tokens
            self._data[agent_name]["call_count"] += 1
            if flex:
                self._data[agent_name]["flex_input_tokens"] += input_tokens
                self._data[agent_name]["flex_output_tokens"] += output_tokens
                self._data[agent_name]["flex_call_count"] += 1
    
    def get_totals(self) -> Dict:
        """Get total token counts across all agents."""
        with self._data_lock:
            total_input = sum(d["input_tokens"] for d in self._data.values())
            total_output = sum(d["output_tokens"] for d in self._data.values())
            return {
                "input_tokens": total_input,
                "output_tokens": total_output,
                "total_tokens": total_input + total_output,
            }

    def get_snapshot(self) -> Dict:
        """
        Return a JSON-serializable snapshot of the tracker state.

        Designed to be returned from a worker process (ProcessPoolExecutor) and
        merged back into the main-process tracker via merge_snapshot(). The
        snapshot is a shallow copy of the per-agent dict; integers and strings
        are pickle/JSON safe.
        """
        with self._data_lock:
            return {
                "data": {
                    agent: {
                        "model": d["model"],
                        "input_tokens": d["input_tokens"],
                        "output_tokens": d["output_tokens"],
                        "call_count": d["call_count"],
                        "flex_input_tokens": d.get("flex_input_tokens", 0),
                        "flex_output_tokens": d.get("flex_output_tokens", 0),
                        "flex_call_count": d.get("flex_call_count", 0),
                    }
                    for agent, d in self._data.items()
                }
            }

    def merge_snapshot(self, snapshot: Optional[Dict]):
        """
        Add the counts from another tracker's snapshot into this one.

        Used by the main process to aggregate per-gene token usage produced in
        worker processes. Safe to call with None or an empty dict (no-op).
        Per-agent entries are summed; the model recorded in this tracker is
        kept (the snapshot's model is used only when creating a new entry).
        """
        if not snapshot:
            return
        data = snapshot.get("data") or {}
        if not data:
            return
        with self._data_lock:
            if self._run_start is None:
                self._run_start = datetime.now()
            for agent_name, entry in data.items():
                if agent_name not in self._data:
                    self._data[agent_name] = {
                        "model": entry.get("model", self._model or "unknown"),
                        "input_tokens": 0,
                        "output_tokens": 0,
                        "call_count": 0,
                        "flex_input_tokens": 0,
                        "flex_output_tokens": 0,
                        "flex_call_count": 0,
                    }
                self._data[agent_name]["input_tokens"] += int(entry.get("input_tokens", 0))
                self._data[agent_name]["output_tokens"] += int(entry.get("output_tokens", 0))
                self._data[agent_name]["call_count"] += int(entry.get("call_count", 0))
                self._data[agent_name]["flex_input_tokens"] += int(entry.get("flex_input_tokens", 0))
                self._data[agent_name]["flex_output_tokens"] += int(entry.get("flex_output_tokens", 0))
                self._data[agent_name]["flex_call_count"] += int(entry.get("flex_call_count", 0))
    
    def calculate_cost(
        self,
        model: str,
        input_tokens: int,
        output_tokens: int,
        flex_input_tokens: int = 0,
        flex_output_tokens: int = 0,
    ) -> Optional[Dict]:
        """
        Calculate cost for given token counts.
        
        Flex tokens (those routed through Vertex AI Flex PayGo) are billed at
        FLEX_DISCOUNT * Standard price. `flex_input_tokens` and
        `flex_output_tokens` are assumed to be SUBSETS of `input_tokens` and
        `output_tokens` respectively (not added on top).
        
        Returns:
            Dict with input_cost, output_cost, total_cost, and the same fields
            for the flex breakdown. None if model not priced.
        """
        pricing = get_pricing(model)
        if pricing is None:
            return None
        
        std_input = max(0, input_tokens - flex_input_tokens)
        std_output = max(0, output_tokens - flex_output_tokens)
        
        std_input_cost = (std_input / 1_000_000) * pricing["input"]
        std_output_cost = (std_output / 1_000_000) * pricing["output"]
        flex_input_cost = (flex_input_tokens / 1_000_000) * pricing["input"] * FLEX_DISCOUNT
        flex_output_cost = (flex_output_tokens / 1_000_000) * pricing["output"] * FLEX_DISCOUNT
        
        input_cost = std_input_cost + flex_input_cost
        output_cost = std_output_cost + flex_output_cost
        
        return {
            "input_cost": input_cost,
            "output_cost": output_cost,
            "total_cost": input_cost + output_cost,
            "std_input_cost": std_input_cost,
            "std_output_cost": std_output_cost,
            "flex_input_cost": flex_input_cost,
            "flex_output_cost": flex_output_cost,
        }
    
    def generate_report(self, run_name: str = "unknown") -> str:
        """
        Generate a text report of token usage.
        
        Args:
            run_name: Name of the run for the report header
            
        Returns:
            Formatted text report
        """
        with self._data_lock:
            abstracts_str = "Yes" if self._use_abstracts else "No"
            lines = [
                "=" * 50,
                "TOKEN USAGE REPORT",
                f"Run: {run_name}",
                f"Date: {self._run_start.strftime('%Y-%m-%d %H:%M:%S') if self._run_start else 'N/A'}",
                f"Model: {self._model or 'N/A'}",
                f"Num papers: {self._num_papers or 'N/A'}",
                f"Use abstracts: {abstracts_str if self._use_abstracts is not None else 'N/A'}",
                "=" * 50,
                "",
                "DETAIL PAR AGENT:",
                "-" * 50,
            ]
            
            total_input = 0
            total_output = 0
            total_flex_input = 0
            total_flex_output = 0
            total_flex_calls = 0
            total_calls = 0
            total_cost = 0.0
            has_pricing = self._model and get_pricing(self._model) is not None
            
            # Sort agents for consistent output
            for agent_name in sorted(self._data.keys()):
                data = self._data[agent_name]
                input_tokens = data["input_tokens"]
                output_tokens = data["output_tokens"]
                call_count = data["call_count"]
                flex_input = data.get("flex_input_tokens", 0)
                flex_output = data.get("flex_output_tokens", 0)
                flex_calls = data.get("flex_call_count", 0)
                model = data["model"]
                
                total_input += input_tokens
                total_output += output_tokens
                total_flex_input += flex_input
                total_flex_output += flex_output
                total_flex_calls += flex_calls
                total_calls += call_count
                
                lines.append(f"{agent_name}:")
                if flex_calls > 0:
                    lines.append(f"  Calls: {call_count} (Flex: {flex_calls})")
                else:
                    lines.append(f"  Calls: {call_count}")
                lines.append(f"  Tokens input:  {input_tokens:,}")
                lines.append(f"  Tokens output: {output_tokens:,}")
                if flex_input > 0 or flex_output > 0:
                    lines.append(f"  Of which Flex: {flex_input:,} in / {flex_output:,} out")
                
                cost_info = self.calculate_cost(model, input_tokens, output_tokens, flex_input, flex_output)
                if cost_info:
                    agent_cost = cost_info["total_cost"]
                    total_cost += agent_cost
                    lines.append(
                        f"  Cost: ${cost_info['input_cost']:.4f} (input) + "
                        f"${cost_info['output_cost']:.4f} (output) = ${agent_cost:.4f}"
                    )
                else:
                    lines.append(f"  Cost: N/A (model not in pricing list)")
                
                lines.append("")
            
            # Totals section
            lines.extend([
                "-" * 50,
                "TOTALS:",
                "-" * 50,
                f"Total calls:         {total_calls:,}",
                f"Total tokens input:  {total_input:,}",
                f"Total tokens output: {total_output:,}",
                f"Total tokens:        {total_input + total_output:,}",
            ])
            
            if total_flex_calls > 0:
                if total_calls > 0:
                    flex_pct = 100.0 * total_flex_calls / total_calls
                else:
                    flex_pct = 0.0
                lines.extend([
                    "",
                    f"Flex PayGo calls:    {total_flex_calls:,} / {total_calls:,} ({flex_pct:.1f}%)",
                    f"Flex tokens input:   {total_flex_input:,}",
                    f"Flex tokens output:  {total_flex_output:,}",
                    f"(Flex tokens billed at {int(FLEX_DISCOUNT * 100)}% of Standard PayGo rates)",
                ])
            
            lines.append("")
            
            if has_pricing:
                lines.append(f"TOTAL COST: ${total_cost:.4f}")
            else:
                lines.append("TOTAL COST: N/A (model not in pricing list)")
            
            lines.append("=" * 50)
            
            return "\n".join(lines)
    
    def write_report(self, filepath: str, run_name: str = "unknown"):
        """
        Write the token usage report to a file.
        
        Args:
            filepath: Path to the output file
            run_name: Name of the run for the report header
        """
        report = self.generate_report(run_name)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"[TokenTracker] Report written to {filepath}")
    
    def append_to_report(self, filepath: str, gene_name: str):
        """
        Append current tracking data to an existing report file.
        Used for incremental updates (e.g., after each gene in webapp).
        
        Args:
            filepath: Path to the output file
            gene_name: Name of the gene just processed
        """
        with self._data_lock:
            if not self._data:
                return
            
            # Check if file exists - if not, add header
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            file_exists = os.path.exists(filepath) and os.path.getsize(filepath) > 0
            
            lines = []
            
            if not file_exists:
                # Add header with model and config info
                abstracts_str = "Yes" if self._use_abstracts else "No"
                lines.extend([
                    "=" * 50,
                    "TOKEN USAGE REPORT",
                    f"Model: {self._model or 'N/A'}",
                    f"Num papers: {self._num_papers or 'N/A'}",
                    f"Use abstracts: {abstracts_str if self._use_abstracts is not None else 'N/A'}",
                    f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                    "=" * 50,
                ])
            
            lines.extend([
                "",
                f"--- Gene: {gene_name} ---",
                f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            ])
            
            total_input = 0
            total_output = 0
            total_cost = 0.0
            
            for agent_name in sorted(self._data.keys()):
                data = self._data[agent_name]
                input_tokens = data["input_tokens"]
                output_tokens = data["output_tokens"]
                flex_input = data.get("flex_input_tokens", 0)
                flex_output = data.get("flex_output_tokens", 0)
                model = data["model"]
                
                total_input += input_tokens
                total_output += output_tokens
                
                cost_info = self.calculate_cost(model, input_tokens, output_tokens, flex_input, flex_output)
                cost_str = f"${cost_info['total_cost']:.4f}" if cost_info else "N/A"
                if cost_info:
                    total_cost += cost_info['total_cost']
                
                flex_suffix = f" (Flex {flex_input:,}/{flex_output:,})" if (flex_input or flex_output) else ""
                lines.append(
                    f"  {agent_name}: {input_tokens:,} in / {output_tokens:,} out = {cost_str}{flex_suffix}"
                )
            
            # Add gene total
            lines.append(f"  ----")
            lines.append(f"  TOTAL: {total_input:,} in / {total_output:,} out = ${total_cost:.4f}")
            lines.append("")
            
            # Update cumulative totals file
            self._update_cumulative_totals(filepath, total_input, total_output, total_cost)
            
            with open(filepath, "a", encoding="utf-8") as f:
                f.write("\n".join(lines))
    
    def _update_cumulative_totals(self, filepath: str, input_tokens: int, output_tokens: int, cost: float):
        """
        Update cumulative totals in a separate file.
        """
        totals_path = filepath.replace(".txt", "_totals.txt")
        
        # Read existing totals
        cumul_input = 0
        cumul_output = 0
        cumul_cost = 0.0
        gene_count = 0
        
        if os.path.exists(totals_path):
            try:
                with open(totals_path, "r", encoding="utf-8") as f:
                    for line in f:
                        if line.startswith("cumul_input="):
                            cumul_input = int(line.split("=")[1].strip())
                        elif line.startswith("cumul_output="):
                            cumul_output = int(line.split("=")[1].strip())
                        elif line.startswith("cumul_cost="):
                            cumul_cost = float(line.split("=")[1].strip())
                        elif line.startswith("gene_count="):
                            gene_count = int(line.split("=")[1].strip())
            except Exception:
                pass
        
        # Update totals
        cumul_input += input_tokens
        cumul_output += output_tokens
        cumul_cost += cost
        gene_count += 1
        
        # Write updated totals
        os.makedirs(os.path.dirname(totals_path), exist_ok=True)
        with open(totals_path, "w", encoding="utf-8") as f:
            f.write(f"cumul_input={cumul_input}\n")
            f.write(f"cumul_output={cumul_output}\n")
            f.write(f"cumul_cost={cumul_cost}\n")
            f.write(f"gene_count={gene_count}\n")
            f.write("\n")
            f.write("=" * 50 + "\n")
            f.write("CUMULATIVE TOTALS\n")
            f.write("=" * 50 + "\n")
            f.write(f"Genes processed: {gene_count}\n")
            f.write(f"Total tokens input:  {cumul_input:,}\n")
            f.write(f"Total tokens output: {cumul_output:,}\n")
            f.write(f"Total tokens:        {cumul_input + cumul_output:,}\n")
            f.write(f"TOTAL COST: ${cumul_cost:.4f}\n")
            f.write("=" * 50 + "\n")


# Convenience function for global access
def get_tracker() -> TokenTracker:
    """Get the global TokenTracker instance."""
    return TokenTracker.get_instance()

