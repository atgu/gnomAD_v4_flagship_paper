#!/usr/bin/env python3
"""
Recompute token usage and cost for an existing run by walking its per-gene
JSON files in ``agent_runs/<run_id>/results/``.

This is a retroactive fix for the multiprocessing bug that left
``token_usage.txt`` at zero whenever ``--n_core > 1`` was used (each worker
had its own singleton ``TokenTracker``, so the main-process report received
nothing).

Per-gene JSONs are unaffected by the bug: every agent call records its
``input_tokens`` and ``output_tokens`` directly in the gene's JSON payload.
This script walks each JSON, collects every dict that contains both
``input_tokens`` and ``output_tokens``, groups by ``agent_name``/``model``
and writes a report identical in shape to ``token_usage.txt`` (plus a
``token_usage_recomputed.json`` with the structured breakdown).

Usage:
    python -m utils.recompute_token_usage run_016
    python -m utils.recompute_token_usage run_016 --output custom_report.txt
    python -m utils.recompute_token_usage run_016 --all-runs     # recompute every run_* folder
    python -m utils.recompute_token_usage run_016 --quiet
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Tuple

# Make the project root importable so we can reuse the shared pricing table
THIS_DIR = Path(__file__).resolve().parent
APP_DIR = THIS_DIR.parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

from services.token_tracker import get_pricing  # noqa: E402

AGENT_RUNS_DIR = APP_DIR / "agent_runs"


def _iter_token_records(node) -> Iterable[Tuple[str, str, int, int]]:
    """
    Recursively yield (agent_name, model, input_tokens, output_tokens) for
    every dict in ``node`` that carries both token counts.

    The structure of a gene JSON nests token records under arbitrary keys
    (``raw.a1``, ``raw.a4``, ``deep_analysis.algorithmic_summary``, etc.),
    so we walk the whole tree instead of hard-coding paths.
    """
    if isinstance(node, dict):
        has_in = "input_tokens" in node
        has_out = "output_tokens" in node
        if has_in and has_out:
            agent = node.get("agent_name") or "unknown_agent"
            model = node.get("model") or "unknown"
            try:
                in_tok = int(node.get("input_tokens") or 0)
                out_tok = int(node.get("output_tokens") or 0)
            except (TypeError, ValueError):
                in_tok = out_tok = 0
            if in_tok or out_tok:
                yield agent, model, in_tok, out_tok
        for value in node.values():
            yield from _iter_token_records(value)
    elif isinstance(node, list):
        for item in node:
            yield from _iter_token_records(item)


def aggregate_run(results_dir: Path) -> Dict:
    """
    Walk every ``*.json`` in ``results_dir`` and return a structured summary:

        {
            "per_agent": {agent_name: {"model": str, "input_tokens": int,
                                       "output_tokens": int, "call_count": int}},
            "per_model": {model: {"input_tokens": int, "output_tokens": int,
                                  "call_count": int, "cost": float|None}},
            "gene_count": int,
            "total_input": int,
            "total_output": int,
            "total_cost": float,
            "has_full_pricing": bool,
        }
    """
    per_agent: Dict[str, Dict] = defaultdict(lambda: {
        "model": "unknown",
        "input_tokens": 0,
        "output_tokens": 0,
        "call_count": 0,
    })
    per_model: Dict[str, Dict] = defaultdict(lambda: {
        "input_tokens": 0,
        "output_tokens": 0,
        "call_count": 0,
    })

    gene_count = 0
    files = sorted(results_dir.glob("*.json"))
    for json_path in files:
        gene_count += 1
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as exc:
            print(f"[WARN] Cannot read {json_path.name}: {exc}", file=sys.stderr)
            continue

        for agent, model, in_tok, out_tok in _iter_token_records(data):
            entry = per_agent[agent]
            entry["model"] = model  # last write wins (all calls of an agent use the same model)
            entry["input_tokens"] += in_tok
            entry["output_tokens"] += out_tok
            entry["call_count"] += 1

            mentry = per_model[model]
            mentry["input_tokens"] += in_tok
            mentry["output_tokens"] += out_tok
            mentry["call_count"] += 1

    total_input = sum(e["input_tokens"] for e in per_agent.values())
    total_output = sum(e["output_tokens"] for e in per_agent.values())
    total_cost = 0.0
    has_full_pricing = True
    for model, entry in per_model.items():
        pricing = get_pricing(model)
        if pricing is None:
            entry["cost"] = None
            has_full_pricing = False
            continue
        cost_in = entry["input_tokens"] / 1_000_000 * pricing["input"]
        cost_out = entry["output_tokens"] / 1_000_000 * pricing["output"]
        entry["cost"] = cost_in + cost_out
        entry["input_cost"] = cost_in
        entry["output_cost"] = cost_out
        total_cost += entry["cost"]

    return {
        "per_agent": dict(per_agent),
        "per_model": dict(per_model),
        "gene_count": gene_count,
        "total_input": total_input,
        "total_output": total_output,
        "total_cost": total_cost,
        "has_full_pricing": has_full_pricing,
    }


def format_report(run_name: str, summary: Dict) -> str:
    sep = "=" * 60
    sub = "-" * 60
    agents_seen = set(summary["per_agent"].keys())
    expected_upstream = {"disease_agent", "algorithmic_summary_agent"}
    missing_upstream = expected_upstream - agents_seen

    lines = [
        sep,
        "TOKEN USAGE REPORT (RECOMPUTED FROM PER-GENE JSONS)",
        f"Run: {run_name}",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Genes scanned: {summary['gene_count']}",
        sep,
        "",
    ]

    if missing_upstream:
        lines.extend([
            "NOTE: Only agents whose token counts are persisted in the per-gene",
            "      JSON payloads can be recovered retroactively. For runs created",
            "      before the multiprocessing token-tracker fix, the following",
            "      upstream agents are NOT included in this report (their tokens",
            "      were only recorded in worker-local trackers that died with",
            "      the worker process):",
            f"        - {', '.join(sorted(missing_upstream))}",
            "      As a result the cost reported below is a LOWER BOUND.",
            "      Rough estimate for the missing agents (Claude Haiku 4.5,",
            "      based on run_001 single-process baseline):",
            "        - disease_agent: ~$0.035 per gene with diseases found",
            "        - algorithmic_summary_agent: ~$0.002 per gene",
            "",
        ])

    lines.extend([
        "DETAIL PAR AGENT:",
        sub,
    ])

    for agent_name in sorted(summary["per_agent"].keys()):
        data = summary["per_agent"][agent_name]
        model = data["model"]
        in_tok = data["input_tokens"]
        out_tok = data["output_tokens"]
        call_count = data["call_count"]

        pricing = get_pricing(model)
        if pricing:
            cost_in = in_tok / 1_000_000 * pricing["input"]
            cost_out = out_tok / 1_000_000 * pricing["output"]
            cost_str = (f"${cost_in:.4f} (input) + ${cost_out:.4f} (output) = "
                        f"${cost_in + cost_out:.4f}")
        else:
            cost_str = "N/A (model not in pricing list)"

        lines.append(f"{agent_name}:")
        lines.append(f"  Model: {model}")
        lines.append(f"  Calls: {call_count:,}")
        lines.append(f"  Tokens input:  {in_tok:,}")
        lines.append(f"  Tokens output: {out_tok:,}")
        lines.append(f"  Cost: {cost_str}")
        lines.append("")

    lines.extend([
        sub,
        "PER MODEL:",
        sub,
    ])
    for model in sorted(summary["per_model"].keys()):
        m = summary["per_model"][model]
        cost = m.get("cost")
        cost_str = f"${cost:.4f}" if cost is not None else "N/A"
        lines.append(f"{model}: {m['call_count']:,} calls / "
                     f"{m['input_tokens']:,} in / {m['output_tokens']:,} out = {cost_str}")

    lines.extend([
        "",
        sub,
        "TOTAUX:",
        sub,
        f"Total tokens input:  {summary['total_input']:,}",
        f"Total tokens output: {summary['total_output']:,}",
        f"Total tokens:        {summary['total_input'] + summary['total_output']:,}",
        "",
    ])

    if summary["has_full_pricing"]:
        lines.append(f"TOTAL COST: ${summary['total_cost']:.4f}")
    else:
        lines.append(f"TOTAL COST (partial): ${summary['total_cost']:.4f} "
                     f"(some models have no listed price)")
    lines.append(sep)
    return "\n".join(lines)


def recompute_single_run(run_id: str, output_filename: str | None = None, quiet: bool = False) -> Dict:
    run_path = AGENT_RUNS_DIR / run_id
    results_dir = run_path / "results"
    if not results_dir.is_dir():
        raise FileNotFoundError(f"No results directory for {run_id}: {results_dir}")

    summary = aggregate_run(results_dir)
    report = format_report(run_id, summary)

    txt_filename = output_filename or "token_usage_recomputed.txt"
    json_filename = (Path(txt_filename).stem + ".json")
    txt_path = run_path / txt_filename
    json_path = run_path / json_filename

    txt_path.write_text(report, encoding="utf-8")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({
            "run_id": run_id,
            "generated_at": datetime.now().isoformat(),
            "gene_count": summary["gene_count"],
            "total_input": summary["total_input"],
            "total_output": summary["total_output"],
            "total_cost": summary["total_cost"],
            "has_full_pricing": summary["has_full_pricing"],
            "per_model": summary["per_model"],
            "per_agent": summary["per_agent"],
        }, f, indent=2)

    if not quiet:
        print(report)
        print(f"\n[INFO] Text report written to {txt_path}")
        print(f"[INFO] JSON breakdown written to {json_path}")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Recompute token usage and cost for an existing run from per-gene JSONs.",
    )
    parser.add_argument("run", nargs="?", help="Run ID (e.g. run_016)")
    parser.add_argument("--all-runs", action="store_true",
                        help="Recompute every run_* folder under agent_runs/")
    parser.add_argument("--output", default=None,
                        help="Output filename (relative to run dir). Default: token_usage_recomputed.txt")
    parser.add_argument("--quiet", action="store_true",
                        help="Only write files, skip stdout report")
    args = parser.parse_args()

    if args.all_runs:
        runs = sorted([d.name for d in AGENT_RUNS_DIR.iterdir() if d.is_dir() and d.name.startswith("run_")])
        if not runs:
            print("No run_* folders found under", AGENT_RUNS_DIR, file=sys.stderr)
            sys.exit(1)
        grand_total = 0.0
        for run_id in runs:
            try:
                s = recompute_single_run(run_id, args.output, quiet=args.quiet)
                grand_total += s["total_cost"]
                if args.quiet:
                    print(f"  {run_id}: {s['gene_count']:,} genes, "
                          f"{s['total_input']:,} in / {s['total_output']:,} out, "
                          f"${s['total_cost']:.2f}")
            except FileNotFoundError as exc:
                print(f"  {run_id}: skipped ({exc})", file=sys.stderr)
        print(f"\nGRAND TOTAL across {len(runs)} runs: ${grand_total:.2f}")
    else:
        if not args.run:
            parser.error("provide a run ID (e.g. run_016) or use --all-runs")
        recompute_single_run(args.run, args.output, quiet=args.quiet)


if __name__ == "__main__":
    main()
