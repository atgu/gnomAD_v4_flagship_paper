"""Utility helpers for the deep analysis pipeline."""
from __future__ import annotations

import copy
from typing import Dict, Iterable, List, Sequence

from services.pubmed_service import (
    fetch_abstracts,
    format_articles_for_llm,
    search_pubmed_query,
)


def prepare_disease_names(a1_results: Sequence[Dict]) -> List[str]:
    names: List[str] = []
    seen = set()
    for entry in a1_results:
        name = (entry.get("disease") or "").strip()
        if not name or name in seen:
            continue
        seen.add(name)
        names.append(name)
    return names


def build_keyword_clause(keywords: Iterable[str]) -> str:
    quoted = [kw if kw.startswith('"') else f'"{kw}"' for kw in keywords]
    return "(" + " OR ".join(quoted) + ")"


def build_or_clause(items: Iterable[str]) -> str:
    cleaned = [item.strip() for item in items if item and item.strip()]
    if not cleaned:
        return '"unknown"'
    quoted = [f'"{term}"' for term in cleaned]
    if len(quoted) == 1:
        return quoted[0]
    return "(" + " OR ".join(quoted) + ")"


def format_disease_list(diseases: Iterable[str]) -> str:
    return "\n".join(f"- {d}" for d in diseases)


def format_association_summaries(a1_results: Sequence[Dict]) -> str:
    if not a1_results:
        return "- No association evidence available."
    lines = []
    for entry in a1_results:
        name = entry.get("disease") or "Unknown disease"
        score = entry.get("association_score")
        justification = (entry.get("justification") or "").strip()
        pmids = entry.get("pmids") or []
        pmid_text = ", ".join(str(p) for p in pmids[:5]) if pmids else "none"
        score_text = f"Level {score}" if score else "Level N/A"
        lines.append(
            f"- {name}: {score_text}; Evidence: {justification or 'No summary'}; PMIDs: {pmid_text}"
        )
    return "\n".join(lines)


def fetch_articles_for_query(query: str, max_results: int, use_abstracts: bool):
    articles = search_pubmed_query(query, max_results=max_results)
    if use_abstracts and articles:
        pmids = [a.get("pmid") for a in articles if a.get("pmid")]
        pmids = [pmid for pmid in pmids if pmid]
        if pmids:
            abstracts = fetch_abstracts(pmids)
            for article in articles:
                pmid = article.get("pmid")
                if pmid:
                    article["abstract"] = abstracts.get(pmid, article.get("abstract"))
    return articles


def articles_to_text(articles: Sequence[Dict], use_abstracts: bool, top_abstracts: int = None) -> str:
    if not articles:
        return "No supporting articles were found for this query."
    return format_articles_for_llm(articles, use_abstracts, top_abstracts=top_abstracts)


def combine_articles(primary, base_by_pmid: Dict[str, Dict], a1_results: Sequence[Dict]):
    merged = list(primary or [])
    existing = {art.get("pmid") for art in merged if art.get("pmid")}

    cited_pmids = {
        str(pmid).strip()
        for entry in a1_results or []
        for pmid in entry.get("pmids", [])
    }

    for pmid in cited_pmids:
        if pmid in existing:
            continue
        article = base_by_pmid.get(pmid)
        if article:
            merged.append(copy.deepcopy(article))
            existing.add(pmid)

    return merged
