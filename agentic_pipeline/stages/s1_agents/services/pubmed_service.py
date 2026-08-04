"""PubMed API client for searching and fetching articles."""
import os
import re
import time
import requests
from xml.etree import ElementTree
from datetime import datetime
from config import BASE_URL_NCBI, NCBI_API_KEY as CONFIG_NCBI_KEY

# No key is embedded here. The upstream working copy carried one in source,
# which is why that key must be treated as compromised and rotated. NCBI
# tolerates anonymous use at a lower rate limit (3 requests/s instead of 10),
# so the agents still run without one, just more slowly.
NCBI_API_KEY = CONFIG_NCBI_KEY or os.environ.get("NCBI_API_KEY") or None
if not NCBI_API_KEY:
    print(
        "Warning: no NCBI key (NCBI_API_KEY). PubMed requests are throttled "
        "to 3/s instead of 10/s. See agentic_pipeline/.env.example.",
        flush=True,
    )
PUBMED_CALL_COOLDOWN = 0.05
VERBOSE_MODE = False  # Global flag to control verbose output

# Optional publication-date ceiling for ALL PubMed esearch queries. When set
# (e.g. "2019/12/31"), every search is restricted to articles published on or
# before this date via datetype=pdat. Default None = no restriction (historical
# behavior unchanged). Used to reconstruct the literature as it stood at a past
# date (e.g. a "2019" run to test prospective discovery potential).
MAX_PUBDATE = None
MIN_PUBDATE = "1000/01/01"


def _apply_pubdate_cap(esearch_params):
    """Restrict an esearch params dict to MAX_PUBDATE if one is configured.

    Only esearch needs this; esummary/efetch fetch by explicit PMID so they
    inherit the date restriction for free. Mutates and returns the dict.
    """
    if MAX_PUBDATE:
        esearch_params["datetype"] = "pdat"
        esearch_params["mindate"] = MIN_PUBDATE
        esearch_params["maxdate"] = MAX_PUBDATE
    return esearch_params


def _filter_articles_by_pubdate_cap(articles):
    """Deterministically drop any article whose displayed publication year is
    after the configured cap.

    The server-side datetype=pdat restriction occasionally lets through records
    whose electronic pub date is <= cap but whose journal-issue date (the value
    shown in `pubdate`) is later. For a clean, defensible "as of <year>" run we
    also enforce the ceiling on the visible year here.
    """
    if not MAX_PUBDATE:
        return articles
    cap_year = int(str(MAX_PUBDATE)[:4])
    kept = []
    for a in articles:
        m = re.search(r"(19|20)\d{2}", str(a.get("pubdate", "")))
        # Drop only when we can parse a year AND it exceeds the cap. Unparseable
        # dates are kept (server-side pdat already restricted them).
        if m and int(m.group(0)) > cap_year:
            continue
        kept.append(a)
    return kept

# Genes whose symbol is ambiguous and needs a more specific PubMed query.
# For those, "X gene"[Title/Abstract] is used instead of "X"[Title/Abstract]
# to cut the noise (common English words, medical abbreviations, disease names)
AMBIGUOUS_GENE_NAMES = {
    # Medical/scientific abbreviations
    "APC", "BAD", "BAX", "BID", "CA2", "CAMP", "EGF", "EGFR", "FOS", "GC", "HR",
    "JUN", "MTOR", "MYC", "NGF", "PIP", "SRC", "TNF", "VHL",
    # Disease and syndrome names
    "CAD", "DMD", "FAP", "SCD",
    # Common English words in scientific writing
    "APP", "ARC", "CAT", "IMPACT", "KIN", "KIT", "MAG", "MAX", "MET",
    "NODAL", "RAN", "REST", "SET", "SON", "TANK", "TUB", "VIM",
    # Two-character genes (far too ambiguous)
    "AR", "C2", "C3", "C5", "C6", "C7", "C9", "CP", "CS", "F2", "F3", "F5",
    "F7", "F8", "F9", "FH", "GK", "HP", "IK", "KL", "KY", "MB", "PC", "SI",
    "TF", "TG", "TH", "XG", "XK",
}


def get_pubmed_gene_term(gene_name: str) -> str:
    """
    Returns the appropriate PubMed search term for a gene.
    
    For ambiguous genes (common English words, medical abbreviations, disease names),
    uses "GENE gene"[Title/Abstract] to reduce noise.
    For other genes, uses "GENE"[Title/Abstract].
    
    Args:
        gene_name: The gene symbol (e.g., "SET", "BRCA1")
        
    Returns:
        The PubMed search term (e.g., '"SET gene"[Title/Abstract]')
    """
    if gene_name.upper() in AMBIGUOUS_GENE_NAMES:
        return f'"{gene_name} gene"[Title/Abstract]'
    else:
        return f'"{gene_name}"[Title/Abstract]'


def make_api_request_with_retry(url, params, timeout=30):
    """Makes an API request with an extended retry mechanism."""
    delays = [0.1] * 10 + [0.5] * 10 + [1] * 10 + [2, 5, 10, 20, 30, 60, 120]
    max_attempts = len(delays) + 1
    for i in range(max_attempts):
        try:
            if PUBMED_CALL_COOLDOWN > 0:
                time.sleep(PUBMED_CALL_COOLDOWN)
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response
        except (requests.exceptions.RequestException, requests.exceptions.HTTPError) as e:
            if VERBOSE_MODE:
                print(f"WARNING: API request failed (attempt {i+1}/{max_attempts}): {e}")
            if i < len(delays):
                time.sleep(delays[i])
            else:
                # Always surface critical errors
                if not VERBOSE_MODE:
                    print("ERROR: All API request attempts failed.")
                else:
                    print("ERROR: All API request attempts failed.")
                return None
    return None


def fetch_abstracts(pmids):
    """
    Fetch abstracts for a list of PMIDs using efetch.
    Returns a dictionary mapping PMID -> abstract text.
    """
    if not pmids:
        return {}
    
    efetch_url = f"{BASE_URL_NCBI}efetch.fcgi"
    
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "rettype": "abstract"
    }
    
    if NCBI_API_KEY:
        params['api_key'] = NCBI_API_KEY
    
    response = make_api_request_with_retry(efetch_url, params)
    if not response:
        return {}
    
    abstracts = {}
    
    try:
        root = ElementTree.fromstring(response.content)
        
        for article in root.findall('.//PubmedArticle'):
            pmid_elem = article.find('.//PMID')
            abstract_elem = article.find('.//Abstract/AbstractText')
            
            if pmid_elem is not None:
                pmid = pmid_elem.text
                if abstract_elem is not None:
                    abstract_text = "".join(abstract_elem.itertext())
                    abstracts[pmid] = abstract_text
                else:
                    abstracts[pmid] = ""
        
        if VERBOSE_MODE:
            print(f"INFO: Fetched abstracts for {len(abstracts)}/{len(pmids)} articles.")
        
    except Exception as e:
        print(f"ERROR: Could not parse abstract XML: {e}")
    
    return abstracts


def search_pubmed(gene_name, keywords, max_results):
    """
    Searches PubMed using a single combined query (always one-search mode).
    Deduplicates results and tracks how each article was found.
    """
    if VERBOSE_MODE:
        print("INFO: Performing combined single search.")
    
    # Use appropriate search term based on gene name ambiguity
    gene_term = get_pubmed_gene_term(gene_name)
    
    if keywords:
        keyword_clause = " OR ".join([f'"{kw}"' for kw in keywords])
        search_term = f'({gene_term}) AND ({keyword_clause})'
    else:
        search_term = gene_term
    
    search_queries = [{"term": search_term, "tag": "one_search", "retmax": max_results}]

    all_results = {}

    for query in search_queries:
        term, tag, retmax = query["term"], query["tag"], query["retmax"]

        esearch_params = {
            "db": "pubmed", "term": term, "retmax": retmax,
            "sort": "relevance", "usehistory": "y", "retmode": "json"
        }
        _apply_pubdate_cap(esearch_params)
        if NCBI_API_KEY:
            esearch_params['api_key'] = NCBI_API_KEY

        search_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esearch.fcgi", esearch_params)
        if not search_response:
            continue

        try:
            id_list = search_response.json().get("esearchresult", {}).get("idlist", [])
            if not id_list:
                continue

            esummary_params = {"db": "pubmed", "id": ",".join(id_list), "retmode": "json"}
            if NCBI_API_KEY:
                esummary_params['api_key'] = NCBI_API_KEY

            summary_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esummary.fcgi", esummary_params)
            if not summary_response:
                continue
            
            summary_data = summary_response.json()
            
            for pmid in id_list:
                if pmid not in all_results:
                    article_info = summary_data.get("result", {}).get(pmid, {})
                    title = article_info.get("title", "No title found")
                    pubdate = article_info.get("pubdate", "No date")
                    all_results[pmid] = {
                        "title": title, 
                        "pmid": pmid, 
                        "pubdate": pubdate,
                        "found_by": [tag],
                        "abstract": None
                    }
                else:
                    if tag not in all_results[pmid]["found_by"]:
                        all_results[pmid]["found_by"].append(tag)
        
        except Exception as e:
            print(f"ERROR: Could not parse PubMed response for term '{term}': {e}")
            continue

    final_results = _filter_articles_by_pubdate_cap(list(all_results.values()))
    if VERBOSE_MODE:
        print(f"INFO: Found a total of {len(final_results)} unique articles for {gene_name}.")
    
    return final_results


def search_pubmed_query(search_term, max_results=50):
    """
    Executes a single PubMed query string and returns summary records.
    """
    if VERBOSE_MODE:
        print(f"INFO: Custom PubMed query: {search_term}")

    esearch_params = {
        "db": "pubmed",
        "term": search_term,
        "retmax": max_results,
        "sort": "relevance",
        "retmode": "json"
    }
    _apply_pubdate_cap(esearch_params)
    if NCBI_API_KEY:
        esearch_params["api_key"] = NCBI_API_KEY

    search_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esearch.fcgi", esearch_params)
    if not search_response:
        return []

    try:
        id_list = search_response.json().get("esearchresult", {}).get("idlist", [])
    except Exception as exc:
        print(f"ERROR: Failed to parse custom query esearch response: {exc}")
        return []

    if not id_list:
        return []

    esummary_params = {"db": "pubmed", "id": ",".join(id_list), "retmode": "json"}
    if NCBI_API_KEY:
        esummary_params["api_key"] = NCBI_API_KEY

    summary_response = make_api_request_with_retry(f"{BASE_URL_NCBI}esummary.fcgi", esummary_params)
    if not summary_response:
        return []

    try:
        summary_data = summary_response.json()
    except Exception as exc:
        print(f"ERROR: Failed to parse custom query esummary response: {exc}")
        return []

    results = []
    for pmid in id_list:
        article_info = summary_data.get("result", {}).get(pmid, {})
        if not article_info:
            continue
        results.append(
            {
                "title": article_info.get("title", "No title found"),
                "pmid": pmid,
                "pubdate": article_info.get("pubdate", "No date"),
                "found_by": ["custom_query"],
                "abstract": None,
            }
        )

    return _filter_articles_by_pubdate_cap(results)


def search_lof_gof_pubmed(gene_name, max_results=10):
    """
    Search PubMed for LoF/GoF related articles for a gene.
    Returns a list of articles with pmid, title, and found_by tags.
    """
    esearch_url = f"{BASE_URL_NCBI}esearch.fcgi"
    esummary_url = f"{BASE_URL_NCBI}esummary.fcgi"
    
    # Use appropriate search term based on gene name ambiguity
    gene_term = get_pubmed_gene_term(gene_name)
    query = f"{gene_term} AND (gain-of-function OR loss-of-function)"
    
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "sort": "relevance"
    }
    _apply_pubdate_cap(params)
    
    if VERBOSE_MODE:
        print(f"INFO: Searching PubMed for LoF/GoF: {query}")
    
    search_response_raw = make_api_request_with_retry(esearch_url, params)
    if not search_response_raw:
        return []
    
    try:
        search_response = search_response_raw.json()
    except Exception as e:
        print(f"ERROR: Failed to parse search response JSON: {e}")
        return []
    
    pmids = search_response.get("esearchresult", {}).get("idlist", [])
    
    if not pmids:
        if VERBOSE_MODE:
            print(f"INFO: No LoF/GoF articles found for {gene_name}")
        return []
    
    if VERBOSE_MODE:
        print(f"INFO: Found {len(pmids)} LoF/GoF articles")
    
    summary_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json"
    }
    
    summary_response_raw = make_api_request_with_retry(esummary_url, summary_params)
    if not summary_response_raw:
        return []
    
    try:
        summary_response = summary_response_raw.json()
    except Exception as e:
        print(f"ERROR: Failed to parse summary response JSON: {e}")
        return []
    
    articles = []
    result_dict = summary_response.get("result", {})
    
    for pmid in pmids:
        if pmid in result_dict:
            article_data = result_dict[pmid]
            articles.append({
                "pmid": pmid,
                "title": article_data.get("title", "No title available"),
                "pubdate": article_data.get("pubdate", "No date"),
                "found_by": ["lof_gof_search"],
                "abstract": None
            })
    
    return _filter_articles_by_pubdate_cap(articles)


def deduplicate_articles(articles):
    """
    Deduplicates a list of articles based on title. If titles are the same,
    keeps the one with the most recent publication date.
    """
    def parse_pubdate(pubdate_str):
        """Tries to parse various PubMed date formats into a datetime object."""
        formats_to_try = [
            "%Y %b %d",
            "%Y %b",
            "%Y",
        ]
        for fmt in formats_to_try:
            try:
                return datetime.strptime(pubdate_str, fmt)
            except ValueError:
                continue
        return None

    unique_articles = {}
    for article in articles:
        title = article.get("title", "").strip().lower()
        if not title:
            continue
        
        current_date = parse_pubdate(article.get("pubdate", ""))

        if title not in unique_articles:
            unique_articles[title] = article
        else:
            existing_article = unique_articles[title]
            existing_date = parse_pubdate(existing_article.get("pubdate", ""))
            
            if current_date and not existing_date:
                unique_articles[title] = article
            elif current_date and existing_date and current_date > existing_date:
                unique_articles[title] = article
    
    return list(unique_articles.values())


def format_articles_for_llm(articles, use_abstracts=False, top_abstracts=None):
    """
    Format a list of articles for LLM input.
    
    Args:
        articles: List of article dictionaries
        use_abstracts: If True, includes abstract text for all articles
        top_abstracts: If set (int), only include abstracts for the first N articles.
                       This overrides use_abstracts for a hybrid mode.
                       Example: top_abstracts=10 means first 10 get abstracts, rest get titles only.
    """
    if not articles:
        return "No articles found."
    
    formatted_lines = []
    for i, article in enumerate(articles):
        title = article.get('title', 'No title')
        pmid = article.get('pmid', 'Unknown')
        found_by = article.get('found_by', [])
        abstract = article.get('abstract', '')
        
        line = f"- {title} (PMID: {pmid}) [Found by: {', '.join(found_by)}]"
        
        # Determine if we should include abstract for this article
        include_abstract = False
        if top_abstracts is not None:
            # Hybrid mode: only top N articles get abstracts
            include_abstract = (i < top_abstracts) and abstract
        elif use_abstracts:
            # All abstracts mode
            include_abstract = bool(abstract)
        
        if include_abstract:
            line += f"\n  Abstract: {abstract}"
        
        formatted_lines.append(line)
    
    return "\n".join(formatted_lines)

