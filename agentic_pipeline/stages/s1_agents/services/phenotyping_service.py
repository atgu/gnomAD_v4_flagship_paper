"""Service for phenotyping contamination analysis."""
import pandas as pd
import os
import time
import json
import subprocess
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils.helpers import load_prompt, sanitize_filename
from services.llm_service import call_llm
from config import APP_ROOT


def ensure_data_file_exists(filename):
    """
    Ensures a data file exists locally, downloading from GCS if needed.
    
    Args:
        filename (str): The filename (e.g., 'code_endpoint_counts_DF13_v3.tsv')
    
    Returns:
        str: Path to the local file
    """
    local_path = os.path.join(APP_ROOT, "data", "phenotyping", filename)
    
    # If file already exists, return its path
    if os.path.exists(local_path):
        return local_path
    
    # Download from public GCS URL
    public_url = f"https://storage.googleapis.com/llm_agents_bucket/{filename}"
    print(f"Downloading {filename} from {public_url}...")
    
    try:
        import requests
        
        response = requests.get(public_url, stream=True)
        response.raise_for_status()
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        
        # Download with progress
        with open(local_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"Successfully downloaded {filename}")
        return local_path
        
    except Exception as e:
        print(f"Error downloading {filename} from public URL: {str(e)}")
        raise FileNotFoundError(f"Could not download {filename} from {public_url}")
    except FileNotFoundError:
        print("Error: 'gsutil' command not found. Please install Google Cloud SDK.")
        raise


def get_endpoint_longname(endpoint_name):
    """
    Retrieves the LONGNAME for a given endpoint from the endpoints definition file.
    
    Args:
        endpoint_name (str): The endpoint name (e.g., 'G6_PARKINSON')
        
    Returns:
        str: The LONGNAME or None if not found
    """
    endpoints_file = os.path.join(APP_ROOT, "data", "phenotyping", "endpoints_definition.xlsx")
    
    try:
        df_endpoints = pd.read_excel(endpoints_file)
        match = df_endpoints[df_endpoints['NAME'] == endpoint_name]
        if not match.empty:
            return match.iloc[0]['LONGNAME']
        else:
            print(f"Warning: Endpoint '{endpoint_name}' not found in endpoints file")
            return None
    except Exception as e:
        print(f"Error reading endpoints file: {e}")
        return None


def call_claude_contamination_check_batch(endpoint_longname, name_en_list, endpoint_name="unknown"):
    """
    Calls Claude to evaluate multiple enrichments in a single request.
    
    Args:
        endpoint_longname (str): Full descriptive name of the GWAS endpoint
        name_en_list (list): List of enriched phenotype names
        endpoint_name (str): Short endpoint name for debug files
        
    Returns:
        list: List of dicts with {'score': int, 'justification': str} or error values
    """
    from parsers.llm_parsers import parse_contamination_response
    
    # Create numbered list of phenotypes
    phenotypes_text = "\n".join([f"{i+1}. {name_en}" for i, name_en in enumerate(name_en_list)])
    
    # Load prompt template
    prompt_template = load_prompt('phenotyping_contamination_prompt')
    
    if not prompt_template:
        print("ERROR: Could not load phenotyping contamination prompt")
        return [{'score': 'NA', 'justification': 'Prompt template not found'} for _ in name_en_list]
    
    prompt = prompt_template.format(
        endpoint_longname=endpoint_longname,
        phenotypes_text=phenotypes_text
    )
    
    try:
        # Use existing LLM service
        response_text = call_llm(
            prompt, model="claude-haiku-4-5", temperature=0.0, agent_name="contamination_check_agent"
        )
        
        if not response_text:
            return [{'score': 'NA', 'justification': 'Empty LLM response'} for _ in name_en_list]
        
        print(f"DEBUG: LLM Response (first 500 chars):\n{response_text[:500]}")
        
        # Parse the response
        results = parse_contamination_response(response_text, name_en_list, endpoint_name)
        
        print(f"DEBUG: Parsed {len(results)} results")
        for i, r in enumerate(results[:3]):  # Show first 3
            print(f"  Result {i+1}: score={r['score']}, justification={r['justification'][:50]}...")
        
        return results
            
    except Exception as e:
        print(f"ERROR: LLM call failed: {e}")
        return [{'score': 'NA', 'justification': f'API error: {str(e)}'} for _ in name_en_list]


def load_and_filter_endpoint_data(endpoint_name, max_rows=None):
    """
    Loads and filters endpoint data from TSV file.
    
    Args:
        endpoint_name (str): The endpoint to filter by
        max_rows (int): Maximum rows to return after filtering
        
    Returns:
        pd.DataFrame: Filtered and enriched data
    """
    # Ensure the large data file is available (download from GCS if needed)
    counts_file = ensure_data_file_exists("code_endpoint_counts_DF13_v3.tsv")
    codes_info_file = os.path.join(APP_ROOT, "data", "phenotyping", "fg_codes_info_v9.tsv")
    
    print(f"Loading counts data for endpoint: {endpoint_name}")
    print(f"Counts file path: {counts_file}")
    
    # Load data in chunks for memory efficiency
    chunk_iter = pd.read_csv(counts_file, sep='\t', chunksize=100000)
    
    filtered_chunks = []
    for chunk in chunk_iter:
        filtered_chunks.append(chunk[chunk['ENDPOINT'] == endpoint_name])
    
    if not filtered_chunks:
        raise ValueError(f"No data found for endpoint '{endpoint_name}'")
    
    df_filtered = pd.concat(filtered_chunks, ignore_index=True)
    print(f"Found {len(df_filtered)} rows for endpoint '{endpoint_name}'")
    
    # Filter by counts_mlogp > 10
    df_filtered['counts_mlogp'] = pd.to_numeric(df_filtered['counts_mlogp'], errors='coerce')
    df_filtered = df_filtered.dropna(subset=['counts_mlogp'])
    df_filtered = df_filtered[df_filtered['counts_mlogp'] > 10]
    
    if df_filtered.empty:
        raise ValueError(f"No data found for endpoint '{endpoint_name}' with counts_mlogp > 10")
    
    print(f"Found {len(df_filtered)} rows after filtering by counts_mlogp")
    
    # Load codes info
    print("Loading code information...")
    df_codes = pd.read_csv(codes_info_file, sep='\t')
    
    # Enrich data
    print("Enriching data with code names...")
    df_enriched = enrich_with_code_names(df_filtered, df_codes)
    
    # Sort and limit
    df_enriched = df_enriched.sort_values(by='counts_mlogp', ascending=False).reset_index(drop=True)
    
    if max_rows and max_rows > 0:
        original_count = len(df_enriched)
        df_enriched = df_enriched.head(max_rows)
        print(f"Limited to top {max_rows} rows (from {original_count} total)")
    
    return df_enriched


def enrich_with_code_names(df_filtered, df_codes):
    """
    Enriches filtered data with human-readable code names using hierarchical merge.
    
    Args:
        df_filtered (pd.DataFrame): Filtered endpoint data
        df_codes (pd.DataFrame): Code information dictionary
        
    Returns:
        pd.DataFrame: Enriched data with name_en column
    """
    # Replace string 'NA' with actual NaN values
    df_filtered.replace('NA', pd.NA, inplace=True)
    df_codes.replace('NA', pd.NA, inplace=True)
    
    # Define merge keys
    merge_keys_base = ['FG_CODE1', 'vocabulary_id']
    merge_keys_c2 = merge_keys_base + ['FG_CODE2']
    merge_keys_c3 = merge_keys_c2 + ['FG_CODE3']
    
    # Split codes dictionary by specificity
    codes_c1 = df_codes[df_codes['FG_CODE2'].isna() & df_codes['FG_CODE3'].isna()]
    codes_c2 = df_codes[df_codes['FG_CODE2'].notna() & df_codes['FG_CODE3'].isna()]
    codes_c3 = df_codes[df_codes['FG_CODE2'].notna() & df_codes['FG_CODE3'].notna()]
    
    # Add temporary unique ID to track rows
    df_filtered['row_id'] = range(len(df_filtered))
    
    # Hierarchical merge - most specific first
    # Stage 1: Match 3 codes
    df_m3 = pd.merge(df_filtered, codes_c3[merge_keys_c3 + ['name_en']], on=merge_keys_c3, how='inner')
    
    # Stage 2: Match 2 codes on remaining rows
    remaining_after_m3 = df_filtered[~df_filtered['row_id'].isin(df_m3['row_id'])]
    df_m2 = pd.merge(remaining_after_m3, codes_c2[merge_keys_c2 + ['name_en']], on=merge_keys_c2, how='inner')
    
    # Stage 3: Match 1 code on remaining rows
    remaining_after_m2 = remaining_after_m3[~remaining_after_m3['row_id'].isin(df_m2['row_id'])]
    df_m1 = pd.merge(remaining_after_m2, codes_c1[merge_keys_base + ['name_en']], on=merge_keys_base, how='inner')
    
    # Combine all matched parts
    df_enriched = pd.concat([df_m3, df_m2, df_m1]).sort_values('row_id').drop(columns=['row_id'])
    
    # Restore 'NA' strings for consistent output
    df_enriched.fillna('NA', inplace=True)
    
    return df_enriched


def _normalize_name(name: str) -> str:
    """
    Normalize phenotype/drug/procedure names for matching:
    - strip leading numbering like '1. ' or '12) '
    - trim whitespace
    - lowercase for comparisons
    """
    if not name:
        return ""
    # Remove leading numbering/punctuation
    name = re.sub(r"^\s*[\dIVXivx]+\s*[.)-]?\s*", "", name)
    return name.strip().lower()


def _canonical_category(cat: str) -> str:
    """
    Map arbitrary category text to one of the canonical keys.
    """
    if not cat:
        return "Other"
    c = cat.strip().lower()
    if "broad" in c:
        return "Broad disease categories"
    if "procedure" in c:
        return "Medical procedures"
    if "drug" in c or "therapy" in c or "treatment" in c:
        return "Drugs/therapies"
    if "symptom" in c or "sign" in c:
        return "Symptoms"
    if "disease" in c or "syndrome" in c:
        # Prefer specific diseases over broad
        if "broad" in c or "category" in c:
            return "Broad disease categories"
        return "Diseases or syndromes"
    return "Other"


def _write_debug_output(agent_key, safe_name, batch_index, prompt, llm_text):
    """
    Write input/output of an agent call for diagnostics.
    """
    runs_dir = os.path.join(APP_ROOT, "data", "runs", "phenotyping")
    os.makedirs(runs_dir, exist_ok=True)
    agent_fs_name = sanitize_filename(agent_key) or "agent"
    fname = os.path.join(
        runs_dir, f"debug_{agent_fs_name}_batch{batch_index:03d}.txt"
    )
    try:
        with open(fname, "w", encoding="utf-8") as f:
            f.write("PROMPT:\n")
            f.write(prompt or "")
            f.write("\n\nLLM_OUTPUT:\n")
            f.write(llm_text or "")
        print(f"[DEBUG] Wrote agent I/O to {fname}")
    except Exception as exc:
        print(f"[DEBUG] Failed to write debug file {fname}: {exc}")


def _detect_parallel_workers(max_cap=8, fallback=4):
    """
    Detect a reasonable number of workers based on CPU count, capped to max_cap.
    """
    try:
        cpu = os.cpu_count() or fallback
    except Exception:
        cpu = fallback
    return max(1, min(max_cap, cpu))


def run_batch_category_agent(
    endpoint_longname,
    name_en_list,
    endpoint_name="unknown",
    batch_size=10,
    parallel_workers=None,
):
    """
    Batch agent that classifies each descriptor into one of:
    Diseases or syndromes, Broad disease categories, Symptoms,
    Drugs/therapies, Medical procedures, Other.
    Outputs a plain-text file with 'descriptor: category' per line.
    """
    prompt_template = load_prompt("phenotyping_batch_category_agent")
    if not prompt_template:
        print("ERROR: phenotyping_batch_category_agent prompt missing")
        return None

    safe_name = sanitize_filename(endpoint_name) or "endpoint"
    runs_dir = os.path.join(APP_ROOT, "data", "runs", "phenotyping")
    os.makedirs(runs_dir, exist_ok=True)
    categories_path = os.path.join(runs_dir, f"{safe_name}_categories.txt")

    if parallel_workers is None:
        parallel_workers = _detect_parallel_workers()

    if parallel_workers is None:
        parallel_workers = _detect_parallel_workers()

    batches = []
    for i in range(0, len(name_en_list), batch_size):
        batches.append((i // batch_size, name_en_list[i : i + batch_size]))

    if not batches:
        return {"categories": []}

    def _process_batch(batch_index, batch):
        descriptors_text = "\n".join(batch)  # raw, one per line
        prompt = prompt_template.format(
            descriptors_text=descriptors_text,
        )

        try:
            llm_text = call_llm(
                prompt,
                model="claude-haiku-4-5",
                temperature=0.0,
                agent_name="phenotype_batch_category_agent",
            )
        except Exception as exc:
            print(f"ERROR: Batch category agent LLM call failed: {exc}")
            llm_text = ""

        _write_debug_output("phenotype_batch_category_agent", safe_name, batch_index, prompt, llm_text)

        parsed = {}
        if llm_text:
            for raw_line in llm_text.splitlines():
                line = raw_line.strip()
                if not line or ":" not in line:
                    continue
                desc, rest = line.split(":", 1)
                desc = desc.strip()
                cat_part = rest.strip()
                justification = ""
                if " - " in cat_part:
                    cat_part, justification = cat_part.split(" - ", 1)
                    cat_part = cat_part.strip()
                    justification = justification.strip()
                cat = _canonical_category(cat_part)
                full_line = f"{desc}: {cat}"
                if justification:
                    full_line = f"{full_line} - {justification}"
                parsed[_normalize_name(desc)] = full_line

        entries = []
        for item in batch:
            norm = _normalize_name(item)
            if norm in parsed:
                entries.append(parsed[norm])
            else:
                entries.append(f"{item}: Other - No justification provided")

        return {"batch_index": batch_index, "entries": entries}

    results = []
    max_workers = parallel_workers if parallel_workers and parallel_workers > 0 else 1
    max_workers = min(max_workers, len(batches))

    if max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_process_batch, idx, batch) for idx, batch in batches]
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:
                    print(f"ERROR: Batch category agent batch failed: {exc}")
    else:
        for idx, batch in batches:
            results.append(_process_batch(idx, batch))

    all_entries = []
    for res in sorted(results, key=lambda r: r.get("batch_index", 0)):
        all_entries.extend(res.get("entries", []))

    # Deduplicate while preserving order
    def _unique(seq):
        seen = set()
        out = []
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    all_entries = _unique(all_entries)

    try:
        with open(categories_path, "w", encoding="utf-8") as f_cat:
            f_cat.write("\n".join(all_entries))
        print(f"Category agent output written to {categories_path}")
    except Exception as exc:
        print(f"ERROR: Failed to write category agent output: {exc}")

    return {"categories": all_entries}


def run_category_split_selection(
    endpoint_longname,
    name_en_list,
    endpoint_name="unknown",
    batch_size=10,
    parallel_workers=None,
):
    """
    Split descriptors by category (using the category file) and run 6 specialized agents,
    then produce unified selected/rejected lists with justifications.
    """
    if parallel_workers is None:
        parallel_workers = _detect_parallel_workers()
    safe_name = sanitize_filename(endpoint_name) or "endpoint"
    runs_dir = os.path.join(APP_ROOT, "data", "runs", "phenotyping")
    categories_path = os.path.join(runs_dir, f"{safe_name}_categories.txt")

    # Load categories (desc: category) from file if present; else classify on the fly
    categories_map = {}
    if os.path.isfile(categories_path):
        try:
            with open(categories_path, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line or ":" not in line:
                        continue
                    desc, rest = line.split(":", 1)
                    desc_norm = _normalize_name(desc)
                    cat_part = rest.strip()
                    if " - " in cat_part:
                        cat_part = cat_part.split(" - ", 1)[0].strip()
                    categories_map[desc_norm] = _canonical_category(cat_part)
        except Exception as exc:
            print(f"WARNING: Failed to read categories file: {exc}")

    # If missing entries, classify quickly via category agent (reusing existing)
    missing = [d for d in name_en_list if _normalize_name(d) not in categories_map]
    if missing:
        cat_result = run_batch_category_agent(
            endpoint_longname,
            missing,
            endpoint_name,
            batch_size=batch_size,
            parallel_workers=parallel_workers,
        )
        # run_batch_category_agent writes file; also fill map
        if cat_result and "categories" in cat_result:
            for line in cat_result["categories"]:
                if ":" not in line:
                    continue
                desc, rest = line.split(":", 1)
                desc_norm = _normalize_name(desc)
                cat_part = rest.strip()
                if " - " in cat_part:
                    cat_part = cat_part.split(" - ", 1)[0].strip()
                categories_map[desc_norm] = _canonical_category(cat_part)

    # Buckets
    buckets = {
        "Diseases or syndromes": [],
        "Broad disease categories": [],
        "Symptoms": [],
        "Drugs/therapies": [],
        "Medical procedures": [],
        "Other": [],
    }
    for d in name_en_list:
        norm = _normalize_name(d)
        cat = categories_map.get(norm, "Other")
        buckets.setdefault(cat, [])
        buckets[cat].append(d)

    # Agent prompts mapping
    prompt_map = {
        "Diseases or syndromes": "phenotyping_filter_diseases",
        "Broad disease categories": "phenotyping_filter_broad",
        "Symptoms": "phenotyping_filter_symptoms",
        "Drugs/therapies": "phenotyping_filter_drugs",
        "Medical procedures": "phenotyping_filter_procedures",
        "Other": "phenotyping_filter_other",
    }

    selected_all = []
    ambiguous_all = []
    rejected_all = []

    for cat_name, items in buckets.items():
        if not items:
            continue
        prompt_name = prompt_map.get(cat_name)
        if not prompt_name:
            continue
        agent_key = f"phenotype_filter_{_normalize_name(cat_name) or 'misc'}"
        sel, amb, rej = _run_category_specific_agent(
            agent_key,
            prompt_name,
            endpoint_longname,
            items,
            batch_size,
            safe_name,
            parallel_workers=parallel_workers,
        )
        selected_all.extend(sel)
        ambiguous_all.extend(amb)
        rejected_all.extend(rej)

    # Deduplicate preserving order; remove rejected that are selected
    def _unique(seq):
        seen = set()
        out = []
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    selected_all = _unique(selected_all)
    ambiguous_all = _unique(ambiguous_all)
    selected_name_set = {
        _normalize_name(line.split(":", 1)[0] if ":" in line else line)
        for line in selected_all
    }
    rejected_all = _unique(
        [r for r in rejected_all if _normalize_name(r.split(':',1)[0]) not in selected_name_set]
    )

    # Write outputs
    selected_path = os.path.join(runs_dir, f"{safe_name}_selected.txt")
    ambiguous_path = os.path.join(runs_dir, f"{safe_name}_ambiguous.txt")
    rejected_path = os.path.join(runs_dir, f"{safe_name}_rejected.txt")
    try:
        with open(selected_path, "w", encoding="utf-8") as f_sel:
            f_sel.write("\n".join(selected_all))
        with open(ambiguous_path, "w", encoding="utf-8") as f_amb:
            f_amb.write("\n".join(ambiguous_all))
        with open(rejected_path, "w", encoding="utf-8") as f_rej:
            f_rej.write("\n".join(rejected_all))
        print(f"Category split filter outputs written to {selected_path}, {ambiguous_path}, and {rejected_path}")
    except Exception as exc:
        print(f"ERROR: Failed to write split filter outputs: {exc}")

    return {"selected": selected_all, "ambiguous": ambiguous_all, "rejected": rejected_all}


def _parse_selected_rejected(output_text):
    """
    Parse output with three sections:
    SELECTED:
    <descriptor>: <justification>
    AMBIGUOUS:
    <descriptor>: <justification>
    REJECTED:
    <descriptor>: <justification>
    Returns (selected_entries, ambiguous_entries, rejected_entries) where each entry is (name, line).
    """
    if not output_text:
        return [], [], []

    selected = []
    ambiguous = []
    rejected = []
    mode = None

    for raw_line in output_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        lower = line.lower()
        if lower.startswith("selected"):
            mode = "selected"
            continue
        if lower.startswith("ambiguous"):
            mode = "ambiguous"
            continue
        if lower.startswith("rejected"):
            mode = "rejected"
            continue

        if mode == "selected":
            if ":" in line:
                phen, _ = line.split(":", 1)
                phen = phen.strip()
            else:
                phen = line
            selected.append((phen, line))
        elif mode == "ambiguous":
            if ":" in line:
                phen, _ = line.split(":", 1)
                phen = phen.strip()
            else:
                phen = line
            ambiguous.append((phen, line))
        elif mode == "rejected":
            if ":" in line:
                phen, _ = line.split(":", 1)
                phen = phen.strip()
            else:
                phen = line
            rejected.append((phen, line))

    return selected, ambiguous, rejected


def run_batch_selection_agent(endpoint_longname, name_en_list, endpoint_name="unknown", batch_size=10):
    """
    Batch agent that classifies phenotypes as exact synonyms, major symptoms,
    or associated drugs for the target phenotype. Produces two plain-text files:
    selected and not-selected.
    """
    prompt_template = load_prompt("phenotyping_batch_filter_agent")
    if not prompt_template:
        print("ERROR: phenotyping_batch_filter_agent prompt missing")
        return None

    safe_name = sanitize_filename(endpoint_name) or "endpoint"
    runs_dir = os.path.join(APP_ROOT, "data", "runs", "phenotyping")
    os.makedirs(runs_dir, exist_ok=True)
    selected_path = os.path.join(runs_dir, f"{safe_name}_selected.txt")
    rejected_path = os.path.join(runs_dir, f"{safe_name}_rejected.txt")

    selected_all = []
    rejected_all = []

    for i in range(0, len(name_en_list), batch_size):
        batch = name_en_list[i : i + batch_size]
        # Send raw phenotypes, one per line, without numbering
        phenotypes_text = "\n".join(batch)
        prompt = prompt_template.format(
            endpoint_longname=endpoint_longname,
            phenotypes_text=phenotypes_text,
        )

        try:
            llm_text = call_llm(
                prompt,
                model="claude-haiku-4-5",
                temperature=0.0,
                agent_name="phenotype_batch_filter_agent",
            )
        except Exception as exc:
            print(f"ERROR: Batch selection agent LLM call failed: {exc}")
            llm_text = ""

        batch_selected_entries, batch_ambiguous_entries, batch_rejected_entries = _parse_selected_rejected(llm_text)

        # Normalize names for comparison
        selected_name_set = {_normalize_name(p) for p, _ in batch_selected_entries}
        ambiguous_name_set = {_normalize_name(p) for p, _ in batch_ambiguous_entries}
        rejected_name_set = {_normalize_name(p) for p, _ in batch_rejected_entries}

        # Everything not explicitly selected is rejected; keep justification if provided
        batch_rejected = []
        for item in batch:
            norm = _normalize_name(item)
            if norm in selected_name_set:
                continue
            if norm in ambiguous_name_set:
                # Use the line from ambiguous_entries
                for p, line in batch_ambiguous_entries:
                    if _normalize_name(p) == norm:
                        # Treat ambiguous as rejected downstream but keep line
                        batch_rejected.append(line)
                        break
                continue
            if norm in rejected_name_set:
                # Use the line from rejected_entries
                for p, line in batch_rejected_entries:
                    if _normalize_name(p) == norm:
                        batch_rejected.append(line)
                        break
            else:
                batch_rejected.append(f"{item}: No justification provided")

        # Store lines (descriptor: justification)
        selected_all.extend([line for _p, line in batch_selected_entries])
        rejected_all.extend(batch_rejected)

    # Deduplicate while preserving order
    def _unique(seq):
        seen = set()
        out = []
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    selected_all = _unique(selected_all)
    # Build a set of normalized names
    selected_name_set = {
        _normalize_name(line.split(":", 1)[0] if ":" in line else line)
        for line in selected_all
    }
    ambiguous_name_set = {
        _normalize_name(line.split(":", 1)[0] if ":" in line else line)
        for line in ambiguous_all
    }
    rejected_all = _unique(
        [
            r
            for r in rejected_all
            if _normalize_name(r.split(":", 1)[0]) not in selected_name_set
            and _normalize_name(r.split(":", 1)[0]) not in ambiguous_name_set
        ]
    )

    ambiguous_path = os.path.join(runs_dir, f"{safe_name}_ambiguous.txt")

    try:
        with open(selected_path, "w", encoding="utf-8") as f_sel:
            f_sel.write("\n".join(selected_all))
        with open(ambiguous_path, "w", encoding="utf-8") as f_amb:
            f_amb.write("\n".join(ambiguous_all))
        with open(rejected_path, "w", encoding="utf-8") as f_rej:
            f_rej.write("\n".join(rejected_all))
        print(f"Selection agent outputs written to {selected_path}, {ambiguous_path}, and {rejected_path}")
    except Exception as exc:
        print(f"ERROR: Failed to write selection agent outputs: {exc}")

    return {"selected": selected_all, "rejected": rejected_all}


def _run_category_specific_agent(
    agent_key,
    prompt_template_name,
    endpoint_longname,
    descriptors,
    batch_size,
    safe_name,
    parallel_workers=None,
):
    """
    Generic runner for category-specific agents that output SELECTED/REJECTED sections.
    Returns lists of lines with justification.
    """
    prompt_template = load_prompt(prompt_template_name)
    if not prompt_template:
        print(f"ERROR: Prompt missing for {agent_key}")
        return [], [], []

    batches = []
    for i in range(0, len(descriptors), batch_size):
        batches.append((i // batch_size, descriptors[i : i + batch_size]))

    if not batches:
        return [], [], []

    def _process_batch(batch_index, batch):
        descriptors_text = "\n".join(batch)
        prompt = prompt_template.format(
            endpoint_longname=endpoint_longname,
            descriptors_text=descriptors_text,
        )

        try:
            llm_text = call_llm(
                prompt,
                model="claude-haiku-4-5",
                temperature=0.0,
                agent_name=agent_key,
            )
        except Exception as exc:
            print(f"ERROR: {agent_key} LLM call failed: {exc}")
            llm_text = ""

        _write_debug_output(agent_key, safe_name, batch_index, prompt, llm_text)

        batch_selected, batch_ambiguous, batch_rejected = _parse_selected_rejected(llm_text)
        selected_lines = [line for _p, line in batch_selected]
        ambiguous_lines = [line for _p, line in batch_ambiguous]
        rejected_lines = [line for _p, line in batch_rejected]

        # Add any missing descriptors as ambiguous with default note
        selected_name_set = {_normalize_name(p) for p, _ in batch_selected}
        ambiguous_name_set = {_normalize_name(p) for p, _ in batch_ambiguous}
        rejected_name_set = {_normalize_name(p) for p, _ in batch_rejected}
        for item in batch:
            norm = _normalize_name(item)
            if norm in selected_name_set or norm in ambiguous_name_set or norm in rejected_name_set:
                continue
            ambiguous_lines.append(f"{item}: No justification provided")

        return {
            "batch_index": batch_index,
            "selected": selected_lines,
            "ambiguous": ambiguous_lines,
            "rejected": rejected_lines,
        }

    results = []
    max_workers = parallel_workers if parallel_workers and parallel_workers > 0 else 1
    max_workers = min(max_workers, len(batches))

    if max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_process_batch, idx, batch) for idx, batch in batches]
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:
                    print(f"ERROR: {agent_key} batch failed: {exc}")
    else:
        for idx, batch in batches:
            results.append(_process_batch(idx, batch))

    selected = []
    ambiguous = []
    rejected = []
    for res in sorted(results, key=lambda r: r.get("batch_index", 0)):
        selected.extend(res.get("selected", []))
        ambiguous.extend(res.get("ambiguous", []))
        rejected.extend(res.get("rejected", []))

    return selected, ambiguous, rejected


def run_full_analysis(analysis_id, endpoint_name, max_rows, analysis_status, analysis_results, batch_size=10):
    """
    Runs the complete analysis pipeline with progress tracking.
    
    Args:
        analysis_id (str): Unique analysis identifier
        endpoint_name (str): Endpoint to analyze
        max_rows (int): Maximum phenotypes to analyze
        analysis_status (dict): Global status dictionary
        analysis_results (dict): Global results dictionary
        batch_size (int): Batch size for LLM calls
    """
    try:
        worker_count = _detect_parallel_workers()

        analysis_status[analysis_id] = {
            'status': 'running',
            'progress': 0,
            'message': 'Starting analysis...'
        }
        
        # Step 1: Load and filter data
        analysis_status[analysis_id]['message'] = 'Loading and filtering data...'
        analysis_status[analysis_id]['progress'] = 20
        
        max_rows_requested = max_rows
        df_enriched = load_and_filter_endpoint_data(endpoint_name, None)
        
        # Step 2: Get endpoint longname
        analysis_status[analysis_id]['message'] = 'Getting endpoint definition...'
        analysis_status[analysis_id]['progress'] = 40
        
        endpoint_longname = get_endpoint_longname(endpoint_name)
        if not endpoint_longname:
            endpoint_longname = endpoint_name
        
        # Step 2a: Compute min_n_cases based on reference endpoint phenotype
        analysis_status[analysis_id]['message'] = 'Computing minimum cases threshold...'
        analysis_status[analysis_id]['progress'] = 43
        # Use n_cases_total for the reference endpoint; fallback to 200 if missing
        df_enriched['n_cases_total_num'] = pd.to_numeric(df_enriched.get('n_cases_total'), errors='coerce')
        endpoint_n_cases = None
        if not df_enriched.empty:
            # df_enriched is already filtered to the endpoint; use the max non-null total cases
            endpoint_n_cases = df_enriched['n_cases_total_num'].dropna().max()
        if pd.isna(endpoint_n_cases) or endpoint_n_cases is None:
            min_n_cases = 200
        else:
            min_n_cases = int(round(0.05 * float(endpoint_n_cases)))
        df_enriched['n_cases_num'] = pd.to_numeric(df_enriched.get('n_cases_with_code'), errors='coerce')
        df_enriched = df_enriched[df_enriched['n_cases_num'] >= min_n_cases]
        df_enriched = df_enriched.drop(columns=['n_cases_num', 'n_cases_total_num'])
        if df_enriched.empty:
            raise ValueError(f"No phenotypes remain after applying min_n_cases={min_n_cases}")

        if max_rows_requested and max_rows_requested > 0:
            original_count = len(df_enriched)
            df_enriched = df_enriched.head(max_rows_requested)
            print(f"Limited to top {max_rows_requested} rows (from {original_count} after min_n_cases filter)")

        # Step 2b: Categorize descriptors, then run category-specific filters
        analysis_status[analysis_id]['message'] = 'Categorizing descriptors...'
        analysis_status[analysis_id]['progress'] = 45
        name_en_list = df_enriched['name_en'].tolist()
        run_batch_category_agent(
            endpoint_longname,
            name_en_list,
            endpoint_name,
            batch_size=batch_size,
            parallel_workers=worker_count,
        )

        analysis_status[analysis_id]['message'] = 'Filtering by category...'
        analysis_status[analysis_id]['progress'] = 48
        run_category_split_selection(
            endpoint_longname,
            name_en_list,
            endpoint_name,
            batch_size=batch_size,
            parallel_workers=worker_count,
        )
        
        # Step 3: Build results from category split (accepted/rejected lists)
        analysis_status[analysis_id]['message'] = 'Finalizing results...'
        analysis_status[analysis_id]['progress'] = 95
        
        # Load selected/rejected with justifications
        safe_name = sanitize_filename(endpoint_name) or "endpoint"
        runs_dir = os.path.join(APP_ROOT, "data", "runs", "phenotyping")
        selected_path = os.path.join(runs_dir, f"{safe_name}_selected.txt")
        ambiguous_path = os.path.join(runs_dir, f"{safe_name}_ambiguous.txt")
        rejected_path = os.path.join(runs_dir, f"{safe_name}_rejected.txt")

        def _load_lines(path):
            if not os.path.isfile(path):
                return []
            lines = []
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    lines.append(line)
            return lines

        selected_lines = _load_lines(selected_path)
        ambiguous_lines = _load_lines(ambiguous_path)
        rejected_lines = _load_lines(rejected_path)

        # Load categories
        categories_path = os.path.join(runs_dir, f"{safe_name}_categories.txt")
        categories_lines = _load_lines(categories_path)
        category_map = {}
        for line in categories_lines:
            if ":" in line:
                desc, rest = line.split(":", 1)
                cat_part = rest.strip()
                if " - " in cat_part:
                    cat_part = cat_part.split(" - ", 1)[0].strip()
                category_map[desc.strip().lower()] = cat_part

        def _parse_line(line):
            if ":" in line:
                name, rest = line.split(":", 1)
                return name.strip(), rest.strip()
            return line.strip(), ""

        status_map = {}
        for line in rejected_lines:
            name, just = _parse_line(line)
            status_map[name.lower()] = ("rejected", just)
        for line in ambiguous_lines:
            name, just = _parse_line(line)
            status_map[name.lower()] = ("ambiguous", just)
        for line in selected_lines:
            name, just = _parse_line(line)
            status_map[name.lower()] = ("accepted", just)

        # Prepare results from df_enriched, keeping counts
        df_enriched['name_key'] = df_enriched['name_en'].str.lower()
        # Sort inside each status by n_cases_with_code (desc)
        def _get_status(row):
            return status_map.get(row['name_key'], ("rejected", ""))

        def _get_category(row):
            return category_map.get(row['name_key'], "Other")

        results = []
        for _, row in df_enriched.iterrows():
            status, justification = _get_status(row)
            category = _get_category(row)
            results.append({
                'name_en': row['name_en'],
                'status': status,
                'category': category,
                'justification': justification,
                'counts_mlogp': float(row['counts_mlogp']),
                'n_cases_with_code': row.get('n_cases_with_code', 'NA'),
                'n_controls_with_code': row.get('n_controls_with_code', 'NA')
            })

        # Split and sort
        rejected_items = [r for r in results if r['status'] == 'rejected']
        ambiguous_items = [r for r in results if r['status'] == 'ambiguous']
        accepted_items = [r for r in results if r['status'] == 'accepted']
        def _sort_key(item):
            try:
                return -float(item.get('n_cases_with_code', 0))
            except Exception:
                return 0
        rejected_items = sorted(rejected_items, key=_sort_key)
        ambiguous_items = sorted(ambiguous_items, key=_sort_key)
        accepted_items = sorted(accepted_items, key=_sort_key)
        results = rejected_items + ambiguous_items + accepted_items
        
        from datetime import datetime
        analysis_results[analysis_id] = {
            'endpoint_name': endpoint_name,
            'endpoint_longname': endpoint_longname,
            'total_results': len(results),
            'results': results,
            'timestamp': datetime.now().isoformat()
        }
        
        analysis_status[analysis_id] = {
            'status': 'completed',
            'progress': 100,
            'message': f'Analysis completed! Found {len(results)} results.'
        }
        
    except Exception as e:
        analysis_status[analysis_id] = {
            'status': 'error',
            'message': f'Error during analysis: {str(e)}'
        }
        print(f"ERROR in analysis {analysis_id}: {e}")

