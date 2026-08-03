"""Service for comparing diseases between Deep Analysis and GenCC."""
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple
from time import perf_counter

from utils.helpers import load_prompt
from services.llm_service import call_llm

# Global flag to control verbose output
VERBOSE_MODE = False

# GenCC evidence order (strongest to weakest)
EVIDENCE_ORDER = [
    'Definitive', 'Strong', 'Moderate', 'Supportive', 'Limited',
    'Disputed', 'Disputed Evidence', 'Refuted', 'Refuted Evidence',
    'No Known Disease Relationship'
]

def _sort_classifications_by_evidence(classifications: list) -> list:
    """
    Sort the classifications by decreasing evidence.
    Unknown classifications are pushed to the end.
    """
    def sort_key(cls):
        try:
            return EVIDENCE_ORDER.index(cls)
        except ValueError:
            return len(EVIDENCE_ORDER)  # Unknown -> last
    
    return sorted(classifications, key=sort_key)


def parse_comparison_score(response: str) -> int:
    """
    Parse LLM response to extract comparison score (1, 2, 3, or 4).
    
    Args:
        response: LLM response text
        
    Returns:
        int: 1 (identical), 2 (similar), 3 (related), 4 (different), or 4 (default if parsing fails)
    """
    if not response:
        return 4
    
    # First, try to find "Answer: X" or ": X" pattern (most reliable)
    answer_match = re.search(r'(?:answer|reply)[:\s]*([1234])\b', response.lower())
    if answer_match:
        return int(answer_match.group(1))
    
    # Try to find a standalone digit 1, 2, 3, or 4 at word boundaries
    match = re.search(r'\b([1234])\b', response.strip())
    if match:
        return int(match.group(1))
    
    # Default to "different" if parsing fails
    return 4


def compare_disease_pair_with_group(
    da_disease: str, 
    gencc_group_name: str, 
    gencc_group_data: Dict,
    model: str = "claude-haiku-4-5"
) -> Dict:
    """
    Compare a DA disease with a GenCC group using cascade approach:
    1. First compare with group name
    2. If score = 4, compare with individual names in the group
    3. Return best match found
    
    Args:
        da_disease: Disease name from Deep Analysis
        gencc_group_name: Group name from GenCC
        gencc_group_data: Dict with 'names' (list of individual names) and 'classifications'
        model: LLM model to use
        
    Returns:
        dict: Comparison result with score and matched gencc_disease name
    """
    # STEP 1: Compare with group name first
    group_result = compare_disease_pair(da_disease, gencc_group_name, model)
    best_score = group_result["score"]
    matched_gencc_name = gencc_group_name
    
    # STEP 2: If score = 4 (completely different), try individual names
    if best_score == 4:
        individual_names = gencc_group_data.get('names', [])
        for individual_name in individual_names:
            # Skip if same as group name (already tried)
            if individual_name == gencc_group_name:
                continue
            
            # Try token matching first (fast)
            tokens_da = _tokenize_disease_name(da_disease)
            tokens_gencc = _tokenize_disease_name(individual_name)
            
            if tokens_da and tokens_gencc:
                intersection = tokens_da & tokens_gencc
                coverage_da = len(intersection) / len(tokens_da) if tokens_da else 0
                coverage_gencc = len(intersection) / len(tokens_gencc) if tokens_gencc else 0
                max_coverage = max(coverage_da, coverage_gencc)
                
                if max_coverage >= 0.3:  # 30% threshold
                    # Found match via token matching
                    return {
                        "da_disease": da_disease,
                        "gencc_disease": gencc_group_name,  # Keep group name for consistency
                        "gencc_individual_name": individual_name,  # Track which individual name matched
                        "score": 1
                    }
            
            # If token matching didn't work, try LLM comparison
            individual_result = compare_disease_pair(da_disease, individual_name, model)
            individual_score = individual_result["score"]
            
            # Keep best score (lower is better: 1=identical, 2=similar, 3=related, 4=different)
            if individual_score < best_score:
                best_score = individual_score
                matched_gencc_name = gencc_group_name  # Still use group name for consistency
                # If we found a good match (score <= 3), we can stop early
                if best_score <= 3:
                    break
    
    return {
        "da_disease": da_disease,
        "gencc_disease": matched_gencc_name,
        "score": best_score
    }


def compare_disease_pair(da_disease: str, gencc_disease: str, model: str = "claude-haiku-4-5") -> Dict:
    """
    Compare a single pair of diseases using asymmetric coverage first, then LLM if needed.
    
    Args:
        da_disease: Disease name from Deep Analysis
        gencc_disease: Disease name from GenCC
        model: LLM model to use
        
    Returns:
        dict: Comparison result with score
    """
    # OPTIMIZATION: Check asymmetric coverage first
    # If >= 30% of tokens from either name are in the other, consider them identical (score=1)
    tokens_da = _tokenize_disease_name(da_disease)
    tokens_gencc = _tokenize_disease_name(gencc_disease)
    
    if tokens_da and tokens_gencc:
        intersection = tokens_da & tokens_gencc
        coverage_da = len(intersection) / len(tokens_da)    # % of DA tokens in GenCC
        coverage_gencc = len(intersection) / len(tokens_gencc)  # % of GenCC tokens in DA
        max_coverage = max(coverage_da, coverage_gencc)
        
        if max_coverage >= 0.3:  # 30% threshold
            return {
                "da_disease": da_disease,
                "gencc_disease": gencc_disease,
                "score": 1  # Identical based on coverage
            }
    
    # Coverage < 30%: call LLM for precise comparison
    prompt_template = load_prompt("gencc_disease_comparison")
    if not prompt_template:
        # Fallback if prompt not found
        prompt = f'Are "{da_disease}" and "{gencc_disease}" the same? Reply: 1 (identical), 2 (similar), or 3 (different).'
    else:
        prompt = prompt_template.format(
            disease_deep_analysis=da_disease,
            disease_gencc=gencc_disease
        )
    
    try:
        response = call_llm(
            prompt,
            model=model,
            temperature=0.0,
            agent_name="gencc_comparison_agent"
        )
        score = parse_comparison_score(response)
        
        # If score is 4 (completely different), retry with Sonnet for better accuracy
        if score == 4:
            try:
                sonnet_response = call_llm(
                    prompt,
                    model="claude-sonnet-4-5",
                    temperature=0.0,
                    agent_name="gencc_comparison_agent"
                )
                sonnet_score = parse_comparison_score(sonnet_response)
                # Keep Sonnet's response (it's more accurate)
                score = sonnet_score
            except Exception as e:
                print(f"WARNING: Sonnet retry failed for '{da_disease}' vs '{gencc_disease}': {e}")
                # Keep the original score of 4 if Sonnet fails
                pass
    except Exception as e:
        print(f"WARNING: LLM comparison failed for '{da_disease}' vs '{gencc_disease}': {e}")
        score = 3  # Default to "different" on error
    
    return {
        "da_disease": da_disease,
        "gencc_disease": gencc_disease,
        "score": score
    }


def _tokenize_disease_name(name: str) -> set:
    """
    Tokenize a disease name into meaningful words.
    Filters out short tokens, numbers, and strips punctuation.
    """
    # Split on whitespace
    raw_tokens = name.lower().strip().split()
    # Strip punctuation from each token
    cleaned_tokens = [t.strip('.,;:()[]{}"\'-') for t in raw_tokens]
    # Keep tokens with length > 2 and not purely numeric
    return {t for t in cleaned_tokens if len(t) > 2 and not t.isdigit()}


def _jaccard_similarity(tokens1: set, tokens2: set) -> float:
    """
    Calculate Jaccard similarity between two sets of tokens.
    Returns |A ∩ B| / |A ∪ B|
    """
    if not tokens1 or not tokens2:
        return 0.0
    
    intersection = tokens1 & tokens2
    union = tokens1 | tokens2
    
    return len(intersection) / len(union) if union else 0.0


def _compute_common_core(names: List[str]) -> str:
    """
    Compute the common core name from a list of disease names.
    
    Returns the intersection of tokens, ordered by how they appear in the shortest name.
    
    Example:
        ["complex hereditary spastic paraplegia", "hereditary spastic paraplegia 75"]
        → "hereditary spastic paraplegia"
    """
    if not names:
        return ""
    
    if len(names) == 1:
        return names[0].title()
    
    # Tokenize all names
    all_token_sets = [_tokenize_disease_name(name) for name in names]
    
    # Compute intersection
    common_tokens = all_token_sets[0]
    for token_set in all_token_sets[1:]:
        common_tokens = common_tokens & token_set
    
    if not common_tokens:
        # No common tokens, fall back to shortest name
        return min(names, key=len).title()
    
    # Use the shortest name as reference for word order
    reference_name = min(names, key=len)
    reference_words = reference_name.lower().strip().split()
    
    # Build the common core name preserving order from reference
    ordered_tokens = [word for word in reference_words if word in common_tokens or 
                      any(word.startswith(t) or t.startswith(word) for t in common_tokens)]
    
    # If ordered_tokens is empty, just join the common tokens
    if not ordered_tokens:
        ordered_tokens = sorted(common_tokens)
    
    return " ".join(ordered_tokens).title()


def _names_are_similar(name1: str, name2: str, threshold: float = 0.3) -> bool:
    """
    Check if two disease names are similar using asymmetric coverage on words.
    
    Two names are similar if at least `threshold` (30%) of tokens from either name
    are present in the other (whichever direction gives the higher coverage).
    
    Args:
        name1: First disease name
        name2: Second disease name
        threshold: Minimum coverage threshold (default: 0.3 = 30%)
        
    Returns:
        True if names are similar
    """
    n1 = name1.lower().strip()
    n2 = name2.lower().strip()
    
    # Exact match
    if n1 == n2:
        return True
    
    # Substring inclusion (original logic)
    if n1 in n2 or n2 in n1:
        return True
    
    # Asymmetric coverage on tokens
    tokens1 = _tokenize_disease_name(name1)
    tokens2 = _tokenize_disease_name(name2)
    
    if not tokens1 or not tokens2:
        return False
    
    intersection = tokens1 & tokens2
    
    # Coverage = what fraction of one set is in the intersection
    coverage1 = len(intersection) / len(tokens1)  # % of tokens1 in tokens2
    coverage2 = len(intersection) / len(tokens2)  # % of tokens2 in tokens1
    
    # Match if at least threshold% of either set is in common
    return max(coverage1, coverage2) >= threshold


def _group_gencc_diseases(gencc_diseases: List[Dict], coverage_threshold: float = 0.3) -> Dict:
    """
    Group GenCC diseases by name similarity using asymmetric coverage on words.
    Returns a dict mapping group_name -> {names: list, classifications: list}
    
    Two diseases are grouped if:
    1. One name contains the other (substring), OR
    2. At least 30% of tokens from one name are in the other (either direction)
    
    Args:
        gencc_diseases: List of disease dicts with 'disease_title' and classifications
        coverage_threshold: Minimum coverage to group (default: 0.3 = 30%)
    """
    if not gencc_diseases:
        return {}
    
    # First, collect all individual diseases with their classifications
    all_gencc_entries = []
    for d in gencc_diseases:
        disease_name = d.get('disease_title', '')
        if not disease_name:
            continue
        
        # Collect classifications
        if 'all_classifications' in d and d['all_classifications']:
            classifications = list(d['all_classifications'])
        else:
            classifications = [d.get('classification_title', 'N/A')]
        
        all_gencc_entries.append({
            'name': disease_name,
            'classifications': set(classifications)
        })
    
    # Group entries by matching names (transitive closure)
    groups = []
    processed_indices = set()
    
    for i, entry1 in enumerate(all_gencc_entries):
        if i in processed_indices:
            continue
        
        # Start a new group
        group = {
            'names': [entry1['name']],
            'classifications': entry1['classifications'].copy()
        }
        processed_indices.add(i)
        
        # Find all entries that match this one (transitive)
        changed = True
        while changed:
            changed = False
            for j, entry2 in enumerate(all_gencc_entries):
                if j in processed_indices:
                    continue
                
                # Check if entry2's name matches any name in the current group
                for existing_name in group['names']:
                    if _names_are_similar(existing_name, entry2['name'], coverage_threshold):
                        # Match found - add to group
                        group['names'].append(entry2['name'])
                        group['classifications'].update(entry2['classifications'])
                        processed_indices.add(j)
                        changed = True
                        break
        
        groups.append(group)
    
    # Build result dict: group_name (common core) -> group data
    result = {}
    for group in groups:
        # Use common core as group identifier (intersection of tokens)
        group_name = _compute_common_core(group['names'])
        result[group_name] = {
            'names': sorted(group['names']),  # All original names in group
            'classifications': _sort_classifications_by_evidence(list(group['classifications']))  # Sorted by evidence
        }
    
    return result


def _parse_llm_grouping_response(response: str) -> List[Tuple[str, List[str]]]:
    """
    Parse the LLM grouping response to extract groups and their diseases.
    
    Returns:
        List of tuples: [(group_title, [disease1, disease2, ...]), ...]
    """
    groups = []
    current_title = None
    current_diseases = []
    
    for line in response.strip().split('\n'):
        line = line.strip()
        if not line:
            continue
            
        if line.upper().startswith('GROUP:'):
            # Save previous group if exists
            if current_title and current_diseases:
                groups.append((current_title, current_diseases))
            # Start new group
            current_title = line[6:].strip().strip('[]')
            current_diseases = []
        elif line.startswith('-'):
            # Disease entry
            disease = line[1:].strip()
            if disease:
                current_diseases.append(disease)
    
    # Don't forget last group
    if current_title and current_diseases:
        groups.append((current_title, current_diseases))
    
    return groups


def _llm_group_gencc_diseases(
    algo_groups: Dict,
    model: str = "claude-haiku-4-5"
) -> Dict:
    """
    Use LLM to semantically group GenCC diseases and assign meaningful titles.
    
    Takes the output of algorithmic grouping and asks LLM to create semantic groups.
    
    Args:
        algo_groups: Output from _group_gencc_diseases (group_name -> {names, classifications})
        model: LLM model to use
        
    Returns:
        Dict mapping group_title -> {names: list, classifications: list}
    """
    if not algo_groups:
        return {}
    
    # Collect all disease names from algorithmic groups
    all_diseases = []
    disease_to_classifications = {}  # Map disease name -> classifications
    
    for group_name, group_data in algo_groups.items():
        for disease_name in group_data['names']:
            if disease_name not in disease_to_classifications:
                all_diseases.append(disease_name)
                disease_to_classifications[disease_name] = group_data['classifications']
    
    # If only 1 disease, no need for LLM grouping
    if len(all_diseases) <= 1:
        return algo_groups
    
    # Create prompt
    diseases_list = '\n'.join([f'- {d}' for d in all_diseases])
    prompt_template = load_prompt("gencc_disease_grouping")
    if not prompt_template:
        # Fallback if prompt not found
        prompt_template = """Group these diseases by clinical similarity. Output format:
GROUP: [title]
- disease1
- disease2

GROUP: [title]
- disease3

Diseases:
{diseases_list}"""
    prompt = prompt_template.format(diseases_list=diseases_list)
    
    try:
        response = call_llm(
            prompt,
            model=model,
            temperature=0.0,
            agent_name="gencc_grouping_agent"
        )
        
        # Parse response
        parsed_groups = _parse_llm_grouping_response(response)
        
        if not parsed_groups:
            # Parsing failed, return original groups
            if VERBOSE_MODE:
                print("WARNING: LLM grouping parsing failed, using algorithmic groups")
            return algo_groups
        
        # Build result dict
        result = {}
        for title, diseases in parsed_groups:
            # Collect all classifications for diseases in this group
            all_classifications = set()
            matched_disease_names = []
            
            for disease in diseases:
                # Find matching disease (case-insensitive)
                disease_lower = disease.lower()
                for orig_disease in all_diseases:
                    if orig_disease.lower() == disease_lower or disease_lower in orig_disease.lower() or orig_disease.lower() in disease_lower:
                        matched_disease_names.append(orig_disease)
                        if orig_disease in disease_to_classifications:
                            all_classifications.update(disease_to_classifications[orig_disease])
                        break
            
            if matched_disease_names:
                result[title] = {
                    'names': matched_disease_names,
                    'classifications': _sort_classifications_by_evidence(list(all_classifications))
                }
        
        # Handle diseases not assigned by LLM (fallback to original)
        assigned_diseases = set()
        for _, group_data in result.items():
            assigned_diseases.update(d.lower() for d in group_data['names'])
        
        for disease in all_diseases:
            if disease.lower() not in assigned_diseases:
                # Disease was missed by LLM, add as its own group
                result[disease.title()] = {
                    'names': [disease],
                    'classifications': disease_to_classifications.get(disease, [])
                }
        
        if VERBOSE_MODE:
            print(f"INFO: LLM grouped {len(all_diseases)} diseases into {len(result)} groups")
        
        return result
        
    except Exception as e:
        print(f"WARNING: LLM grouping failed: {e}, using algorithmic groups")
        return algo_groups


def compare_diseases_batch(
    da_diseases: List[str],
    gencc_diseases: List[Dict],
    model: str = "claude-haiku-4-5",
    max_workers: int = 10
) -> Dict:
    """
    Compare all pairs of diseases between Deep Analysis and GenCC.
    First groups GenCC diseases, then compares DA diseases with GenCC groups.
    
    Args:
        da_diseases: List of disease names from Deep Analysis
        gencc_diseases: List of dicts with 'disease_title' and 'classification_title' from GenCC
        model: LLM model to use
        max_workers: Maximum number of parallel workers
        
    Returns:
        dict: Categorized comparison results
    """
    if not da_diseases or not gencc_diseases:
        return {
            "matched": [],
            "gencc_only": [],
            "da_only": da_diseases.copy() if da_diseases else [],
            "comparison_matrix": []
        }
    
    # STEP 1: Group GenCC diseases algorithmically first
    if VERBOSE_MODE:
        print(f"INFO: Algorithmic grouping of {len(gencc_diseases)} GenCC diseases...")
    algo_groups = _group_gencc_diseases(gencc_diseases)
    
    # STEP 1b: Use LLM to refine grouping with semantic titles
    if VERBOSE_MODE:
        print(f"INFO: LLM semantic grouping of {len(algo_groups)} algorithmic groups...")
    gencc_groups = _llm_group_gencc_diseases(algo_groups, model=model)
    gencc_group_names = sorted(list(gencc_groups.keys()))
    
    if not gencc_group_names:
        return {
            "matched": [],
            "gencc_only": [],
            "da_only": da_diseases.copy(),
            "comparison_matrix": []
        }
    
    if VERBOSE_MODE:
        print(f"INFO: Created {len(gencc_group_names)} GenCC groups from {len(gencc_diseases)} entries")
    
    # STEP 2: Generate comparison tasks: DA disease vs GenCC group (with cascade to individual names)
    comparison_tasks = []
    for da_disease in da_diseases:
        for gencc_group_name in gencc_group_names:
            gencc_group_data = gencc_groups[gencc_group_name]
            comparison_tasks.append((da_disease, gencc_group_name, gencc_group_data))
    
    if not comparison_tasks:
        return {
            "matched": [],
            "gencc_only": [{
                "name": group_name,
                "classification": groups_data['classifications'][0] if groups_data['classifications'] else "N/A",
                "all_classifications": groups_data['classifications']
            } for group_name, groups_data in gencc_groups.items()],
            "da_only": da_diseases.copy(),
            "comparison_matrix": []
        }
    
    # STEP 3: Execute comparisons in parallel (DA disease vs GenCC group with cascade)
    if VERBOSE_MODE:
        print(f"INFO: Comparing {len(comparison_tasks)} pairs ({len(da_diseases)} DA × {len(gencc_group_names)} GenCC groups) with cascade matching...")
    start_time = perf_counter()
    
    comparison_matrix = []
    workers = min(max_workers, len(comparison_tasks))
    
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(compare_disease_pair_with_group, da, gc_name, gc_data, model): (da, gc_name)
            for da, gc_name, gc_data in comparison_tasks
        }
        
        for future in as_completed(futures):
            try:
                result = future.result()
                comparison_matrix.append(result)
            except Exception as e:
                da, gc = futures[future]
                print(f"WARNING: Comparison task failed for '{da}' vs '{gc}': {e}")
                comparison_matrix.append({
                    "da_disease": da,
                    "gencc_disease": gc,
                    "score": 3  # Default to different on error
                })
    
    elapsed = perf_counter() - start_time
    if VERBOSE_MODE:
        print(f"INFO: Completed {len(comparison_matrix)} comparisons in {elapsed:.2f} seconds")
    
            # STEP 4: Categorize results
    # Group matches by gencc_group_name (already grouped, but may have multiple DA matches)
    matched_by_gencc = {}  # Maps gencc_group_name -> {da_names: set, classifications: list, scores: list}
    matched_da_names = set()
    matched_gencc_group_names = set()
    
    for comp in comparison_matrix:
        if comp["score"] in [1, 2, 3]:  # Identical, similar, or related (not 4=completely different)
            gencc_group_name = comp["gencc_disease"]  # This is actually the group name
            da_name = comp["da_disease"]
            
            # Get all classifications for this group
            group_data = gencc_groups.get(gencc_group_name, {})
            classifications = group_data.get('classifications', ["N/A"])
            
            # Group by gencc_group_name
            if gencc_group_name not in matched_by_gencc:
                matched_by_gencc[gencc_group_name] = {
                    "da_names": set(),
                    "classifications": classifications,
                    "scores": []
                }
            
            matched_by_gencc[gencc_group_name]["da_names"].add(da_name)
            matched_by_gencc[gencc_group_name]["scores"].append(comp["score"])
            matched_da_names.add(da_name)
            matched_gencc_group_names.add(gencc_group_name)
    
    # Convert grouped matches to list format
    matched = []
    for gencc_group_name, match_data in matched_by_gencc.items():
        # Use the first DA name as primary (sorted for consistency)
        da_names_list = sorted(list(match_data["da_names"]))
        primary_da_name = da_names_list[0]
        
        # Get best score (1 is better than 2)
        best_score = min(match_data["scores"])
        
        classifications = match_data["classifications"]
        
        # Get original disease names in this group
        group_data = gencc_groups.get(gencc_group_name, {})
        diseases_in_group = group_data.get('names', [])
        
        matched.append({
            "da_name": primary_da_name,
            "da_names": da_names_list,  # Include all matching DA names
            "gencc_name": gencc_group_name,  # Group name (LLM title or common core)
            "gencc_diseases_in_group": diseases_in_group,  # Original GenCC disease names
            "gencc_classification": classifications[0] if classifications else "N/A",  # Primary classification
            "gencc_all_classifications": classifications,  # All classifications
            "score": best_score  # Store best score (1, 2, or 3)
        })
    
    # Diseases only in GenCC (not matched)
    gencc_only = []
    for group_name in gencc_group_names:
        if group_name not in matched_gencc_group_names:
            group_data = gencc_groups[group_name]
            gencc_only.append({
                "name": group_name,
                "diseases_in_group": group_data.get('names', []),  # Original GenCC disease names
                "classification": group_data['classifications'][0] if group_data['classifications'] else "N/A",  # Primary classification
                "all_classifications": group_data['classifications']  # All classifications
            })
    
    # Diseases only in Deep Analysis (not matched with score 1 or 2)
    da_only = [
        {"name": name}
        for name in da_diseases
        if name not in matched_da_names
    ]
    
    return {
        "matched": matched,
        "gencc_only": gencc_only,
        "da_only": da_only,
        "comparison_matrix": comparison_matrix  # Full matrix for debugging/analysis
    }

