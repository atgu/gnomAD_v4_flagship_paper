"""Parsers for LLM responses."""
import re
import xml.etree.ElementTree as ET


INVALID_DISEASE_NAMES = {
    "disease",
    "disease name",
    "conflicting evidence summary",
    "overall summary",
    "summary",
    # Legacy placeholders (for backward compatibility)
    "no diseases identified",
    "no diseases qualify",
    "no diseases detected",
    "no diseases meet criteria",
    "no data available",
    "none identified",
    "no qualifying diseases",
    "no diseases meet inclusion criteria",
    "*no diseases identified*",
    "**no diseases identified**",
    "*no diseases qualify*",
    "*no diseases detected*",
    "*no diseases meet criteria*",
    "*no data available*",
    "*none identified*",
    "*no qualifying diseases*",
}


INHERITANCE_SCORE_LABELS = {
    1: "dominance",
    2: "incomplete dominance mostly dominant",
    3: "incomplete dominance",
    4: "co-dominant",
    5: "incomplete dominance mostly recessive",
    6: "recessive",
    7: "unknown-inheritance",
    8: "conflicting evidence",
}

INHERITANCE_LABEL_TO_SCORE = {
    label: score for score, label in INHERITANCE_SCORE_LABELS.items()
}


def _parse_markdown_table(response_text, min_columns):
    """
    Parses a Markdown table response into a list of rows (each row is a list of column strings).
    """
    rows = []
    if not response_text:
        return rows

    for line in response_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if set(stripped) <= set("|-: "):
            continue

        parts = [col.strip() for col in stripped.split("|")]
        parts = [p for p in parts if p != ""]
        if len(parts) < min_columns:
            continue
        rows.append(parts[:min_columns])

    return rows


def parse_llm_response(response_text, articles_retrieved):
    """Parses the XML-like response from the LLM."""
    if not response_text:
        return None
    
    try:
        level_match = re.search(r'<LEVEL>(.*?)</LEVEL>', response_text, re.DOTALL)
        just_match = re.search(r'<JUSTIFICATION>(.*?)</JUSTIFICATION>', response_text, re.DOTALL)
        pmids_match = re.search(r'<PMIDS>(.*?)</PMIDS>', response_text, re.DOTALL)

        level = int(level_match.group(1).strip()) if level_match else "N/A"
        justification = just_match.group(1).strip() if just_match else "No justification provided."
        
        pmids_str = pmids_match.group(1).strip() if pmids_match else ""
        pmids = []
        if pmids_str and "no article" not in pmids_str.lower():
            raw_pmids = [p.strip() for p in pmids_str.split(',') if p.strip()]
            for p_raw in raw_pmids:
                match = re.search(r'\d+', p_raw)
                if match:
                    pmids.append(match.group(0))

        return {"level": level, "justification": justification, "pmids": pmids}

    except (AttributeError, ValueError) as e:
        print(f"ERROR: Could not parse LLM response: {e}\nResponse was:\n{response_text}")
        return None


def parse_agreement_analysis_response(response_text):
    """Returns the raw agreement analysis response (no XML parsing needed)."""
    if not response_text:
        return None
    
    return response_text.strip()


def parse_phenotype_synonyms_response(response_text):
    """Parses the XML response from the phenotype synonyms LLM."""
    if not response_text:
        return []
    
    try:
        # Try to extract synonyms using regex
        synonyms_match = re.search(r'<synonyms>(.*?)</synonyms>', response_text, re.DOTALL | re.IGNORECASE)
        
        if not synonyms_match:
            print("WARNING: Could not find <synonyms> tags in response")
            return []
        
        synonyms_content = synonyms_match.group(1).strip()
        
        if not synonyms_content:
            print("INFO: No synonyms returned (empty content)")
            return []
        
        # Extract individual synonym tags
        synonym_matches = re.findall(r'<synonym>(.*?)</synonym>', synonyms_content, re.DOTALL | re.IGNORECASE)
        
        synonyms = [s.strip() for s in synonym_matches if s.strip()]
        
        print(f"INFO: Parsed {len(synonyms)} phenotype synonyms")
        return synonyms
        
    except Exception as e:
        print(f"ERROR: Failed to parse phenotype synonyms response: {e}")
        return []


def _extract_pmids(pmids_text, limit=None):
    if not pmids_text or pmids_text.lower() == "none":
        return []
    pmids = []
    for token in pmids_text.split(","):
        token = token.strip()
        if not token:
            continue
        match = re.search(r'\d+', token)
        if match:
            pmids.append(match.group(0))
        else:
            pmids.append(token)
    if limit:
        return pmids[:limit]
    return pmids


def parse_deep_analysis_a1(response_text):
    rows = _parse_markdown_table(response_text, 4)
    diseases = []
    for row in rows:
        disease, score, justification, pmids_text = row
        disease = disease.strip()
        disease_lower = disease.lower()
        
        # Skip invalid/placeholder names (but keep "None" with justification)
        if not disease or disease_lower in INVALID_DISEASE_NAMES:
            continue
        
        # Handle "None" case: keep it with justification for "no phenotype" detection
        if disease_lower == "none":
            diseases.append({
                "disease": "None",
                "association_score": None,
                "justification": justification.strip() if justification else None,
                "pmids": [],
            })
            continue
        
        try:
            score_val = int(score)
        except ValueError:
            score_val = None
        diseases.append(
            {
                "disease": disease,
                "association_score": score_val,
                "justification": justification,
                "pmids": _extract_pmids(pmids_text, limit=5),
            }
        )
    return diseases


def parse_deep_analysis_a2(response_text):
    rows = _parse_markdown_table(response_text, 4)
    results = []
    for row in rows:
        disease, level, justification, pmids_text = row
        disease = disease.strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue
        try:
            level_val = int(level)
        except ValueError:
            level_val = None
        results.append(
            {
                "disease": disease,
                "penetrance_level": level_val,
                "justification": justification,
                "pmids": _extract_pmids(pmids_text, limit=5),
            }
        )
    return results


def parse_deep_analysis_a3(response_text):
    rows = _parse_markdown_table(response_text, 4)
    results = []
    for row in rows:
        disease, inheritance_value, justification, pmids_text = row
        disease = disease.strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue

        raw_value = (inheritance_value or "").strip().lower()
        inheritance_score = None
        inheritance_label = None

        if raw_value:
            try:
                inheritance_score = int(raw_value)
                inheritance_label = INHERITANCE_SCORE_LABELS.get(
                    inheritance_score, "unknown-inheritance"
                )
            except ValueError:
                inheritance_label = raw_value
                inheritance_score = INHERITANCE_LABEL_TO_SCORE.get(
                    raw_value, None
                )

        if not inheritance_label and inheritance_score:
            inheritance_label = INHERITANCE_SCORE_LABELS.get(
                inheritance_score, "unknown-inheritance"
            )
        if not inheritance_label:
            inheritance_label = raw_value or "unknown-inheritance"

        results.append(
            {
                "disease": disease,
                "inheritance": inheritance_label,
                "inheritance_score": inheritance_score,
                "justification": justification,
                "pmids": _extract_pmids(pmids_text, limit=5),
            }
        )
    return results


def parse_deep_analysis_a4(response_text):
    rows = _parse_markdown_table(response_text, 5)
    results = []
    for row in rows:
        disease, mechanism, confidence, justification, pmids_text = row
        disease = disease.strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue
        try:
            confidence_val = int(confidence)
        except ValueError:
            confidence_val = None
        results.append(
            {
                "disease": disease,
                "mechanism": mechanism,
                "confidence_score": confidence_val,
                "justification": justification,
                "pmids": _extract_pmids(pmids_text, limit=3),
            }
        )
    return results


def parse_onset_severity_response(response_text):
    rows = _parse_markdown_table(response_text, 4)
    results = {}
    for row in rows:
        disease, onset_score_text, severity_text, justification = row[:4]
        disease = disease.strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue

        try:
            onset_score = int(onset_score_text.strip())
        except ValueError:
            onset_score = None

        try:
            severity_val = int(severity_text.strip())
        except ValueError:
            severity_val = None

        results[disease] = {
            "onset_score": onset_score,
            "severity_score": severity_val,
            "justification": justification.strip(),
        }
    return results


def parse_pathogenicity_response(response_text):
    rows = _parse_markdown_table(response_text, 2)
    mapping = {}
    for row in rows:
        disease, classification = row[:2]
        disease = disease.strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue
        cls_clean = (classification or "").strip().lower()
        if cls_clean not in {"pathogenic", "protective", "neutral"}:
            cls_clean = "neutral"
        mapping[disease] = cls_clean
    return mapping




def parse_algorithmic_summary_response(response_text):
    if not response_text:
        return None
    try:
        summary_match = re.search(r"<SUMMARY>(.*?)</SUMMARY>", response_text, re.DOTALL | re.IGNORECASE)
        summary_text = summary_match.group(1).strip() if summary_match else response_text.strip()
        summary_text = re.sub(r"\s+", " ", summary_text)
        return summary_text or None
    except Exception as exc:
        print(f"ERROR: Failed to parse algorithmic summary response: {exc}")
        return None


def parse_novelty_assessment_response(response_text, all_available_pmids):
    """Parses the XML response from the novelty assessment LLM."""
    if not response_text:
        return None
    
    try:
        # Extract verdict
        verdict_match = re.search(r'<verdict>(.*?)</verdict>', response_text, re.DOTALL | re.IGNORECASE)
        verdict = verdict_match.group(1).strip() if verdict_match else "Unknown"
        
        # Validate verdict
        valid_verdicts = ["Novel", "Hypothesized", "Existing", "Established"]
        if verdict not in valid_verdicts:
            print(f"WARNING: Invalid verdict '{verdict}', defaulting to 'Unknown'")
            verdict = "Unknown"
        
        # Extract justification
        just_match = re.search(r'<justification>(.*?)</justification>', response_text, re.DOTALL | re.IGNORECASE)
        justification = just_match.group(1).strip() if just_match else "No justification provided."
        
        # Extract PMIDs
        pmids_match = re.search(r'<pmids>(.*?)</pmids>', response_text, re.DOTALL | re.IGNORECASE)
        pmids = []
        
        if pmids_match:
            pmids_content = pmids_match.group(1).strip()
            # Extract individual PMID tags
            pmid_matches = re.findall(r'<pmid>(.*?)</pmid>', pmids_content, re.DOTALL | re.IGNORECASE)
            
            for pmid in pmid_matches:
                pmid = pmid.strip()
                # Extract digits only
                digit_match = re.search(r'\d+', pmid)
                if digit_match:
                    extracted_pmid = digit_match.group(0)
                    # Validate against available PMIDs
                    if extracted_pmid in all_available_pmids:
                        pmids.append(extracted_pmid)
                    else:
                        print(f"WARNING: PMID {extracted_pmid} not in available PMIDs list")
        
        # Limit to 10 PMIDs
        pmids = pmids[:10]
        
        print(f"INFO: Parsed novelty assessment - Verdict: {verdict}, PMIDs: {len(pmids)}")
        
        return {
            "verdict": verdict,
            "justification": justification,
            "pmids": pmids
        }
        
    except Exception as e:
        print(f"ERROR: Failed to parse novelty assessment response: {e}")
        return None


# ---------------------------------------------------------------------------
# Probabilistic mode parsers
# ---------------------------------------------------------------------------

PENETRANCE_PROBA_KEYS = ["mendelian", "high", "moderate", "complex"]
INHERITANCE_PROBA_KEYS = ["dominant", "inc_dom", "incomplete", "codominant", "inc_rec", "recessive"]
ONSET_PROBA_KEYS = ["prenatal", "neonatal", "infancy", "childhood", "adolescence", "adulthood", "late"]
SEVERITY_PROBA_KEYS = ["lethal", "severe", "moderate", "mild", "verymild"]


def _parse_proba_value(val_str):
    """Parse a probability value (0-100) from string, return 0 if invalid."""
    try:
        val = int(val_str.strip())
        if 0 <= val <= 100:
            return val
        return 0
    except (ValueError, AttributeError):
        return 0


def _normalize_distribution(dist_dict):
    """Normalize a distribution dict so probabilities sum to 100."""
    total = sum(dist_dict.values())
    if total == 0:
        return dist_dict
    if total == 100:
        return dist_dict
    # Normalize
    return {k: round(v * 100 / total) for k, v in dist_dict.items()}


def parse_penetrance_proba_response(response_text):
    """
    Parse probabilistic penetrance response.
    
    Expected format:
    | Disease | P_mendelian | P_high | P_moderate | P_complex | Justification | PMIDs |
    
    Returns:
        Dict mapping disease names to distribution dicts with justification and pmids
        e.g., {"Huntington": {"distribution": {...}, "justification": "...", "pmids": [...]}}
    """
    rows = _parse_markdown_table(response_text, 7)
    results = {}
    
    for row in rows:
        disease = row[0].strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue
        
        dist = {
            "mendelian": _parse_proba_value(row[1]),
            "high": _parse_proba_value(row[2]),
            "moderate": _parse_proba_value(row[3]),
            "complex": _parse_proba_value(row[4]),
        }
        
        # Normalize to ensure sum is 100
        dist = _normalize_distribution(dist)
        
        justification = row[5].strip() if len(row) > 5 else None
        pmids = _extract_pmids(row[6], limit=5) if len(row) > 6 else []
        
        results[disease] = {
            "distribution": dist,
            "justification": justification,
            "pmids": pmids,
        }
    
    return results


def parse_inheritance_proba_response(response_text):
    """
    Parse probabilistic inheritance response.
    
    Expected format:
    | Disease | P_dominant | P_inc_dom | P_incomplete | P_codominant | P_inc_rec | P_recessive | Justification | PMIDs |
    
    Returns:
        Dict mapping disease names to distribution dicts with justification and pmids
        e.g., {"Huntington": {"distribution": {...}, "justification": "...", "pmids": [...]}}
    """
    rows = _parse_markdown_table(response_text, 9)
    results = {}
    
    for row in rows:
        disease = row[0].strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue
        
        dist = {
            "dominant": _parse_proba_value(row[1]),
            "inc_dom": _parse_proba_value(row[2]),
            "incomplete": _parse_proba_value(row[3]),
            "codominant": _parse_proba_value(row[4]),
            "inc_rec": _parse_proba_value(row[5]),
            "recessive": _parse_proba_value(row[6]),
        }
        
        # Normalize to ensure sum is 100
        dist = _normalize_distribution(dist)
        
        justification = row[7].strip() if len(row) > 7 else None
        pmids = _extract_pmids(row[8], limit=5) if len(row) > 8 else []
        
        results[disease] = {
            "distribution": dist,
            "justification": justification,
            "pmids": pmids,
        }
    
    return results


def parse_onset_severity_proba_response(response_text):
    """
    Parse probabilistic onset/severity response.
    
    Expected format:
    | Disease | P_prenatal | P_neonatal | P_infancy | P_childhood | P_adolescence | P_adulthood | P_late | P_lethal | P_severe | P_moderate | P_mild | P_verymild | Justification |
    
    Returns:
        Dict mapping disease names to dicts with "onset" and "severity" distributions plus justification
        e.g., {"Huntington": {"onset": {...}, "severity": {...}, "justification": "..."}}
    """
    rows = _parse_markdown_table(response_text, 14)
    results = {}
    
    for row in rows:
        disease = row[0].strip()
        if not disease or disease.lower() in INVALID_DISEASE_NAMES:
            continue
        
        onset_dist = {
            "prenatal": _parse_proba_value(row[1]),
            "neonatal": _parse_proba_value(row[2]),
            "infancy": _parse_proba_value(row[3]),
            "childhood": _parse_proba_value(row[4]),
            "adolescence": _parse_proba_value(row[5]),
            "adulthood": _parse_proba_value(row[6]),
            "late": _parse_proba_value(row[7]),
        }
        
        severity_dist = {
            "lethal": _parse_proba_value(row[8]),
            "severe": _parse_proba_value(row[9]),
            "moderate": _parse_proba_value(row[10]),
            "mild": _parse_proba_value(row[11]),
            "verymild": _parse_proba_value(row[12]),
        }
        
        # Normalize each distribution
        onset_dist = _normalize_distribution(onset_dist)
        severity_dist = _normalize_distribution(severity_dist)
        
        justification = row[13].strip() if len(row) > 13 else None
        
        results[disease] = {
            "onset": onset_dist,
            "severity": severity_dist,
            "justification": justification,
        }
    
    return results


# ---------------------------------------------------------------------------
# Phenotyping parsers
# ---------------------------------------------------------------------------

def parse_contamination_response(response_text, name_en_list, endpoint_name="unknown"):
    """
    Parses the contamination check response from the phenotyping LLM.
    More flexible parsing based on endpoints/process_endpoint.py
    
    Expected format:
        PHENOTYPE: [number]. [name_en]
        SCORE: [0-3]
        JUSTIFICATION: [text]
    """
    if not response_text:
        return [{'score': 'NA', 'justification': 'Empty response'} for _ in name_en_list]
    
    try:
        results = []
        lines = response_text.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Look for numbered phenotype (e.g., "1. clozapine; systemic" or "PHENOTYPE: 1. clozapine")
            if (line and (line[0].isdigit() and '. ' in line)) or line.startswith('PHENOTYPE:'):
                current_score = None
                current_justification = None
                
                # Extract phenotype info
                if line.startswith('PHENOTYPE:'):
                    current_phenotype = line.replace('PHENOTYPE:', '').strip()
                else:
                    current_phenotype = line
                
                # Look for SCORE and JUSTIFICATION in the next few lines
                j = i + 1
                while j < len(lines) and j < i + 5:  # Look ahead max 5 lines
                    next_line = lines[j].strip()
                    
                    if next_line.startswith('SCORE:'):
                        score_text = next_line.replace('SCORE:', '').strip()
                        try:
                            score = int(score_text)
                            if score in [0, 1, 2, 3]:
                                current_score = score
                            else:
                                current_score = 'NA'
                        except ValueError:
                            current_score = 'NA'
                    
                    elif next_line.startswith('JUSTIFICATION:'):
                        current_justification = next_line.replace('JUSTIFICATION:', '').strip()
                        break  # Found justification, move to next phenotype
                    
                    j += 1
                
                # Add result if we found both score and justification
                if current_score is not None and current_justification is not None:
                    results.append({'score': current_score, 'justification': current_justification})
                else:
                    results.append({'score': 'NA', 'justification': 'Incomplete parsing - missing score or justification'})
                
                i = j  # Skip to after justification
            
            i += 1
        
        # Ensure we have results for all input phenotypes
        while len(results) < len(name_en_list):
            results.append({'score': 'NA', 'justification': 'LLM response parsing failed'})
        
        # Truncate if we somehow got too many results
        results = results[:len(name_en_list)]
        
        return results
        
    except Exception as e:
        print(f"ERROR: Failed to parse contamination response: {e}")
        return [{'score': 'NA', 'justification': f'Parse error: {str(e)}'} for _ in name_en_list]


