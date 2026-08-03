"""Service for generating agreement analysis between LLM and LOEUF scores.

New pipeline (v2):
1. If agreement >= 0.7: canned response
2. If agreement < 0.7: calculate explanations based on direction
3. If no explanations found: canned response based on direction
4. If explanations found: call LLM with explanations
"""

import os
from typing import Dict, List, Optional, Tuple
from utils.helpers import load_prompt
from services.llm_service import call_llm
from parsers.llm_parsers import parse_agreement_analysis_response
from config import APP_ROOT

# Paths to explanation data files
DATA_DIR = os.path.join(APP_ROOT, "data")
SYN_DATA_PATH = os.path.join(DATA_DIR, "julia_syn.tsv")
TISSUE_SPEC_PATH = os.path.join(DATA_DIR, "gtex_tissue_specificity.tsv")
MOUSE_LETHAL_PATH = os.path.join(DATA_DIR, "embryonic_lethal_mouse_genes_with_human.csv")

# Cached data (loaded once)
_syn_data = None
_tissue_spec_data = None
_mouse_lethal_genes = None


def _load_syn_data():
    """Load synonymous variant data."""
    global _syn_data
    if _syn_data is not None:
        return _syn_data
    
    if not os.path.exists(SYN_DATA_PATH):
        print(f"WARNING: Syn data not found at {SYN_DATA_PATH}")
        _syn_data = {}
        return _syn_data
    
    import pandas as pd
    try:
        df = pd.read_csv(SYN_DATA_PATH, sep='\t')
        _syn_data = {row['gene']: {'obs_syn': row['obs_syn'], 'exp_syn': row['exp_syn']} 
                     for _, row in df.iterrows()}
        print(f"INFO: Loaded syn data for {len(_syn_data)} genes")
    except Exception as e:
        print(f"ERROR loading syn data: {e}")
        _syn_data = {}
    return _syn_data


def _load_tissue_spec_data():
    """Load tissue specificity data."""
    global _tissue_spec_data
    if _tissue_spec_data is not None:
        return _tissue_spec_data
    
    if not os.path.exists(TISSUE_SPEC_PATH):
        print(f"WARNING: Tissue spec data not found at {TISSUE_SPEC_PATH}")
        _tissue_spec_data = {}
        return _tissue_spec_data
    
    import pandas as pd
    try:
        df = pd.read_csv(TISSUE_SPEC_PATH, sep='\t')
        _tissue_spec_data = {}
        for _, row in df.iterrows():
            gene = row['gene']
            _tissue_spec_data[gene] = {
                'testis_rank': row.get('testis_rank'),
                'ovary_rank': row.get('ovary_rank')
            }
        print(f"INFO: Loaded tissue spec data for {len(_tissue_spec_data)} genes")
    except Exception as e:
        print(f"ERROR loading tissue spec data: {e}")
        _tissue_spec_data = {}
    return _tissue_spec_data


def _load_mouse_lethal_genes():
    """Load mouse embryonic lethal genes."""
    global _mouse_lethal_genes
    if _mouse_lethal_genes is not None:
        return _mouse_lethal_genes
    
    if not os.path.exists(MOUSE_LETHAL_PATH):
        print(f"WARNING: Mouse lethal data not found at {MOUSE_LETHAL_PATH}")
        _mouse_lethal_genes = set()
        return _mouse_lethal_genes
    
    import pandas as pd
    try:
        df = pd.read_csv(MOUSE_LETHAL_PATH)
        _mouse_lethal_genes = set(df['HumanSymbol'].dropna().unique())
        print(f"INFO: Loaded {len(_mouse_lethal_genes)} mouse embryonic lethal genes")
    except Exception as e:
        print(f"ERROR loading mouse lethal data: {e}")
        _mouse_lethal_genes = set()
    return _mouse_lethal_genes


def _gamma_upper_bound(obs, exp, quantile=0.95):
    """Calculate upper bound of obs/exp using gamma distribution."""
    from scipy.stats import gamma
    if exp is None or exp <= 0 or obs is None:
        return None
    return gamma.ppf(quantile, obs + 1) / exp


def _gamma_lower_bound(obs, exp, quantile=0.05):
    """Calculate lower bound of obs/exp using gamma distribution."""
    from scipy.stats import gamma
    if exp is None or exp <= 0 or obs is None or obs <= 0:
        return None
    return gamma.ppf(quantile, obs) / exp


def _get_agreement_label(agreement_score: float) -> str:
    """Convert agreement score to label."""
    if agreement_score >= 0.7:
        return "Strong agreement"
    elif agreement_score >= 0.5:
        return "Mild disagreement"
    elif agreement_score >= 0.4:
        return "Moderate disagreement"
    elif agreement_score >= 0.2:
        return "Strong disagreement"
    elif agreement_score >= 0.05:
        return "Very strong disagreement"
    else:
        return "Extreme disagreement"


def _parse_syn_raw_value(raw_value: str) -> Tuple[float, float, float, str]:
    """
    Parse the raw_value from Synonymous Depletion explanation.
    
    Returns:
        (lof_ratio, syn_ratio, fold_diff, bias_level)
        bias_level: "explains_all", "explains_much", "explains_partial", "explains_little", "explains_none"
    """
    import re
    
    # Extract ratios: "LoF obs/exp = 0.25, Syn obs/exp = 0.74 — ..."
    lof_match = re.search(r'LoF obs/exp = ([\d.]+)', raw_value)
    syn_match = re.search(r'Syn obs/exp = ([\d.]+)', raw_value)
    
    lof_ratio = float(lof_match.group(1)) if lof_match else None
    syn_ratio = float(syn_match.group(1)) if syn_match else None
    
    # Calculate fold difference
    fold_diff = None
    if lof_ratio and syn_ratio and lof_ratio > 0:
        R = lof_ratio / syn_ratio
        fold_diff = 1 / R if R > 0 else 999
    
    # Determine bias level from raw_value content
    if "similar depletion" in raw_value or "likely explains" in raw_value:
        bias_level = "explains_all"
    elif "explains much" in raw_value:
        bias_level = "explains_much"
    elif "partially explains" in raw_value:
        bias_level = "explains_partial"
    elif "only small part" in raw_value:
        bias_level = "explains_little"
    else:
        bias_level = "explains_none"
    
    return lof_ratio, syn_ratio, fold_diff, bias_level


def _generate_syn_explanation_text(lof_ratio: float, syn_ratio: float, fold_diff: float, bias_level: str) -> Tuple[str, bool]:
    """
    Generate explanation text for Synonymous Depletion.
    
    Returns:
        (explanation_text, bias_explains_all)
        bias_explains_all: True if we should stop here (no need for other explanations)
    """
    if bias_level == "explains_all":
        # R >= 0.8: Syn explains most → STOP
        text = (
            f"Synonymous variants show similar depletion to LoF variants "
            f"(Syn obs/exp = {syn_ratio:.2f}, LoF obs/exp = {lof_ratio:.2f}), "
            f"suggesting a regional mutation model bias. "
            f"This bias likely explains most of the observed LOEUF constraint; "
            f"the disagreement with loss-of-function disease literature is therefore likely technical in origin."
        )
        return text, True
    
    elif bias_level == "explains_much":
        # R 0.6-0.8: Syn explains much but not all
        text = (
            f"Synonymous variants show depletion (Syn obs/exp = {syn_ratio:.2f}), "
            f"suggesting a regional mutation model bias. "
            f"However, LoF variants are {fold_diff:.1f}x more depleted (LoF obs/exp = {lof_ratio:.2f}), "
            f"indicating that the bias does not fully explain the constraint."
        )
        return text, False
    
    elif bias_level == "explains_partial":
        # R 0.4-0.6: Syn explains partially
        text = (
            f"Synonymous variants show depletion (Syn obs/exp = {syn_ratio:.2f}), "
            f"suggesting a regional mutation model bias. "
            f"However, LoF variants are {fold_diff:.1f}x more depleted (LoF obs/exp = {lof_ratio:.2f}), "
            f"indicating that functional constraint contributes beyond the mutation bias."
        )
        return text, False
    
    elif bias_level == "explains_little":
        # R 0.2-0.4: Syn explains only small part
        text = (
            f"Synonymous variants show depletion (Syn obs/exp = {syn_ratio:.2f}), "
            f"which could suggest a regional mutation model bias. "
            f"However, LoF variants are {fold_diff:.1f}x more depleted (LoF obs/exp = {lof_ratio:.2f}), "
            f"indicating that this bias explains only a small part of the constraint."
        )
        return text, False
    
    else:
        # R < 0.2: Syn doesn't explain
        text = (
            f"Synonymous variants show slight depletion (Syn obs/exp = {syn_ratio:.2f}), "
            f"suggesting a possible mutation model bias. "
            f"However, LoF variants are {fold_diff:.1f}x more depleted (LoF obs/exp = {lof_ratio:.2f}), "
            f"ruling out the bias as an explanation for the strong LOEUF constraint."
        )
        return text, False


# Explanation phrase templates for algorithmic generation
# Note: "Synonymous Depletion" is handled specially with cascade logic

# Hypothesis templates (used when syn doesn't explain all)
HYPOTHESIS_TEMPLATES_LOEUF_SEVERE = {
    "Reproductive Specificity (Testis)": (
        "One hypothesis is that the constraint reflects selection for male fertility: "
        "this gene is highly expressed in testis ({raw_value}), "
        "a phenotype underrepresented in disease genetics studies."
    ),
    "Reproductive Specificity (Ovary)": (
        "One hypothesis is that the constraint reflects selection for female fertility: "
        "this gene is highly expressed in ovary ({raw_value}), "
        "a phenotype underrepresented in disease genetics studies."
    ),
    "Mouse Embryonic Lethal": (
        "One hypothesis is that this gene is essential for early embryonic development: "
        "the mouse ortholog causes embryonic lethality, "
        "a phenotype rarely captured in human loss-of-function disease literature."
    ),
}

# Understudied gene template (used when no other explanation found)
UNDERSTUDIED_TEMPLATE = (
    "This gene is likely understudied in the medical literature. "
    "The strong LOEUF constraint suggests deleterious effects not yet characterized by research, "
    "possibly related to fertility or early embryonic development—phenotypes difficult to observe clinically."
)

# Legacy templates (kept for backward compatibility)
EXPLANATION_TEMPLATES_LOEUF_SEVERE = {
    "Reproductive Specificity (Testis)": "The gene shows high testis-specific expression ({raw_value}), suggesting purifying selection for male fertility not captured in loss-of-function disease literature.",
    "Reproductive Specificity (Ovary)": "The gene shows high ovary-specific expression ({raw_value}), suggesting purifying selection for female fertility not captured in loss-of-function disease literature.",
    "Mouse Embryonic Lethal": "The mouse ortholog causes embryonic lethality, suggesting {gene_name} may be essential for early human development.",
}

EXPLANATION_TEMPLATES_LITERATURE_SEVERE = {
    "Gain-of-Function Mechanism": "The primary mechanism is gain-of-function ({raw_value}), which is not captured by LOEUF as it measures only loss-of-function constraint.",
    "Dominant-Negative Mechanism": "The primary mechanism is dominant-negative ({raw_value}), which is not captured by LOEUF as it measures only loss-of-function constraint.",
    "Late Onset": "The disease has late onset ({raw_value}), reducing selective pressure before reproductive age and thus weakening LOEUF signal.",
    "Recessive Inheritance": "The inheritance is recessive ({raw_value}), meaning heterozygotes are unaffected and LOEUF (which detects heterozygote selection) shows weaker constraint.",
    "Incomplete Penetrance": "The disease shows incomplete penetrance ({raw_value}), reducing selective pressure and weakening the LOEUF signal.",
    "Low Statistical Power": "Limited statistical power ({raw_value}) may make LOEUF unreliable for this gene.",
    "Synonymous Enrichment": "Synonymous variants are enriched ({raw_value}), suggesting potential technical bias inflating the LOEUF score.",
}


def _generate_algorithmic_explanation(
    gene_name: str,
    agreement_label: str,
    direction: str,
    primary: Dict,
    secondary: Dict = None,
    all_explanations: List[Dict] = None
) -> str:
    """
    Generate explanation algorithmically without LLM.
    
    New cascade logic for LOEUF more severe:
    1. Check Synonymous Depletion first
       - If bias explains all → STOP
       - If bias doesn't explain all → continue
    2. Look for other explanations (Repro, Mouse) as hypotheses
    3. If none found → understudied gene message
    
    Args:
        gene_name: Gene symbol
        agreement_label: e.g., "Moderate disagreement"
        direction: "loeuf_more_severe" or "literature_more_severe"
        primary: Primary explanation dict with 'name' and 'raw_value'
        secondary: Optional secondary explanation dict
        all_explanations: Full list of valid explanations (for cascade logic)
    
    Returns:
        Generated explanation text
    """
    # Build explanation parts
    parts = []
    
    # Header: agreement label + direction
    if direction == "loeuf_more_severe":
        direction_text = f"LOEUF indicates stronger constraint than the loss-of-function disease literature predicts for {gene_name}."
    else:
        direction_text = f"Disease literature suggests {gene_name} is more important for health than LOEUF constraint indicates."
    
    parts.append(f"{agreement_label}.")
    parts.append(direction_text)
    
    # === LOEUF MORE SEVERE: Use new cascade logic ===
    if direction == "loeuf_more_severe":
        all_expl = all_explanations or []
        
        # Step 1: Find Synonymous Depletion
        syn_explanation = next((e for e in all_expl if e['name'] == 'Synonymous Depletion'), None)
        
        # Step 2: Find other explanations (Repro, Mouse)
        other_explanations = [e for e in all_expl if e['name'] in [
            'Reproductive Specificity (Testis)',
            'Reproductive Specificity (Ovary)',
            'Mouse Embryonic Lethal'
        ]]
        
        if syn_explanation:
            # Parse syn data
            raw_value = syn_explanation.get('raw_value', '')
            lof_ratio, syn_ratio, fold_diff, bias_level = _parse_syn_raw_value(raw_value)
            
            if lof_ratio is not None and syn_ratio is not None:
                # Generate syn explanation
                syn_text, bias_explains_all = _generate_syn_explanation_text(
                    lof_ratio, syn_ratio, fold_diff, bias_level
                )
                parts.append(syn_text)
                
                if bias_explains_all:
                    # Bias explains everything → STOP here
                    return " ".join(parts)
                
                # Bias doesn't explain all → look for hypotheses
                if other_explanations:
                    hypotheses_added = []
                    for other in other_explanations:
                        template = HYPOTHESIS_TEMPLATES_LOEUF_SEVERE.get(other['name'])
                        if template:
                            hypothesis = template.format(
                                raw_value=other.get('raw_value', ''),
                                gene_name=gene_name
                            )
                            hypotheses_added.append(hypothesis)
                    
                    # Add all hypotheses with proper connectors
                    for i, hyp in enumerate(hypotheses_added):
                        if i == 0:
                            parts.append(hyp)
                        else:
                            # Replace "One hypothesis" with "Another hypothesis" for subsequent
                            hyp_modified = hyp.replace("One hypothesis", "Another hypothesis")
                            parts.append(hyp_modified)
                else:
                    # No other explanation found → understudied
                    parts.append(UNDERSTUDIED_TEMPLATE)
            else:
                # Couldn't parse syn data, fallback to old logic
                parts.append(f"Synonymous variants show depletion ({raw_value}), suggesting possible mutation bias.")
        
        elif other_explanations:
            # No syn data, but have other explanations
            hypotheses_added = []
            for other in other_explanations:
                template = HYPOTHESIS_TEMPLATES_LOEUF_SEVERE.get(other['name'])
                if template:
                    hypothesis = template.format(
                        raw_value=other.get('raw_value', ''),
                        gene_name=gene_name
                    )
                    hypotheses_added.append(hypothesis)
            
            # Add all hypotheses with proper connectors
            for i, hyp in enumerate(hypotheses_added):
                if i == 0:
                    parts.append(hyp)
                else:
                    hyp_modified = hyp.replace("One hypothesis", "Another hypothesis")
                    parts.append(hyp_modified)
        
        else:
            # No explanations at all → understudied
            parts.append(UNDERSTUDIED_TEMPLATE)
        
        return " ".join(parts)
    
    # === LITERATURE MORE SEVERE: Use original logic ===
    else:
        templates = EXPLANATION_TEMPLATES_LITERATURE_SEVERE
        
        if primary:
            primary_name = primary.get('name', '')
            raw_value = primary.get('raw_value', '')
            
            if primary_name in templates:
                template = templates[primary_name]
                explanation = template.format(
                    raw_value=raw_value,
                    gene_name=gene_name
                )
                parts.append(explanation)
        
        # Add secondary explanation if present
        if secondary and secondary.get('name') in templates:
            template = templates[secondary['name']]
            explanation = template.format(
                raw_value=secondary.get('raw_value', ''),
                gene_name=gene_name
            )
            parts.append(f"Additionally, {explanation[0].lower()}{explanation[1:]}")
        
        return " ".join(parts)


def _calculate_explanations_loeuf_severe(
    gene_name: str,
    exp_missense: float = None,
    obs_lof: float = None,
    exp_lof: float = None
) -> List[Dict]:
    """
    Calculate explanations for when LOEUF is more severe than literature.
    Returns list of explanations with name, strength, raw_value.
    """
    explanations = []
    
    # Load data
    syn_data = _load_syn_data()
    tissue_spec = _load_tissue_spec_data()
    mouse_lethal = _load_mouse_lethal_genes()
    
    # 1. Synonymous Depletion (with LoF/Syn ratio comparison)
    if gene_name in syn_data:
        obs_syn = syn_data[gene_name]['obs_syn']
        exp_syn = syn_data[gene_name]['exp_syn']
        upper_bound = _gamma_upper_bound(obs_syn, exp_syn)
        
        # Activation: syn upper bound < 1.0 (there is depletion)
        if upper_bound is not None and upper_bound < 1.0:
            # Calculate simple ratios
            syn_ratio = obs_syn / exp_syn if exp_syn and exp_syn > 0 else None
            lof_ratio = obs_lof / exp_lof if obs_lof is not None and exp_lof and exp_lof > 0 else None
            
            if syn_ratio is not None and lof_ratio is not None and syn_ratio > 0:
                # R = LoF ratio / Syn ratio
                R = lof_ratio / syn_ratio
                
                # Strength based on how much LoF is explained by Syn
                if R >= 0.8:
                    strength = 1.0
                    explanation_text = "similar depletion, mutation bias likely explains constraint"
                elif R >= 0.6:
                    strength = 0.75
                    fold_diff = 1 / R if R > 0 else 999
                    explanation_text = f"LoF {fold_diff:.1f}x more depleted, mutation bias explains much but not all"
                elif R >= 0.4:
                    strength = 0.5
                    fold_diff = 1 / R if R > 0 else 999
                    explanation_text = f"LoF {fold_diff:.1f}x more depleted, mutation bias partially explains"
                elif R >= 0.2:
                    strength = 0.25
                    fold_diff = 1 / R if R > 0 else 999
                    explanation_text = f"LoF {fold_diff:.1f}x more depleted, mutation bias explains only small part"
                else:
                    strength = 0.1
                    fold_diff = 1 / R if R > 0 else 999
                    explanation_text = f"LoF {fold_diff:.1f}x more depleted, mutation bias alone unlikely to explain"
                
                explanations.append({
                    'name': 'Synonymous Depletion',
                    'strength': strength,
                    'raw_value': f"LoF obs/exp = {lof_ratio:.2f}, Syn obs/exp = {syn_ratio:.2f} — {explanation_text}"
                })
            else:
                # Fallback if we don't have LoF data: use old logic based on upper bound
                if upper_bound <= 0.7:
                    strength = 1.0
                elif upper_bound <= 0.8:
                    strength = 0.75
                elif upper_bound <= 0.9:
                    strength = 0.5
                else:
                    strength = 0.25
                explanations.append({
                    'name': 'Synonymous Depletion',
                    'strength': strength,
                    'raw_value': f"obs/exp syn upper bound = {upper_bound:.2f} (depleted)"
                })
    
    # 2. Reproductive Specificity (Testis)
    if gene_name in tissue_spec:
        testis_rank = tissue_spec[gene_name].get('testis_rank')
        if testis_rank is not None and testis_rank <= 5:
            if testis_rank == 1:
                strength = 1.0
            elif testis_rank == 2:
                strength = 0.8
            elif testis_rank == 3:
                strength = 0.6
            elif testis_rank == 4:
                strength = 0.4
            else:  # rank 5
                strength = 0.2
            explanations.append({
                'name': 'Reproductive Specificity (Testis)',
                'strength': strength,
                'raw_value': f"testis expression rank {testis_rank} out of 54 tissues"
            })
    
    # 3. Reproductive Specificity (Ovary)
    if gene_name in tissue_spec:
        ovary_rank = tissue_spec[gene_name].get('ovary_rank')
        if ovary_rank is not None and ovary_rank <= 5:
            if ovary_rank == 1:
                strength = 1.0
            elif ovary_rank == 2:
                strength = 0.8
            elif ovary_rank == 3:
                strength = 0.6
            elif ovary_rank == 4:
                strength = 0.4
            else:  # rank 5
                strength = 0.2
            explanations.append({
                'name': 'Reproductive Specificity (Ovary)',
                'strength': strength,
                'raw_value': f"ovary expression rank {ovary_rank} out of 54 tissues"
            })
    
    # 4. Mouse Embryonic Lethal
    if gene_name in mouse_lethal:
        explanations.append({
            'name': 'Mouse Embryonic Lethal',
            'strength': 1.0,
            'raw_value': "present in MGI embryonic lethal gene list"
        })
    
    # Sort by strength descending
    explanations.sort(key=lambda x: x['strength'], reverse=True)
    return explanations


def _calculate_explanations_literature_severe(
    gene_name: str,
    gof_score: int = None,
    dn_score: int = None,
    lof_score: int = None,
    onset_score: int = None,
    inheritance_score: int = None,
    penetrance_score: int = None,
    exp_missense: float = None
) -> List[Dict]:
    """
    Calculate explanations for when Literature is more severe than LOEUF.
    Returns list of explanations with name, strength, raw_value.
    """
    explanations = []
    
    # Load syn data
    syn_data = _load_syn_data()
    
    # 1. GoF mechanism
    if gof_score is not None and gof_score < 5:
        if gof_score == 1:
            strength = 1.0
            evidence = "strong evidence"
        elif gof_score == 2:
            strength = 0.75
            evidence = "reliable evidence"
        elif gof_score == 3:
            strength = 0.5
            evidence = "moderate evidence"
        elif gof_score == 4:
            strength = 0.25
            evidence = "weak evidence"
        else:
            strength = 0.0
            evidence = "no evidence"
        
        if strength > 0:
            explanations.append({
                'name': 'Gain-of-Function Mechanism',
                'strength': strength,
                'raw_value': f"GoF score = {gof_score} ({evidence} for GoF)"
            })
    
    # 2. DN mechanism
    if dn_score is not None and dn_score < 5:
        if dn_score == 1:
            strength = 1.0
            evidence = "strong evidence"
        elif dn_score == 2:
            strength = 0.75
            evidence = "reliable evidence"
        elif dn_score == 3:
            strength = 0.5
            evidence = "moderate evidence"
        elif dn_score == 4:
            strength = 0.25
            evidence = "weak evidence"
        else:
            strength = 0.0
            evidence = "no evidence"
        
        if strength > 0:
            explanations.append({
                'name': 'Dominant-Negative Mechanism',
                'strength': strength,
                'raw_value': f"DN score = {dn_score} ({evidence} for DN)"
            })
    
    # 3. Late Onset
    onset_labels = {
        1: "Prenatal", 2: "Neonatal", 3: "Infancy", 4: "Childhood",
        5: "Adolescence", 6: "Adulthood", 7: "Late-onset",
        8: "Variable", 9: "Unknown"
    }
    if onset_score is not None and onset_score in [5, 6, 7]:
        if onset_score == 7:
            strength = 1.0
        elif onset_score == 6:
            strength = 0.67
        elif onset_score == 5:
            strength = 0.33
        else:
            strength = 0.0
        
        if strength > 0:
            label = onset_labels.get(onset_score, "Unknown")
            explanations.append({
                'name': 'Late Onset',
                'strength': strength,
                'raw_value': f"onset score = {onset_score} ({label})"
            })
    
    # 4. Recessive Inheritance
    inheritance_labels = {
        1: "Autosomal dominant", 2: "Co-dominant", 3: "Incomplete dominance",
        4: "Digenic", 5: "Semidominant", 6: "Autosomal recessive",
        7: "Unknown", 8: "Conflicting"
    }
    if inheritance_score is not None and inheritance_score in [5, 6]:
        if inheritance_score == 6:
            strength = 1.0
        elif inheritance_score == 5:
            strength = 0.67
        else:
            strength = 0.0
        
        if strength > 0:
            label = inheritance_labels.get(inheritance_score, "Unknown")
            explanations.append({
                'name': 'Recessive Inheritance',
                'strength': strength,
                'raw_value': f"inheritance score = {inheritance_score} ({label})"
            })
    
    # 5. Incomplete Penetrance
    penetrance_labels = {
        1: "Complete penetrance", 2: "High penetrance",
        3: "Reduced penetrance", 4: "Low penetrance"
    }
    if penetrance_score is not None and penetrance_score in [2, 3, 4]:
        if penetrance_score == 4:
            strength = 1.0
        elif penetrance_score == 3:
            strength = 0.67
        elif penetrance_score == 2:
            strength = 0.33
        else:
            strength = 0.0
        
        if strength > 0:
            label = penetrance_labels.get(penetrance_score, "Unknown")
            explanations.append({
                'name': 'Incomplete Penetrance',
                'strength': strength,
                'raw_value': f"penetrance score = {penetrance_score} ({label})"
            })
    
    # 6. Low Expected Variants
    if exp_missense is not None and exp_missense < 30:
        strength = max(0, (30 - exp_missense) / 30)
        if strength > 0.1:  # Only include if meaningful
            explanations.append({
                'name': 'Low Statistical Power',
                'strength': strength,
                'raw_value': f"expected missense variants = {exp_missense:.1f}"
            })
    
    # 7. Synonymous Enrichment
    if gene_name in syn_data:
        obs_syn = syn_data[gene_name]['obs_syn']
        exp_syn = syn_data[gene_name]['exp_syn']
        lower_bound = _gamma_lower_bound(obs_syn, exp_syn)
        if lower_bound is not None and lower_bound >= 1.0:
            if lower_bound >= 1.5:
                strength = 1.0
            elif lower_bound >= 1.3:
                strength = 0.75
            elif lower_bound >= 1.1:
                strength = 0.5
            else:
                strength = 0.25
            explanations.append({
                'name': 'Synonymous Enrichment',
                'strength': strength,
                'raw_value': f"obs/exp lower bound = {lower_bound:.2f} (enriched)"
            })
    
    # Sort by strength descending
    explanations.sort(key=lambda x: x['strength'], reverse=True)
    return explanations


def _format_explanation_section(explanation: Dict) -> str:
    """Format an explanation for the prompt."""
    if explanation is None:
        return "No secondary explanation found."
    return f"**{explanation['name']}**\n- Raw value: {explanation['raw_value']}"


def generate_agreement_analysis(
    model: str,
    temperature: float,
    gene_name: str,
    level: int,
    loeuf_score: float,
    agreement_score: float,
    expected_lof: float,
    loeuf_data: Dict,
    gof_level_data: Dict,
    worst_gof_disease: Dict = None,
    worst_lof_disease: Dict = None,
    worst_dn_disease: Dict = None,
    worst_disease_inheritance: str = 'Unknown'
) -> Optional[Dict]:
    """
    Generates agreement analysis between LLM assessment and LOEUF scores.
    
    New pipeline (v2):
    1. If agreement >= 0.7: return canned response
    2. Calculate explanations based on direction
    3. If no explanations: return canned response
    4. If explanations found: call LLM
    
    Returns:
        Dict with keys: text, agreement_label, explanations (list of names)
    """
    agreement_direction = loeuf_data.get('agreement_direction', 'unknown')
    agreement_score_val = float(agreement_score) if agreement_score is not None else 0.0
    
    # === Case 1: Strong Agreement ===
    if agreement_score_val >= 0.7:
        print(f"INFO: Strong agreement for {gene_name} (score={agreement_score_val:.2f})")
        return {
            "text": "The literature-based assessment and the gnomAD LOEUF score are in strong agreement, reinforcing the current understanding of the gene's constraint.",
            "agreement_label": "Strong agreement",
            "explanations": []
        }
    
    # === Case 2: Disagreement - Calculate explanations ===
    agreement_label = _get_agreement_label(agreement_score_val)
    
    # Get scores from gof_level_data
    gof_score = gof_level_data.get('gof_score')
    dn_score = gof_level_data.get('dn_score')
    lof_score = gof_level_data.get('lof_score')
    
    # Get other scores from loeuf_data or other sources
    onset_score = loeuf_data.get('onset_score')
    inheritance_score = loeuf_data.get('inheritance_score')
    penetrance_score = loeuf_data.get('penetrance_score')
    exp_missense = loeuf_data.get('exp_missense') or expected_lof
    obs_lof = loeuf_data.get('obs')
    exp_lof = loeuf_data.get('exp')
    
    if agreement_direction == 'likelihood_more_pathogenic':
        # LOEUF more severe than literature
        explanations = _calculate_explanations_loeuf_severe(
            gene_name, exp_missense, obs_lof=obs_lof, exp_lof=exp_lof
        )
        prompt_name = 'agreement_analysis_loeuf_severe'
        
        # Canned response if no explanations
        canned_response = (
            f"{gene_name} is more constrained than the literature seems to predict. "
            f"It is likely that this gene is understudied compared to its importance for human health. "
            f"This gene could also be an unknown embryonic lethal gene or impact human fertility."
        )
        
    elif agreement_direction == 'similar':
        # LOEUF and literature mostly agree
        print(f"INFO: Similar agreement direction for {gene_name}, returning simple agreement message")
        return {
            "text": f"LOEUF and literature seem to mostly agree for {gene_name}.",
            "agreement_label": agreement_label,
            "explanations": []
        }
        
    else:  # prior_more_pathogenic or unknown
        # Literature more severe than LOEUF
        explanations = _calculate_explanations_literature_severe(
            gene_name=gene_name,
            gof_score=gof_score,
            dn_score=dn_score,
            lof_score=lof_score,
            onset_score=onset_score,
            inheritance_score=inheritance_score,
            penetrance_score=penetrance_score,
            exp_missense=exp_missense
        )
        prompt_name = 'agreement_analysis_literature_severe'
        
        # Canned response if no explanations
        canned_response = (
            f"{gene_name} is considered important for health by the literature, "
            f"but LOEUF values do not reflect constraint, despite the disease caused being penetrant, "
            f"dominant and early onset. The disagreement is surprising."
        )
    
    # Filter explanations with strength > 0.2
    valid_explanations = [e for e in explanations if e['strength'] > 0.2]
    
    # === Case 3: No explanations found ===
    if len(valid_explanations) == 0:
        print(f"INFO: No explanations found for {gene_name}, using canned response")
        return {
            "text": canned_response,
            "agreement_label": agreement_label,
            "explanations": []  # Empty = no explanation found
        }
    
    # === Case 4: Explanations found ===
    primary = valid_explanations[0]
    secondary = valid_explanations[1] if len(valid_explanations) > 1 else None
    
    # If primary is GoF or DN and lof_score is 5 (no LoF evidence), don't give secondary
    # The mechanism explanation is sufficient and clear
    if primary['name'] in ['Gain-of-Function Mechanism', 'Dominant-Negative Mechanism']:
        if lof_score == 5:
            secondary = None
            print(f"INFO: Primary is {primary['name']} with lof_score=5, skipping secondary")
    
    print(f"INFO: Found {len(valid_explanations)} explanations for {gene_name}")
    print(f"  Primary: {primary['name']} (strength={primary['strength']:.2f})")
    if secondary:
        print(f"  Secondary: {secondary['name']} (strength={secondary['strength']:.2f})")
    
    # Build explanation names list
    explanation_names = [primary['name']]
    if secondary:
        explanation_names.append(secondary['name'])
    
    # Determine direction for algorithmic generation
    direction = "loeuf_more_severe" if agreement_direction == 'likelihood_more_pathogenic' else "literature_more_severe"
    
    # Generate explanation algorithmically (default, no LLM)
    algorithmic_text = _generate_algorithmic_explanation(
        gene_name=gene_name,
        agreement_label=agreement_label,
        direction=direction,
        primary=primary,
        secondary=secondary,
        all_explanations=valid_explanations
    )
    
    print(f"INFO: Agreement analysis completed algorithmically for {gene_name}")
    return {
        "text": algorithmic_text,
        "agreement_label": agreement_label,
        "explanations": explanation_names
    }


def generate_dispo_analysis(gene_name: str, dispo_raw: float) -> Dict:
    """
    Generate Discovery Potential analysis text based on raw DisPo value
    and existing hypothesis data (mouse lethal, reproductive expression).

    Thresholds on raw signed disagreement:
      >12  very strong   → check hypotheses
      6-12 strong        → check hypotheses
      3-6  moderate      → brief note
      <3   low           → no actionable signal
    """
    if dispo_raw is None:
        return {"level": "unknown", "label": "N/A", "text": "No Discovery Potential data available."}

    abs_val = abs(dispo_raw)

    if abs_val >= 12:
        level, label = "very_strong", "Very strong"
    elif abs_val >= 6:
        level, label = "strong", "Strong"
    elif abs_val >= 3:
        level, label = "moderate", "Moderate"
    else:
        level, label = "low", "Low"

    if dispo_raw < 0:
        return {
            "level": "low",
            "label": "Low (literature > LOEUF)",
            "text": (
                f"Low discovery potential for {gene_name}. "
                "The disease literature suggests stronger pathogenicity than LOEUF constraint indicates. "
                "There is no evidence of an undiscovered constraining role."
            ),
        }

    if level == "low":
        return {
            "level": "low",
            "label": "Low",
            "text": (
                f"Low discovery potential for {gene_name}. "
                "LOEUF constraint and disease literature are broadly consistent; "
                "no strong signal of undiscovered pathogenic roles."
            ),
        }

    if level == "moderate":
        return {
            "level": "moderate",
            "label": "Moderate",
            "text": (
                f"Moderate discovery potential for {gene_name}. "
                "LOEUF suggests somewhat more constraint than the loss-of-function disease literature predicts, "
                "but the signal is not strong enough to draw specific hypotheses."
            ),
        }

    # strong or very_strong  → generate hypotheses
    explanations = _calculate_explanations_loeuf_severe(gene_name)

    mouse_expl = [e for e in explanations if e['name'] == 'Mouse Embryonic Lethal']
    repro_expl = [e for e in explanations
                  if e['name'].startswith('Reproductive Specificity')]

    parts = [f"{label} discovery potential for {gene_name}."]
    parts.append(
        "LOEUF indicates substantially more constraint than the loss-of-function disease literature predicts, "
        "suggesting undiscovered pathogenic roles."
    )

    hypotheses = []
    if mouse_expl:
        hypotheses.append(
            "The mouse ortholog causes embryonic lethality, suggesting this gene may be essential "
            "for early embryonic development — a phenotype rarely captured in human loss-of-function disease literature."
        )
    if repro_expl:
        for r in repro_expl:
            tissue = "testis" if "Testis" in r['name'] else "ovary"
            hypotheses.append(
                f"This gene is highly expressed in {tissue} ({r['raw_value']}), "
                f"suggesting possible selection for {'male' if tissue == 'testis' else 'female'} fertility "
                "— a phenotype underrepresented in disease genetics studies."
            )

    if hypotheses:
        for i, h in enumerate(hypotheses):
            prefix = "One hypothesis: " if i == 0 else "Another hypothesis: "
            parts.append(prefix + h)
    else:
        parts.append(
            "No specific hypothesis (mouse lethality, reproductive expression) was identified. "
            "The disease literature for this gene may not be fully developed yet "
            "(recent findings, limited number of publications), "
            "and the high constraint score warrants further research."
        )

    return {"level": level, "label": label, "text": " ".join(parts)}


def generate_agreement_analysis_with_llm(
    model: str,
    temperature: float,
    gene_name: str,
    agreement_label: str,
    direction: str,
    primary: Dict,
    secondary: Dict = None,
    canned_response: str = ""
) -> Optional[str]:
    """
    Generate agreement analysis using LLM (kept for testing/comparison).
    
    Args:
        model: LLM model to use
        temperature: LLM temperature
        gene_name: Gene symbol
        agreement_label: e.g., "Moderate disagreement"
        direction: "loeuf_more_severe" or "literature_more_severe"
        primary: Primary explanation dict
        secondary: Optional secondary explanation dict
        canned_response: Fallback response if LLM fails
    
    Returns:
        Generated explanation text or None
    """
    prompt_name = 'agreement_analysis_loeuf_severe' if direction == "loeuf_more_severe" else 'agreement_analysis_literature_severe'
    
    # Build explanation names list
    explanation_names = [primary['name']]
    if secondary:
        explanation_names.append(secondary['name'])
    
    # Load prompt
    prompt_template = load_prompt(prompt_name)
    if not prompt_template:
        print(f"ERROR: Could not load prompt {prompt_name}")
        return canned_response
    
    # Format prompt
    agreement_score_pct = 50  # Placeholder, not critical for LLM
    primary_section = f"**Primary: {primary['name']}**\n- {primary['raw_value']}"
    secondary_section = f"**Secondary: {secondary['name']}**\n- {secondary['raw_value']}" if secondary else "No secondary explanation."
    
    prompt = prompt_template.format(
        gene_name=gene_name,
        agreement_score=agreement_score_pct,
        agreement_label=agreement_label,
        primary_explanation_section=primary_section,
        secondary_explanation_section=secondary_section
    )
    
    # Call LLM
    print(f"INFO: Calling LLM for agreement analysis ({prompt_name}) for {gene_name}")
    response = call_llm(prompt, model, temperature, agent_name="agreement_analysis_agent")
    
    if response:
        cleaned = response.strip()
        if cleaned:
            print(f"INFO: LLM agreement analysis completed for {gene_name}")
            return cleaned
    
    print(f"WARNING: LLM call failed for {gene_name}")
    return canned_response
