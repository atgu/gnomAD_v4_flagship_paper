"""GoF p-value calculation service."""
import pandas as pd
import numpy as np
import math
from scipy.stats import gamma

# Global flag to control verbose output
VERBOSE_MODE = False


def get_first_element(array_str):
    if pd.isna(array_str) or array_str == '' or array_str == 'null' or array_str == 'None':
        return np.nan
    try:
        elements = array_str.strip('[]').split(',')
        return float(elements[0])
    except (ValueError, IndexError, AttributeError):
        return np.nan


def get_second_element(array_str):
    if pd.isna(array_str) or array_str == '' or array_str == 'null' or array_str == 'None':
        return np.nan
    try:
        elements = array_str.strip('[]').split(',')
        return float(elements[1])
    except (ValueError, IndexError, AttributeError):
        return np.nan


def generate_ci_high3(obs, exp, alpha=0.05):
    if pd.isna(obs) or pd.isna(exp) or exp == 0:
        return np.nan
    return gamma.ppf(1 - alpha, a=obs + 1, scale=1) / exp


def gamma_pvalue_fast(mis_obs, mis_exp, lof_obs, lof_exp):
    if pd.isna(mis_obs) or pd.isna(mis_exp) or pd.isna(lof_obs) or pd.isna(lof_exp) or mis_exp <= 0 or lof_exp <= 0:
        return np.nan
    a = int(mis_obs + 1)
    b = int(lof_obs + 1)
    lam1, lam2 = mis_exp, lof_exp
    
    log_p = math.log(lam1 / (lam1 + lam2))
    log_q = math.log(lam2 / (lam1 + lam2))
    log_sum = None
    for k in range(a):
        log_binom = (math.lgamma(b + k) - math.lgamma(b) - math.lgamma(k + 1))
        log_term = log_binom + k * log_p + b * log_q
        if log_sum is None:
            log_sum = log_term
        else:
            if log_sum > log_term:
                log_sum = log_sum + math.log1p(math.exp(log_term - log_sum))
            else:
                log_sum = log_term + math.log1p(math.exp(log_sum - log_term))
    return math.exp(log_sum) if log_sum is not None else np.nan


def load_gof_data_from_raw(gene, gof_df):
    """
    Simplified version of load_from_raw to get obs/exp for LoF, AM_95, and ESM_95
    from a pre-loaded dataframe.
    """
    loeuf_all = gof_df.copy()

    loeuf_all['obs_lof'] = loeuf_all['linear__new_loftee_80__adj_r'].apply(get_first_element)
    loeuf_all['exp_lof'] = loeuf_all['linear__new_loftee_80__adj_r'].apply(get_second_element)
    loeuf_all['obs_am_95'] = loeuf_all['linear__am_per_95__adj_r'].apply(get_first_element)
    loeuf_all['exp_am_95'] = loeuf_all['linear__am_per_95__adj_r'].apply(get_second_element)
    loeuf_all['obs_esm_95'] = loeuf_all['linear__esm_per_95__adj_r'].apply(get_first_element)
    loeuf_all['exp_esm_95'] = loeuf_all['linear__esm_per_95__adj_r'].apply(get_second_element)

    candidates = loeuf_all[loeuf_all['gene'].str.upper() == gene.upper()]
    if candidates.empty:
        print(f"INFO: Gene '{gene}' not found for GoF calculation.")
        return None

    mane_candidates = candidates[candidates['mane_select'] == True]
    
    def pick_with_valid_lof(df):
        valid = df[df['obs_lof'].notna() & df['exp_lof'].notna() & (df['exp_lof'] != 0)]
        if valid.empty: return None
        valid = valid.copy()
        valid['lof_high_tmp'] = valid.apply(lambda r: generate_ci_high3(r['obs_lof'], r['exp_lof']), axis=1)
        return valid.loc[valid['lof_high_tmp'].idxmin()]

    picked = pick_with_valid_lof(mane_candidates)
    if picked is None:
        picked = pick_with_valid_lof(candidates)
    if picked is None:
        picked = mane_candidates.iloc[0] if not mane_candidates.empty else candidates.iloc[0]

    row = picked.copy()
    epsilon = 1.0
    for prefix in ['am_95', 'esm_95', 'lof']:
        obs_col, exp_col = f'obs_{prefix}', f'exp_{prefix}'
        if not pd.isna(row[obs_col]): row[obs_col] += epsilon
        if not pd.isna(row[exp_col]): row[exp_col] += epsilon
            
    return {
        'obs_lof': row.get('obs_lof'), 'exp_lof': row.get('exp_lof'),
        'obs_am': row.get('obs_am_95'), 'exp_am': row.get('exp_am_95'),
        'obs_esm': row.get('obs_esm_95'), 'exp_esm': row.get('exp_esm_95')
    }


def compute_gof_pvalues(gene_name, gof_raw_data_df):
    """
    Computes GoF p-values for AM and ESM vs LoF using the gamma method.
    """
    if VERBOSE_MODE:
        print(f"DEBUG: Starting GoF calculation for gene: {gene_name}")
    
    if gof_raw_data_df is not None:
        gof_data = load_gof_data_from_raw(gene_name, gof_raw_data_df)
    else:
        print(f"ERROR: GoF data not available. Small data file was not pre-loaded.")
        return {"gof_am": None, "gof_esm": None}
        
    if VERBOSE_MODE:
        print(f"DEBUG: GoF data loaded: {gof_data}")
    
    if not gof_data:
        if VERBOSE_MODE:
            print(f"DEBUG: No GoF data found for {gene_name}")
        return {'gof_am_pvalue': None, 'gof_esm_pvalue': None}

    lof_obs, lof_exp = gof_data['obs_lof'], gof_data['exp_lof']
    if VERBOSE_MODE:
        print(f"DEBUG: LoF obs/exp: {lof_obs}/{lof_exp}")
        print(f"DEBUG: AM obs/exp: {gof_data.get('obs_am')}/{gof_data.get('exp_am')}")
        print(f"DEBUG: ESM obs/exp: {gof_data.get('obs_esm')}/{gof_data.get('exp_esm')}")
    
    am_pvalue = gamma_pvalue_fast(
        mis_obs=gof_data.get('obs_am'),
        mis_exp=gof_data.get('exp_am'),
        lof_obs=lof_obs,
        lof_exp=lof_exp
    )
    
    esm_pvalue = gamma_pvalue_fast(
        mis_obs=gof_data.get('obs_esm'),
        mis_exp=gof_data.get('exp_esm'),
        lof_obs=lof_obs,
        lof_exp=lof_exp
    )
    
    if VERBOSE_MODE:
        print(f"DEBUG: Calculated p-values - AM: {am_pvalue}, ESM: {esm_pvalue}")
    return {
        'gof_am_pvalue': am_pvalue, 
        'gof_esm_pvalue': esm_pvalue,
        'obs_lof': lof_obs,
        'exp_lof': lof_exp,
        'obs_am': gof_data.get('obs_am'),
        'exp_am': gof_data.get('exp_am'),
        'obs_esm': gof_data.get('obs_esm'),
        'exp_esm': gof_data.get('exp_esm')
    }

