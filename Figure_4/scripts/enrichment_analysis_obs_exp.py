#!/usr/bin/env python3
"""
Enrichment analysis with p-values from obs_exp_for_loeuf_missense.tsv.
Uses obs_missense_avg/exp_missense_avg vs obs_p_misannot_80/exp_p_misannot_80.
"""

import sys, subprocess, importlib.util, os
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from typing import Set
from sklearn.neighbors import NearestNeighbors
import argparse
import math

# ensure scipy
if importlib.util.find_spec('scipy') is None:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'scipy', '--quiet'])
from scipy.stats import fisher_exact, chi2

# Command-line arguments
parser = argparse.ArgumentParser(description='Enrichment analysis using obs_exp_for_loeuf_missense.tsv')
parser.add_argument('--data-file', type=str, default='../data/enrichment/obs_exp_for_loeuf_missense.tsv',
                   help='Input data file path')
parser.add_argument('--syn-filter', type=float, default=0.3, 
                   help='Synonymous ratio filter (0.1 = strict, 0.2 = default, 0.3 = wide)')
parser.add_argument('--include-protein-folding', action='store_true',
                   help='Add Protein folding benchmark (GO:0006457, ~261 genes)')
args = parser.parse_args()

# Synonymous ratio filter settings
SYN_FILTER_THRESHOLD = args.syn_filter
if SYN_FILTER_THRESHOLD == 0.1:
    SYN_LOW = 0.9
    SYN_HIGH = 1.1
elif SYN_FILTER_THRESHOLD == 0.2:
    SYN_LOW = 0.8
    SYN_HIGH = 1.2
else:
    # Derive bounds from threshold
    SYN_LOW = 1.0 - SYN_FILTER_THRESHOLD
    SYN_HIGH = 1.0 + SYN_FILTER_THRESHOLD

print("ENRICHMENT ANALYSIS (obs_exp_for_loeuf_missense.tsv)")
print("=" * 100)
print(f"Synonymous filter: threshold={SYN_FILTER_THRESHOLD}, range=[{SYN_LOW:.1f}, {SYN_HIGH:.1f}]")

# Load data
DATA_FILE = args.data_file
ORIGINAL_FILE = '../data/enrichment/mis_lof_with_p_values_no_syn_filter.tsv'
GOF_FILE = '../data/enrichment/g2p_pure_gof_with_p_values.tsv'
HGMD_FILE = '../data/enrichment/gofonly_HGMD2019_genes.tsv'
CANCER_FILE = '../data/enrichment/cancer_genes.tsv'
HI_FILE = '../data/enrichment/HI_table_with_category.tsv'
CHANNELS_FILE = '../data/enrichment/channels.tsv'
SCORES_FILE = '../data/enrichment/gene_scores.tsv'

print("Loading main data...")
df = pd.read_table(DATA_FILE, engine='python')
df.columns = df.columns.str.replace('"', '', regex=False).str.strip()

# Normalize gene column name
if 'gene_symbol' in df.columns:
    df['gene'] = df['gene_symbol'].astype(str).str.strip()
elif 'gene' not in df.columns:
    raise ValueError("Neither 'gene' nor 'gene_symbol' column found in input file")

# Load original file for syn_ratio
print("Loading original file for syn_ratio...")
try:
    original_df = pd.read_table(ORIGINAL_FILE, engine='python')
    original_df.columns = original_df.columns.str.replace('"', '', regex=False).str.strip()
    original_df['gene'] = original_df['gene'].astype(str).str.strip()
    
    if 'syn_ratio' in original_df.columns:
        syn_ratio_df = original_df[['gene', 'syn_ratio']].copy()
        df = df.merge(syn_ratio_df, on='gene', how='left')
        print(f"   syn_ratio merged for {df['syn_ratio'].notna().sum()} genes")
    else:
        print("   WARNING: 'syn_ratio' column not found in original file")
        df['syn_ratio'] = np.nan
except Exception as e:
    print(f"   WARNING: Error loading original file: {e}")
    print("   WARNING: No synonymous filtering applied")
    df['syn_ratio'] = np.nan

gof = pd.read_table(GOF_FILE, engine='python')
gof.columns = gof.columns.str.replace('"', '', regex=False).str.strip()

print(f"{len(df)} genes in main dataset")

SYN_FILTER_APPLIED = False

if 'syn_ratio' in df.columns and df['syn_ratio'].notna().any():
    print(f"Applying synonymous filter: ratio [{SYN_LOW:.1f}, {SYN_HIGH:.1f}]")
    before_filter = len(df)
    df = df[(df['syn_ratio'] >= SYN_LOW) & (df['syn_ratio'] <= SYN_HIGH)]
    after_filter = len(df)
    SYN_FILTER_APPLIED = True
    print(f"   Genes before filter: {before_filter}")
    print(f"   Genes after filter: {after_filter}")
    print(f"   Genes dropped: {before_filter - after_filter}")
else:
    print("WARNING: 'syn_ratio' unavailable — no synonymous filtering applied")

# Gamma p-value helper
def gamma_pvalue_fast(mis_obs, mis_exp, lof_obs, lof_exp):
    """P(mis_rate < lof_rate) via closed form (integer shapes)"""
    if pd.isna(mis_obs) or pd.isna(mis_exp) or pd.isna(lof_obs) or pd.isna(lof_exp):
        return np.nan
    if mis_exp <= 0 or lof_exp <= 0:
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

print("Computing p_missense_avg_vs_lof_80...")
df['p_missense_avg_vs_lof_80'] = df.apply(
    lambda row: gamma_pvalue_fast(
        row['obs_missense_avg'], row['exp_missense_avg'],
        row['obs_p_misannot_80'], row['exp_p_misannot_80']
    ), axis=1
)
print(f"   P-values computed for {df['p_missense_avg_vs_lof_80'].notna().sum()} genes")

# Per-gene p-value table (same genes as enrichment), sorted ascending by p-value
_syn_suffix = f"_syn{SYN_FILTER_THRESHOLD}" if SYN_FILTER_APPLIED else "_synNone"
_gene_pval_cols = ['gene', 'p_missense_avg_vs_lof_80']
_extra = [c for c in ['obs_missense_avg', 'exp_missense_avg', 'obs_p_misannot_80', 'exp_p_misannot_80', 'syn_ratio'] if c in df.columns]
_gene_pval_out = f"../data/enrichment/gene_pvalues_p_missense_avg_vs_lof_80{_syn_suffix}.tsv"
_gene_pval_df = df[_gene_pval_cols + _extra].copy()
_gene_pval_df = _gene_pval_df.sort_values('p_missense_avg_vs_lof_80', ascending=True, na_position='last')
_gene_pval_df.to_csv(_gene_pval_out, sep='\t', index=False)
print(f"   Per-gene p-value table saved: {_gene_pval_out} ({len(_gene_pval_df)} genes, ascending p-value)")

# Filtered export: genes with obs/exp p_misannot_80 ratio > 0.8
if 'obs_p_misannot_80' in df.columns and 'exp_p_misannot_80' in df.columns:
    _valid_exp = df['exp_p_misannot_80'] > 0
    _lof_ratio = np.where(_valid_exp, df['obs_p_misannot_80'].values / df['exp_p_misannot_80'].values, np.nan)
    _mask_gt08 = _lof_ratio > 0.8
    _df_gt08 = df.loc[_mask_gt08, _gene_pval_cols + _extra].copy()
    _df_gt08 = _df_gt08.sort_values('p_missense_avg_vs_lof_80', ascending=True, na_position='last')
    _gene_pval_gt08_out = f"../data/enrichment/gene_pvalues_p_missense_avg_vs_lof_80{_syn_suffix}_lof_ratio_gt_0.8.tsv"
    _df_gt08.to_csv(_gene_pval_gt08_out, sep='\t', index=False)
    print(f"   Filtered table (obs/exp p_misannot_80 > 0.8): {_gene_pval_gt08_out} ({len(_df_gt08)} genes)")
else:
    print("   WARNING: obs/exp p_misannot_80 columns missing — lof_ratio > 0.8 table not written")

methods = ['p_missense_avg_vs_lof_80']
print(f"{len(methods)} p-value metric(s) for plots: {methods}")

OMIM_MODE = 2  # 1 = by inheritance mode; 2 = pooled (raw/filtered)

DISTRIBUTION_SIDE = "low"  # or "high" or "both"

# Calculate LOF constraint ratios for matching
def calculate_lof_ratios(df):
    """Compute LOF obs/exp ratios for matching."""
    print("Computing LOF ratios for matching (obs_p_misannot_80/exp_p_misannot_80)...")
    
    if 'obs_p_misannot_80' in df.columns and 'exp_p_misannot_80' in df.columns:
        valid_exp_lof = df['exp_p_misannot_80'] > 0
        df.loc[valid_exp_lof, 'lof_ratio'] = df.loc[valid_exp_lof, 'obs_p_misannot_80'] / df.loc[valid_exp_lof, 'exp_p_misannot_80']
        df['matching_ratio'] = df['lof_ratio']  # Use LOF ratio for matching
        print(f"   LOF ratio computed for {valid_exp_lof.sum()} genes")
    else:
        print("   ERROR: obs_p_misannot_80/exp_p_misannot_80 columns not found")
        df['lof_ratio'] = np.nan
        df['matching_ratio'] = np.nan
    
    return df

df = calculate_lof_ratios(df)

# Replace matching_ratio with LOEUF-filtered score for KNN matching
print("Loading LOEUF-filtered scores for matching...")
try:
    scores_df = pd.read_csv(SCORES_FILE, sep='\t', usecols=['GENE_ID', 'LOEUF-filtered'])
    scores_df.rename(columns={'GENE_ID': 'gene', 'LOEUF-filtered': 'loeuf_filter'}, inplace=True)
    scores_df['gene'] = scores_df['gene'].astype(str).str.strip()
    df['gene'] = df['gene'].astype(str).str.strip()
    df = df.merge(scores_df, on='gene', how='left')
    df['matching_ratio'] = df['loeuf_filter']
    print(f"   LOEUF-filtered available for {df['matching_ratio'].notna().sum()} genes")
except Exception as e:
    print(f"   ERROR loading/merging LOEUF scores: {e}")

# Load oncogenes
def load_oncogenes_data():
    try:
        cancer_df = pd.read_csv(CANCER_FILE, sep='\t', on_bad_lines='skip')
        oncogenes_df = cancer_df[cancer_df['Gene Type'] == 'ONCOGENE'].copy()
        oncogene_names = set(oncogenes_df['Hugo Symbol'].tolist())
        for aliases in oncogenes_df['Gene Aliases'].dropna():
            if pd.notna(aliases) and str(aliases) != 'nan':
                alias_list = str(aliases).split(',')
                oncogene_names.update(alias.strip() for alias in alias_list if alias.strip())
        print(f"Oncogenes loaded: {len(oncogenes_df)} primary genes, {len(oncogene_names)} with aliases")
        return oncogene_names, oncogenes_df
    except Exception as e:
        print(f"ERROR loading oncogenes: {e}")
        return set(), pd.DataFrame()

def load_tsg_data():
    try:
        cancer_df = pd.read_csv(CANCER_FILE, sep='\t', on_bad_lines='skip')
        tsg_df = cancer_df[cancer_df['Gene Type'] == 'TSG'].copy()
        tsg_names = set(tsg_df['Hugo Symbol'].tolist())
        for aliases in tsg_df['Gene Aliases'].dropna():
            if pd.notna(aliases) and str(aliases) != 'nan':
                alias_list = str(aliases).split(',')
                tsg_names.update(alias.strip() for alias in alias_list if alias.strip())
        print(f"TSGs loaded: {len(tsg_df)} primary genes, {len(tsg_names)} with aliases")
        return tsg_names, tsg_df
    except Exception as e:
        print(f"ERROR loading TSG: {e}")
        return set(), pd.DataFrame()

# Load GOF_not_LOF benchmark
def load_gof_benchmark():
    g2p_all = set(gof[gof['confidence'] != 'limited']['gene'].astype(str).str.strip())
    hgmd_set = set()
    if os.path.exists(HGMD_FILE):
        hgmd = pd.read_table(HGMD_FILE, engine='python')
        hgmd.columns = hgmd.columns.str.replace('"', '', regex=False).str.strip()
        gene_col = 'gene' if 'gene' in hgmd.columns else hgmd.columns[0]
        hgmd_set = set(hgmd[gene_col].astype(str).str.strip())
    gof_set = g2p_all.union(hgmd_set)
    print(f"GOF_not_LOF: G2P ({len(g2p_all)}) + HGMD ({len(hgmd_set)}) = {len(gof_set)} genes")
    return gof_set

# Load HI genes
def load_hi_genes():
    try:
        hi_df = pd.read_csv(HI_FILE, sep='\t')
        hi_genes = set(hi_df[hi_df['Category'].isin(['severe', 'new_HI'])]['Gene'].astype(str).str.strip())
        print(f"HI genes: {len(hi_genes)} (severe + new_HI)")
        return hi_genes
    except Exception as e:
        print(f"ERROR loading HI: {e}")
        return set()

# Load Channels genes
def load_channels_genes():
    try:
        channels_df = pd.read_csv(CHANNELS_FILE, sep='\t')
        channels_genes = set(channels_df['Approved symbol'].astype(str).str.strip())
        print(f"Channel genes: {len(channels_genes)} ion channel genes")
        return channels_genes, channels_df
    except Exception as e:
        print(f"ERROR loading Channels: {e}")
        return set(), pd.DataFrame()

# Other benchmark loaders
def load_or_genes():
    try:
        or_file = '../data/enrichment/olfactory_receptors.tsv'
        with open(or_file, 'r') as f:
            or_genes = set(line.strip() for line in f if line.strip())
        print(f"OR genes: {len(or_genes)} olfactory receptor genes")
        return or_genes
    except Exception as e:
        print(f"ERROR loading OR: {e}")
        return set()

def load_kinases_genes():
    try:
        kinases_file = '../data/enrichment/kinases.tsv'
        with open(kinases_file, 'r') as f:
            kinases_genes = set(line.strip() for line in f if line.strip())
        print(f"Kinase genes: {len(kinases_genes)}")
        return kinases_genes
    except Exception as e:
        print(f"ERROR loading Kinases: {e}")
        return set()

def load_all_omim_raw_genes():
    try:
        omim_file = '../data/enrichment/omim.tsv'
        omim_df = pd.read_csv(omim_file, sep='\t')
        all_genes = set()
        for genes_str in omim_df['genes']:
            if pd.notna(genes_str):
                genes = genes_str.split('|')
                for gene in genes:
                    if gene.strip() and gene.strip() != 'NA':
                        all_genes.add(gene.strip())
        print(f"All OMIM raw genes: {len(all_genes)} (unfiltered)")
        return all_genes
    except Exception as e:
        print(f"ERROR loading All OMIM raw: {e}")
        return set()

def load_all_omim_genes():
    try:
        omim_file = '../data/enrichment/omim.tsv'
        omim_df = pd.read_csv(omim_file, sep='\t')
        all_omim_genes = set()
        for genes_str in omim_df['genes']:
            if pd.notna(genes_str):
                genes = genes_str.split('|')
                for gene in genes:
                    if gene.strip() and gene.strip() != 'NA':
                        all_omim_genes.add(gene.strip())
        print(f"All OMIM genes: {len(all_omim_genes)} (before filtering)")
        oncogenes, _ = load_oncogenes_data()
        channels_genes, _ = load_channels_genes()
        gof_genes = load_gof_benchmark()
        helicases_genes = load_helicases_genes()
        kinases_genes = load_kinases_genes()
        dimer_genes = load_dimer_genes()
        genes_to_exclude = oncogenes.union(channels_genes).union(gof_genes).union(helicases_genes).union(kinases_genes).union(dimer_genes)
        filtered_genes = all_omim_genes - genes_to_exclude
        print(f"All OMIM genes (filtered): {len(filtered_genes)}")
        return filtered_genes
    except Exception as e:
        print(f"ERROR loading All OMIM: {e}")
        return set()

def load_dimer_genes():
    try:
        dimer_file = '../data/enrichment/dimeres.tsv'
        dimer_df = pd.read_csv(dimer_file, sep='\t')
        dimer_genes = set(dimer_df['SYMBOL'].dropna().unique())
        print(f"Dimer genes: {len(dimer_genes)}")
        return dimer_genes
    except Exception as e:
        print(f"ERROR loading Dimer: {e}")
        return set()

def load_helicases_genes():
    try:
        helicases_file = '../data/enrichment/helicases.tsv'
        helicases_df = pd.read_csv(helicases_file, sep='\t')
        helicases_genes = set(helicases_df['SYMBOL'].dropna().unique())
        print(f"Helicase genes: {len(helicases_genes)}")
        return helicases_genes
    except Exception as e:
        print(f"ERROR loading Helicases: {e}")
        return set()

def load_protein_folding_genes():
    try:
        pf_file = '../data/enrichment/protein_folding_GO0006457.txt'
        with open(pf_file, 'r') as f:
            genes = set(line.strip() for line in f if line.strip())
        print(f"Protein folding genes: {len(genes)} (GO:0006457)")
        return genes
    except Exception as e:
        print(f"ERROR loading Protein folding: {e}")
        return set()

def match_benchmark_with_controls(df, bench_set, ratio_col, n_matches=3):
    """Match benchmark genes with controls by similar constraint ratios"""
    
    print(f"Matching by {ratio_col} (n_matches={n_matches})...")
    
    # Filter valid ratios
    valid_ratio_mask = df[ratio_col].notna() & np.isfinite(df[ratio_col])
    df_valid = df[valid_ratio_mask].copy()
    
    if len(df_valid) == 0:
        print("   ERROR: No valid ratio for matching")
        return df
    
    # Identify benchmark genes in valid set
    bench_mask = df_valid['gene'].astype(str).str.strip().isin(bench_set)
    bench_genes = df_valid[bench_mask]
    other_genes = df_valid[~bench_mask]
    
    if len(bench_genes) == 0:
        print("   ERROR: No benchmark genes with valid ratio")
        return df_valid
    
    if len(other_genes) == 0:
        print("   ERROR: No control genes with valid ratio")
        return df_valid
    
    print(f"   Benchmark: {len(bench_genes)} genes, Controls: {len(other_genes)} genes")
    
    # Fit KNN model
    X_other = other_genes[[ratio_col]].values
    X_bench = bench_genes[[ratio_col]].values
    
    # Use min of n_matches and available controls
    actual_n_matches = min(n_matches, len(other_genes))
    
    nbrs = NearestNeighbors(n_neighbors=actual_n_matches, metric='euclidean')
    nbrs.fit(X_other)
    
    # Find matches for each benchmark gene
    distances, indices = nbrs.kneighbors(X_bench)
    
    # Collect matched control genes
    matched_control_indices = set()
    for gene_matches in indices:
        matched_control_indices.update(gene_matches)
    
    matched_controls = other_genes.iloc[list(matched_control_indices)]
    
    print(f"   {len(matched_controls)} matched control genes")
    
    # Combine benchmark + matched controls
    result_df = pd.concat([bench_genes, matched_controls], ignore_index=True)
    
    print(f"   Final dataset: {len(result_df)} genes ({len(bench_genes)} benchmark + {len(matched_controls)} controls)")
    
    return result_df

def calculate_enrichment_complete_matched(df, bench_set, bench_name, thresh, methods_list, use_matching=True):
    """Calculate enrichment for ALL methods with optional matching"""
    
    if DISTRIBUTION_SIDE == "low":
        print(f"\nFull analysis for {bench_name} (p < {thresh}) - Matching: {use_matching}")
    elif DISTRIBUTION_SIDE == "high":
        print(f"\nFull analysis for {bench_name} (p > {1-thresh:.3f}) - Matching: {use_matching}")
    else:
        print(f"\nFull analysis for {bench_name} (p < {thresh} and p > {1-thresh:.3f}) - Matching: {use_matching}")
    
    # Apply matching if requested
    if use_matching and 'matching_ratio' in df.columns:
        df_analysis = match_benchmark_with_controls(df, bench_set, 'matching_ratio', n_matches=3)
    else:
        df_analysis = df.copy()
    
    # Add benchmark flag
    in_bench_mask = df_analysis['gene'].astype(str).str.strip().isin(bench_set)
    
    results = {}
    
    for col in methods_list:
        if col not in df_analysis.columns:
            print(f"   WARNING: {col}: missing column")
            continue
            
        # Filter valid values
        valid_mask = df_analysis[col].notna()
        df_valid = df_analysis[valid_mask]
        in_bench_valid = in_bench_mask[valid_mask]
        
        # Distribution-side filter
        if DISTRIBUTION_SIDE == "low":
            sig_mask = df_valid[col] < thresh
        elif DISTRIBUTION_SIDE == "high":
            sig_mask = df_valid[col] > (1 - thresh)
        else:  # both — default to low for this iteration
            sig_mask = df_valid[col] < thresh
            
        n_sig = sig_mask.sum()
        
        # MODIFICATION: Always include the method, even if n_sig = 0
        if n_sig == 0:
            # Set enrichment to 0 instead of skipping
            results[col] = {
                'enrichment': 0.0,
                'fisher_p': 1.0,
                'bench_sig': 0,
                'total_sig': 0,
                'bench_total': in_bench_valid.sum()
            }
            print(f"   {col}: 0/{n_sig} genes (enrichment = 0.0)")
            continue
            
        n_bench_sig = (sig_mask & in_bench_valid).sum()
        bench_prop_sig = n_bench_sig / n_sig
        bench_prop_all_valid = in_bench_valid.mean()
        
        enrichment = bench_prop_sig / bench_prop_all_valid if bench_prop_all_valid > 0 else 0.0
        
        # Fisher exact test
        a = n_bench_sig
        b = in_bench_valid.sum() - a
        c = n_sig - a
        d = len(df_valid) - a - b - c
        
        if a + b > 0 and c + d > 0:
            fisher_p = fisher_exact([[a, b], [c, d]], alternative='two-sided')[1]
        else:
            fisher_p = 1.0
        
        results[col] = {
            'enrichment': enrichment,
            'fisher_p': fisher_p,
            'bench_sig': n_bench_sig,
            'total_sig': n_sig,
            'bench_total': in_bench_valid.sum()
        }
        
        print(f"   {col}: {n_bench_sig}/{n_sig} genes (enrichment = {enrichment:.2f})")
    
    return results

def calculate_fisher_order(results_by_bench, metric_key):
    """Order benchmarks by the single p-value metric."""
    method_key = 'p_missense_avg_vs_lof_80'
    
    fisher_values = {}
    for bench_name, bench_results in results_by_bench.items():
        if method_key in bench_results:
            if metric_key == "pvalue":
                fisher_values[bench_name] = 1 - bench_results[method_key]["fisher_p"]
            else:
                fisher_values[bench_name] = bench_results[method_key]["enrichment"]
        else:
            fisher_values[bench_name] = 0.0
    
    # Sort descending by metric
    sorted_benchmarks = sorted(fisher_values.items(), key=lambda x: x[1], reverse=True)
    return [bench for bench, _ in sorted_benchmarks]

def save_results_to_tsv(results_by_bench, thresh, matched, background_counts, final_benchmark_counts, output_file):
    """Write enrichment results to TSV for ggplot."""
    rows = []
    
    for bench_name, bench_results in results_by_bench.items():
        for method, result in bench_results.items():
            row = {
                'benchmark': bench_name,
                'method': method,
                'enrichment': result['enrichment'],
                'fisher_p': result['fisher_p'],
                'neg_log10_fisher_p': -math.log10(result['fisher_p']) if result['fisher_p'] > 0 else 0,
                'bench_sig': result['bench_sig'],
                'total_sig': result['total_sig'],
                'bench_total': result['bench_total'],
                'matched': 'matched' if matched else 'unmatched',
                'threshold': thresh,
                'background_size': background_counts.get(bench_name, 0) if background_counts else 0,
                'benchmark_gene_count': final_benchmark_counts.get(bench_name, 0) if final_benchmark_counts else 0,
                'syn_filter': SYN_FILTER_THRESHOLD if SYN_FILTER_APPLIED else None,
                'clean_background': CLEAN_BACKGROUND
            }
            rows.append(row)
    
    results_df = pd.DataFrame(rows)
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved results to {output_file}")

def create_generic_barplot(results_by_bench, metric_key, thresh, matched, show_only_enrichment=False, background_counts=None, final_benchmark_counts=None):
    """Generic bar plot: enrichment or -log10 p-value for n benchmarks / one method."""
    plot_type = "enrichment" if metric_key == "enrichment" else "pvalue"

    # Union of methods across benchmarks
    all_methods = set()
    for res in results_by_bench.values():
        all_methods.update(res.keys())
    
    # Fixed method order (single method)
    custom_order = ['p_missense_avg_vs_lof_80']
    methods_sorted = [m for m in custom_order if m in all_methods]
    
    # Reorder benchmarks
    fisher_order = calculate_fisher_order(results_by_bench, metric_key)
    results_by_bench = {bench: results_by_bench[bench] for bench in fisher_order}

    x = np.arange(len(methods_sorted))
    n_bench = len(results_by_bench)
    width = 0.8 / n_bench

    fig, ax = plt.subplots(figsize=(max(10, len(methods_sorted)*0.8), 6))

    # tab20-like distinct colors
    colors = plt.cm.get_cmap('tab20', n_bench)

    for i, (bench_name, bench_results) in enumerate(results_by_bench.items()):
        values = []
        for m in methods_sorted:
            if m in bench_results:
                if metric_key == "pvalue":
                    v = bench_results[m]["fisher_p"]
                else:
                    v = bench_results[m][metric_key]
            else:
                v = 0.0
            # pvalue mode → -log10
            if metric_key == "pvalue":
                v = -math.log10(v) if v > 0 else 0
            values.append(v)
        # Final gene count if available
        if final_benchmark_counts and bench_name in final_benchmark_counts:
            gene_count = final_benchmark_counts[bench_name]
        else:
            gene_count = benchmark_gene_counts.get(bench_name, 0)
        legend_label = f"{bench_name} ({gene_count})"
        bars = ax.bar(x + (i - n_bench/2)*width + width/2, values, width, label=legend_label, color=colors(i))
        # Bar annotations
        for bar, m in zip(bars, methods_sorted):
            height = bar.get_height()
            if metric_key == "enrichment":
                info = results_by_bench[bench_name][m]
                if show_only_enrichment:
                    label = f"{info['fisher_p']:.1e}"
                else:
                    label = f"{info['bench_sig']}/{info['total_sig']}\n(p={info['fisher_p']:.1e})"
            else:  # pvalue plot: show enrichment + counts
                info = results_by_bench[bench_name][m]
                if show_only_enrichment:
                    label = f"{info['enrichment']:.2f}"
                else:
                    label = f"enr={info['enrichment']:.2f}\n{info['bench_sig']}/{info['total_sig']}"
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.05, label, ha='center', va='bottom', fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels([])
    ylabel = "Enrichment ratio" if metric_key == "enrichment" else "-log10(p-value)"
    ax.set_ylabel(ylabel)
    match_txt = "_lof_matched" if matched else ""
    
    # Title from distribution mode and background size
    total_background = 0
    if background_counts:
        total_background = sum(background_counts.values())
        avg_background = total_background // len(background_counts) if background_counts else 0
        background_info = f" (bg: ~{avg_background:,})"
    else:
        background_info = ""
    
    if DISTRIBUTION_SIDE == "low":
        ax.set_title(f"{plot_type.capitalize()} comparison (p < {thresh}){match_txt}{background_info}")
        p_threshold = str(thresh).replace('.', '')
    elif DISTRIBUTION_SIDE == "high":
        ax.set_title(f"{plot_type.capitalize()} comparison (p > {1-thresh:.3f}){match_txt}{background_info}")
        p_threshold = f"{str(1-thresh).replace('.', '')}_high"
    else:  # both
        ax.set_title(f"{plot_type.capitalize()} comparison (p < {thresh} and p > {1-thresh:.3f}){match_txt}{background_info}")
        p_threshold = str(thresh).replace('.', '')

    if metric_key == "pvalue":
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=1)
    else:
        ax.axhline(y=1, color='red', linestyle='--', linewidth=1)
    ax.legend()
    plt.tight_layout()
    suffix = "_enrichment_only" if show_only_enrichment else ""
    background_suffix = "_clean_bg" if CLEAN_BACKGROUND else "_full_bg"
    pvalue_suffix = "_gamma"
    syn_filter_suffix = f"_syn{SYN_FILTER_THRESHOLD}" if SYN_FILTER_APPLIED else ""
    out_name = f"../figures/{plot_type}_comparison_generic{match_txt}_p{p_threshold}{suffix}{background_suffix}{pvalue_suffix}{syn_filter_suffix}.png"
    plt.savefig(out_name, dpi=300)
    plt.close()
    print(f"Saved {out_name}")

def create_clean_background_for_benchmark(df, current_benchmark_genes, all_benchmark_genes):
    """Build a clean background excluding all benchmark genes."""
    
    dataset_genes = set(df['gene'].astype(str).str.strip())
    
    all_benchmark_genes_union = set()
    for bench_name, genes in all_benchmark_genes.items():
        all_benchmark_genes_union.update(genes)
    
    all_benchmark_genes_in_dataset = all_benchmark_genes_union.intersection(dataset_genes)
    
    clean_background_genes = dataset_genes - all_benchmark_genes_in_dataset
    print(f"   Shared background: {len(clean_background_genes)} genes (all benchmarks excluded)")
    
    # Keep clean background + current benchmark genes only
    clean_background_mask = df['gene'].astype(str).str.strip().isin(clean_background_genes)
    benchmark_mask = df['gene'].astype(str).str.strip().isin(current_benchmark_genes)
    combined_mask = clean_background_mask | benchmark_mask
    
    df_clean = df[combined_mask].copy()
    
    print(f"   Dataset for benchmark: {len(df_clean)} genes (background + benchmark)")
    return df_clean

# ---------------------------------------------------------------------------
# General configuration (edit benchmark loaders here)
# ---------------------------------------------------------------------------
if OMIM_MODE == 1:
    BENCHMARK_LOADERS = {
        "Oncogenes": lambda: load_oncogenes_data()[0],
        "TSG": lambda: load_tsg_data()[0],
        "GOF genes": lambda: load_gof_benchmark(),
        "Channels": lambda: load_channels_genes()[0],
        "HI": lambda: load_hi_genes(),
        "Kinases": lambda: load_kinases_genes(),
        "OMIM filtered": lambda: load_all_omim_genes(),
        "OMIM": lambda: load_all_omim_raw_genes(),
        "Dimer": lambda: load_dimer_genes(),
        "Helicases": lambda: load_helicases_genes(),
    }
else:
    BENCHMARK_LOADERS = {
        "Oncogenes": lambda: load_oncogenes_data()[0],
        "TSG": lambda: load_tsg_data()[0],
        "GOF genes": lambda: load_gof_benchmark(),
        "Channels": lambda: load_channels_genes()[0],
        "HI": lambda: load_hi_genes(),
        "Kinases": lambda: load_kinases_genes(),
        "OMIM filtered": lambda: load_all_omim_genes(),
        "OMIM": lambda: load_all_omim_raw_genes(),
        "Dimer": lambda: load_dimer_genes(),
        "Helicases": lambda: load_helicases_genes(),
    }
if args.include_protein_folding:
    BENCHMARK_LOADERS["Protein folding"] = lambda: load_protein_folding_genes()

CLEAN_BACKGROUND = True

# Per threshold: which matplotlib plots to generate
THRESHOLD_CONFIG = {
    0.1: {"enrichment": False, "pvalue": True},
}

benchmark_sets = {}
benchmark_gene_counts = {}
for bench_name, loader_fn in BENCHMARK_LOADERS.items():
    gene_set = set(loader_fn())
    benchmark_sets[bench_name] = gene_set
    benchmark_gene_counts[bench_name] = len(gene_set)
    print(f"{bench_name}: {len(gene_set)} genes loaded")

# Results cache for optional line plots
all_results_by_benchmark = {bench_name: {} for bench_name in benchmark_sets.keys()}

for thresh, plot_dict in THRESHOLD_CONFIG.items():
    print(f"\nProcessing threshold p < {thresh}")
    
    if DISTRIBUTION_SIDE == "both":
        distribution_modes = ["low", "high"]
    else:
        distribution_modes = [DISTRIBUTION_SIDE]
    
    for distribution_mode in distribution_modes:
        original_mode = DISTRIBUTION_SIDE
        DISTRIBUTION_SIDE = distribution_mode
        
        for matched in [True, False]:
            print("MATCHED" if matched else "UNMATCHED")
            results_by_bench = {}
            for bench_name, gene_set in benchmark_sets.items():
                print(f"Analyzing {bench_name}...")
                
                if CLEAN_BACKGROUND:
                    all_benchmarks = {name: genes for name, genes in benchmark_sets.items()}
                    df_clean_for_benchmark = create_clean_background_for_benchmark(df, gene_set, all_benchmarks)
                    print(f"   Clean shared background enabled")
                else:
                    df_clean_for_benchmark = df.copy()
                    print(f"   Full dataset used as background")
                
                res = calculate_enrichment_complete_matched(df_clean_for_benchmark, gene_set, bench_name, thresh, methods, use_matching=matched)
                results_by_bench[bench_name] = res
                
                # Cache for line plots (matched=True, low mode only)
                if matched and distribution_mode == "low":
                    all_results_by_benchmark[bench_name][thresh] = res
            
            background_counts = {}
            final_benchmark_counts = {}
            
            if CLEAN_BACKGROUND:
                all_benchmark_genes_union = set()
                for name, genes in benchmark_sets.items():
                    all_benchmark_genes_union.update(genes)
                
                dataset_genes = set(df['gene'].astype(str).str.strip())
                all_benchmark_genes_in_dataset = all_benchmark_genes_union.intersection(dataset_genes)
                clean_background_genes = dataset_genes - all_benchmark_genes_in_dataset
                background_size = len(clean_background_genes)
                
                for bench_name in benchmark_sets.keys():
                    background_counts[bench_name] = background_size
                
                for bench_name, gene_set in benchmark_sets.items():
                    benchmark_genes_in_dataset = gene_set.intersection(dataset_genes)
                    final_benchmark_counts[bench_name] = len(benchmark_genes_in_dataset)
                    
                print(f"   Shared background: {background_size:,} genes (all benchmarks excluded)")
                print(f"   Final genes per benchmark: {final_benchmark_counts}")
            else:
                total_genes = len(df)
                for bench_name in benchmark_sets.keys():
                    background_counts[bench_name] = total_genes
                    final_benchmark_counts[bench_name] = benchmark_gene_counts[bench_name]
            
            match_txt = "_matched" if matched else "_unmatched"
            p_threshold = str(thresh).replace('.', '')
            syn_suffix = f"_syn{SYN_FILTER_THRESHOLD}" if SYN_FILTER_APPLIED else ""
            tsv_output = f"../data/enrichment/enrichment_results_p{p_threshold}{match_txt}{syn_suffix}.tsv"
            save_results_to_tsv(results_by_bench, thresh, matched, background_counts, final_benchmark_counts, tsv_output)
            
            # Optional matplotlib outputs
            if plot_dict.get("enrichment", False):
                create_generic_barplot(results_by_bench, "enrichment", thresh, matched, show_only_enrichment=True, background_counts=background_counts, final_benchmark_counts=final_benchmark_counts)
            if plot_dict.get("pvalue", False):
                create_generic_barplot(results_by_bench, "pvalue", thresh, matched, show_only_enrichment=True, background_counts=background_counts, final_benchmark_counts=final_benchmark_counts)
        
        DISTRIBUTION_SIDE = original_mode

print("\nAnalysis complete (obs_exp_for_loeuf_missense.tsv configuration).")
print("All plots use metric p_missense_avg_vs_lof_80.")
print("Matched and unmatched LOF-constraint versions were generated.")
print(f"Background mode: {'clean' if CLEAN_BACKGROUND else 'full'}")
