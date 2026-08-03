#!/usr/bin/env python3
"""
Merge monte_carlo_min.tsv with fetal_gene_expression_tissue_with_symbols.csv

Usage:
    python merge_monte_carlo_with_fetal.py run_016

Output:
    monte_carlo_min_with_fetal.tsv in the same directory as monte_carlo_min.tsv
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# Project root path
project_root = Path(__file__).parent.parent.parent.parent


def main():
    parser = argparse.ArgumentParser(
        description='Merge monte_carlo_min.tsv with the fetal expression data'
    )
    parser.add_argument(
        'run',
        type=str,
        help='Run name (e.g. run_016)'
    )
    parser.add_argument(
        '--input',
        type=str,
        default='monte_carlo_min.tsv',
        help='Name of the input TSV file (default: monte_carlo_min.tsv)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='monte_carlo_min_with_fetal.tsv',
        help='Name of the output TSV file (default: monte_carlo_min_with_fetal.tsv)'
    )
    parser.add_argument(
        '--fetal-file',
        type=str,
        default=None,
        help='Path to the fetal expression CSV. Falling back to the '
             'PEPPER_FETAL_EXPRESSION environment variable, then to the legacy '
             'location app/fetal_gene_expression_tissue_with_symbols.csv.'
    )

    args = parser.parse_args()

    # Paths. An absolute --input / --output is used as is (pathlib drops the
    # prefix), which is what makes it possible to work outside the original
    # working tree.
    run_dir = project_root / 'app' / 'agent_runs' / args.run
    input_path = run_dir / args.input
    output_path = run_dir / args.output
    if args.fetal_file:
        fetal_path = Path(args.fetal_file).expanduser()
    elif os.environ.get('PEPPER_FETAL_EXPRESSION'):
        fetal_path = Path(os.environ['PEPPER_FETAL_EXPRESSION']).expanduser()
    else:
        fetal_path = project_root / 'app' / 'fetal_gene_expression_tissue_with_symbols.csv'

    # Check the files
    if not input_path.exists():
        print(f"Error: {input_path} does not exist", file=sys.stderr)
        sys.exit(1)

    if not fetal_path.exists():
        print(f"Error: {fetal_path} does not exist", file=sys.stderr)
        sys.exit(1)

    # Load the data
    print(f"Loading {input_path}...", file=sys.stderr)
    mc_df = pd.read_csv(input_path, sep='\t')
    print(f"  {len(mc_df)} genes in monte_carlo_min", file=sys.stderr)

    print(f"Loading {fetal_path}...", file=sys.stderr)
    fetal_df = pd.read_csv(fetal_path)
    print(f"  {len(fetal_df)} rows in fetal expression (before aggregation)", file=sys.stderr)

    # Normalise gene_symbol to upper case
    mc_df['gene_symbol'] = mc_df['gene_symbol'].str.upper()
    fetal_df['gene_symbol'] = fetal_df['gene_symbol'].str.upper()

    # Drop the RowID column when present (gene_symbol is the key we keep)
    if 'RowID' in fetal_df.columns:
        fetal_df = fetal_df.drop(columns=['RowID'])

    # Collapse duplicates by gene_symbol, averaging the numeric columns
    # (one row per gene symbol)
    numeric_cols = fetal_df.select_dtypes(include='number').columns.tolist()
    other_cols = [c for c in fetal_df.columns if c not in numeric_cols and c != 'gene_symbol']
    if other_cols:
        # Non-numeric columns are dropped from the aggregation to avoid ambiguity
        fetal_df_agg = (
            fetal_df[['gene_symbol'] + numeric_cols]
            .groupby('gene_symbol', as_index=False)
            .mean(numeric_only=True)
        )
    else:
        fetal_df_agg = (
            fetal_df
            .groupby('gene_symbol', as_index=False)
            .mean(numeric_only=True)
        )

    print(f"  {len(fetal_df_agg)} unique gene_symbol in fetal expression (after aggregation)", file=sys.stderr)

    # Left join on gene_symbol
    print("Merge (left join on gene_symbol)...", file=sys.stderr)
    merged_df = mc_df.merge(fetal_df_agg, on='gene_symbol', how='left')

    # Stats
    n_matched = merged_df[fetal_df.columns[1]].notna().sum() if len(fetal_df.columns) > 1 else 0
    print(f"  {n_matched}/{len(mc_df)} genes with fetal data", file=sys.stderr)

    # Save
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"Results written to {output_path}", file=sys.stderr)
    print(f"  Columns: {len(merged_df.columns)}", file=sys.stderr)


if __name__ == '__main__':
    main()
