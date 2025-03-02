#!/usr/bin/env python3

import pandas as pd
import argparse

def format_gwas(input_file: str, output_file: str, default_N: str = "NA"):
    df = pd.read_csv(input_file, sep='\t', dtype=str, engine='python')

    column_map = {
        'SNP': ['rsid'],
        'beta': ['beta'],
        'se': ['standard_error'],
        'EAF': ['effect_allele_frequency'],
        'A1': ['effect_allele'],
        'A2': ['other_allele'],
        'Pvalue': ['p_value']
    }

    out_cols = ['SNP', 'beta', 'se', 'EAF', 'A1', 'A2', 'Pvalue', 'N']
    out_df = pd.DataFrame(columns=out_cols)

    for final_col, possible_cols in column_map.items():
        found_col = None
        for c in possible_cols:
            if c in df.columns:
                found_col = c
                break

        if found_col is not None:
            out_df[final_col] = df[found_col]
        else:
            out_df[final_col] = 'NA'

    out_df["N"] = default_N
    out_df.to_csv(output_file, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Convert GWAS result files into a unified format with columns: SNP, beta, se, EAF, A1, A2, Pvalue, N."
    )
    parser.add_argument("--input", required=True, help="Path to input GWAS file.")
    parser.add_argument("--output", required=True, help="Path to output formatted file.")
    parser.add_argument("--N", default="NA", help="Value to fill in for N column.")
    args = parser.parse_args()

    format_gwas(args.input, args.output, args.N)


if __name__ == "__main__":
    main()
