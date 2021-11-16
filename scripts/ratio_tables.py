#!/usr/bin/python3
import pandas as pd
import numpy as np
tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
tpm_df = tpm_df.loc[:, ~tpm_df.columns.str.startswith('yAS')]  # drop KO and normal strains
tpm_df = tpm_df.replace(0, 0.0000001)  # so that ratios can be computed instead of dividing by 0

# Create column of IP/input ratio for all relevant samples
conditions = ['H2O2', 'Heat', 'Stat', 'Normal']
for cond in conditions:
    for i in [1,2,3]:
        tpm_df[f'{cond}_{i}_ratio'] = tpm_df[f'{cond}_IP_{i}'] / tpm_df[f'{cond}_input_{i}']
    tpm_df[f'{cond}_avg_ratio'] = tpm_df.filter(items=[f'{cond}_1_ratio', f'{cond}_2_ratio', f'{cond}_3_ratio']).mean(axis=1)

# Create ratios of ratios columns (ratio of stress ratio vs normal ratio)
for stress in ['H2O2', 'Heat', 'Stat']:
    tpm_df[f'{stress}_avg_ratio_vs_normal'] = tpm_df[f'{stress}_avg_ratio'] / tpm_df['Normal_avg_ratio']

final_df = tpm_df.filter(regex='ratio|gene_')  # keep ratio columns and gene_id/gene_name columns
final_df.to_csv(snakemake.output.ratio_df, sep='\t', index=False)
