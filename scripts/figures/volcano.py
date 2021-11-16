#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess as sp
import os

wt_ko_df = pd.read_csv(snakemake.input.deseq_wt_ko, sep='\t')
biotype_df = pd.read_csv(snakemake.input.biotype_df, sep='\t')
stress_normal_dir = snakemake.input.deseq_ip_stress_normal_dir
input_stress_normal_dir = snakemake.input.deseq_input_stress_normal_dir
output_dir = snakemake.output.volcano_output_dir
biotype_colors_dict = snakemake.params.biotype_colors

# Create output dir of volcano plots
sp.call(f'mkdir -p {output_dir}', shell=True)

# Merge biotype column to the wt_ko DESeq df
wt_ko_df = wt_ko_df.merge(biotype_df, how='left', on='gene_id')
wt_ko_df['gene_biotype'] = wt_ko_df['gene_biotype'].str.lstrip()  # remove white space before string

# Generate a volcano plot for KO compared to WT (hue: gene_biotype)
ft.volcano(wt_ko_df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            'Log2 of the fold change between\nKO and WT samples',
            '-log10(adjusted p-value)',
            'Abundance of all RNAs\nin sla1 KO compared to WT',
            biotype_colors_dict, f'{output_dir}/knockout-wild_type_biotype.svg',
            alpha=0.7)


# Generate a volcano plot for the different stress ip vs normal ip comparisons and
# also stress ip vs other stress ip (hue: gene_biotype)
for file in os.listdir(stress_normal_dir):
    if '.tsv' in file:
        cond1, cond2 = file.split('-')
        cond2 = cond2.split('_v2.')[0]
        df = pd.read_csv(f'{stress_normal_dir}/{file}', sep='\t')
        df = df.merge(biotype_df, how='left', on='gene_id')
        df['gene_biotype'] = df['gene_biotype'].str.lstrip()  # remove white space before string
        ft.volcano(df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            f'Log2 of the fold change between\n{cond1} and {cond2} samples',
            '-log10(adjusted p-value)',
            f'Abundance of all RNAs\nin {cond1} compared to {cond2}',
            biotype_colors_dict, f'{output_dir}/{cond1}-{cond2}_biotype.svg',
            alpha=0.7)


# Generate a volcano plot for the different stress inputs vs normal input comparisons (hue: gene_biotype)
for file in os.listdir(input_stress_normal_dir):
    if '.tsv' in file:
        cond1, cond2 = file.split('-')
        cond2 = cond2.split('_v2.')[0]
        df = pd.read_csv(f'{input_stress_normal_dir}/{file}', sep='\t')
        df = df.merge(biotype_df, how='left', on='gene_id')
        df['gene_biotype'] = df['gene_biotype'].str.lstrip()  # remove white space before string
        ft.volcano(df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            f'Log2 of the fold change between\n{cond1} and {cond2} samples',
            '-log10(adjusted p-value)',
            f'Abundance of all RNAs\nin {cond1} compared to {cond2}',
            biotype_colors_dict, f'{output_dir}/{cond1}-{cond2}_biotype.svg',
            alpha=0.7)
