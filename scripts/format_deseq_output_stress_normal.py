#!/usr/bin/python3
import pandas as pd
import os
import numpy as np
import subprocess as sp

""" Generate a dataframe for each relevant DESeq2 comparison with only genes with
    non-null values and log2 foldchange, pval and padj columns"""
input_dir = snakemake.input.stress_normal_deseq_dir
output_dir = snakemake.output.output_dir
sp.call(f'mkdir -p {output_dir}', shell=True)
cols = ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']

for file in os.listdir(input_dir):
    name = file.split('.')[0]
    if '.csv' in file:
        df = pd.read_csv(f'{input_dir}/{file}', names=cols, header=0)
        # Drop NA and keep only relevant columns
        df = df.dropna(subset=['log2FoldChange', 'pvalue', 'padj'])
        # Add -log10(padj) column
        df['-log_padj'] = -np.log10(df['padj'])
        df = df.drop(columns=['baseMean', 'lfcSE', 'stat'])
        df.to_csv(f'{output_dir}/{name}_v2.tsv', index=False, sep='\t')
