#!/usr/bin/python3
import pandas as pd
import os
import numpy as np

""" Generate a dataframe for each relevant DESeq2 comparison with only genes with
    non-null values and log2 foldchange, pval and padj columns"""
cols = ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
WT_vs_KO = pd.read_csv(snakemake.input.WT_vs_KO)
WT_vs_KO.columns = cols

# Drop NA and keep only relevant columns
WT_vs_KO = WT_vs_KO.dropna(subset=['log2FoldChange', 'pvalue', 'padj'])
WT_vs_KO = WT_vs_KO[['gene_id', 'log2FoldChange', 'pvalue', 'padj']]

# Add -log10(padj) column
WT_vs_KO['-log_padj'] = -np.log10(WT_vs_KO['padj'])
WT_vs_KO.to_csv(snakemake.output.WT_vs_KO_clean, index=False, sep="\t")
