#!/usr/bin/python3
import pandas as pd
import subprocess as sp
""" Generate a dataframe of gene biotype per gene_id from a gtf and merge that
    column to an abundance df. Add an average TPM column per sample type."""

cols = ['chr', 'source', 'feature', 'start', 'end', 'score1', 'strand', 'score2', 'features']
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', names=cols)
gtf = gtf[gtf['feature'] == 'gene']
output = snakemake.output.gene_biotype_df
# Get gene biotype for all genes in gtf from the features column
gtf_df = gtf[['features']]
gtf_df = gtf_df['features'].str.split(';', expand=True)
gtf_df = gtf_df.iloc[:, 0:3]
gtf_df.columns = ['gene_id', 'gene_name', 'gene_biotype']
gtf_df = gtf_df.drop(columns=['gene_name'])
gtf_df['gene_id'] = gtf_df['gene_id'].str.replace("gene_id ", "")
gtf_df['gene_id'] = gtf_df['gene_id'].str.replace('"', '')
gtf_df['gene_biotype'] = gtf_df['gene_biotype'].str.replace("gene_biotype ", "")
gtf_df['gene_biotype'] = gtf_df['gene_biotype'].str.replace('"', '')

# Tweak for problematic gene in gtf
gtf_df = gtf_df.replace('gene_id SPNCRNA.1568', 'mRNA')

gtf_df.to_csv(output, index=False, sep="\t")
sp.call(f"sed -iE 's/gene_id.SPNCRNA.1568/mRNA/g' {output}", shell=True)
