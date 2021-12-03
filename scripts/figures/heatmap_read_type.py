#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

heatmap_outputs = snakemake.output.heatmap

for i, avg_df in enumerate(snakemake.input.dfs):
    df = pd.read_csv(avg_df, sep='\t')
    df = df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
    cond = avg_df.split('/')[2]
    cond = cond.split('_avg')[0]
    
    # Multiply all values by 10^6 to get a better representation on the heatmap
    df = df * 1000000
    # Add pseudocount of 0.001 so that values of 0 become 0.001 (this way, the log10 can be applied on df)
    df = df + 0.001
    # Create simpler column names and get log10 of all values
    cols = []
    for col in list(df.columns):
        simple_col = col.split('_')[1]
        cols.append(simple_col)
        df[col+'_log10'] = np.log10(df[col])
    # Keep only log10 columns
    df = df.filter(like='_log10', axis=1)
    
    # Create heatmap (not clustered)
    output_path = [path for path in heatmap_outputs if cond in path][0]
    ft.heatmap_fixed(df, 'plasma', '\nlog10(average RPM)', cols,
                output_path, col_cluster=False, row_cluster=False, method='weighted')

