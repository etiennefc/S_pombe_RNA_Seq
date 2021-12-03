#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

conds = ['Normal', 'H2O2', 'Heat', 'Stat', 'yAS113_yAS99'] 
input_dfs = snakemake.input.dfs
heatmap_outputs = snakemake.output.heatmap

for cond in conds:
    if cond == "yAS113_yAS99":
        ko, wt = cond.split('_')
        ko_path = [path for path in input_dfs if ko in path][0]
        wt_path = [path for path in input_dfs if wt in path][0]
        ko_df = pd.read_csv(ko_path, sep='\t')
        wt_df = pd.read_csv(wt_path, sep='\t')
        ko_df = ko_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
        wt_df = wt_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
        
        # Multiply all values by 10^6 to get a better representation on the heatmap
        ko_df, wt_df = ko_df * 1000000, wt_df * 1000000

        # Add pseudocount of 0.001 so that values of 0 become 0.001 (this way, the log10 can be applied on df)
        ko_df = ko_df + 0.001
        wt_df = wt_df + 0.001

        # Divide value of ko_df by wt_df
        div_df = ko_df.div(wt_df)

        # Create simpler column names and get log10 of all values
        cols = []
        for col in list(div_df.columns):
            simple_col = col.split('_')[1]
            cols.append(simple_col)
            div_df[col+'_log10'] = np.log10(div_df[col])
        
        div_df = div_df.filter(like='_log10', axis=1)
        
        # Create heatmap (not clustered for isotypes nor read_type)
        output_path = [path for path in heatmap_outputs if cond in path][0]
        ft.heatmap_fixed(div_df, 'plasma', '\nlog10(average KO/WT RPM ratio)', cols,
                        output_path, col_cluster=False, row_cluster=False)
        
    else:
        ip_input_paths = [path for path in input_dfs if cond in path]
        ip_path = [path for path in ip_input_paths if 'IP' in path][0]
        input_path = [path for path in ip_input_paths if 'input' in path][0]
        ip_df = pd.read_csv(ip_path, sep='\t')
        input_df = pd.read_csv(input_path, sep='\t')
        ip_df = ip_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
        input_df = input_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
        
        # Multiply all values by 10^6 to get a better representation on the heatmap
        ip_df, input_df = ip_df * 1000000, input_df * 1000000

        # Add pseudocount of 0.001 so that values of 0 become 0.001 (this way, the log10 can be applied on df)
        ip_df = ip_df + 0.001
        input_df = input_df + 0.001

        # Divide value of ip_df by input_df
        div_df = ip_df.div(input_df)

        # Create simpler column names and get log10 of all values
        cols = []
        for col in list(div_df.columns):
            simple_col = col.split('_')[1]
            cols.append(simple_col)
            div_df[col+'_log10'] = np.log10(div_df[col])

        div_df = div_df.filter(like='_log10', axis=1)

        # Create heatmap (not clustered for isotypes nor read_type)
        output_path = [path for path in heatmap_outputs if cond in path][0]
        ft.heatmap_fixed(div_df, 'plasma', f'\nlog10(average {cond} IP/input RPM ratio)', cols,
                        output_path, col_cluster=False, row_cluster=False)


