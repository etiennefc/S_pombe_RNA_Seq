#!/usr/bin/python3
import pandas as pd

# Load all dfs per condition
normal_input, normal_ip, h2o2_input, h2o2_ip, heat_input, heat_ip, stat_input, stat_ip, yas99, yas113 = [], [], [], [], [], [], [], [], [], []
for df_path in snakemake.input.normalized_dfs:
    if 'R1' in df_path:  # select only reads R1 (not R2)
        if 'Normal_input' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            normal_input.append(df)
        elif 'Normal_IP' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            normal_ip.append(df)
        elif 'H2O2_input' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            h2o2_input.append(df)
        elif 'H2O2_IP' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            h2o2_ip.append(df)
        elif 'Heat_input' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            heat_input.append(df)
        elif 'Heat_IP' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            heat_ip.append(df)
        elif 'Stat_input' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            stat_input.append(df)
        elif 'Stat_IP' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            stat_ip.append(df)
        elif 'yAS99' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            yas99.append(df)            
        else:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            yas113.append(df)

# Concat dfs per condition
normal_input_concat, normal_ip_concat, h2o2_input_concat, h2o2_ip_concat, heat_input_concat = pd.concat(normal_input), pd.concat(normal_ip), pd.concat(h2o2_input), pd.concat(h2o2_ip), pd.concat(heat_input)
heat_ip_concat, stat_input_concat, stat_ip_concat, yas99_concat, yas113_concat = pd.concat(heat_ip), pd.concat(stat_input), pd.concat(stat_ip), pd.concat(yas99), pd.concat(yas113)

# Get an average value per condition
normal_input_avg = normal_input_concat.groupby('isotype_anticodon').mean()
normal_ip_avg = normal_ip_concat.groupby('isotype_anticodon').mean()
h2o2_input_avg = h2o2_input_concat.groupby('isotype_anticodon').mean()
h2o2_ip_avg = h2o2_ip_concat.groupby('isotype_anticodon').mean()
heat_input_avg = heat_input_concat.groupby('isotype_anticodon').mean()
heat_ip_avg = heat_ip_concat.groupby('isotype_anticodon').mean()
stat_input_avg = stat_input_concat.groupby('isotype_anticodon').mean()
stat_ip_avg = stat_ip_concat.groupby('isotype_anticodon').mean()
yas99_avg = yas99_concat.groupby('isotype_anticodon').mean()
yas113_avg = yas113_concat.groupby('isotype_anticodon').mean()

# Select normalized columns only and isotype_anticodon
normal_input_avg = normal_input_avg.filter(like='_norm', axis=1).reset_index()
normal_ip_avg = normal_ip_avg.filter(like='_norm', axis=1).reset_index()
h2o2_input_avg = h2o2_input_avg.filter(like='_norm', axis=1).reset_index()
h2o2_ip_avg = h2o2_ip_avg.filter(like='_norm', axis=1).reset_index()
heat_input_avg = heat_input_avg.filter(like='_norm', axis=1).reset_index()
heat_ip_avg = heat_ip_avg.filter(like='_norm', axis=1).reset_index()
stat_input_avg = stat_input_avg.filter(like='_norm', axis=1).reset_index()
stat_ip_avg = stat_ip_avg.filter(like='_norm', axis=1).reset_index()
yas99_avg = yas99_avg.filter(like='_norm', axis=1).reset_index()
yas113_avg = yas113_avg.filter(like='_norm', axis=1).reset_index()


normal_input_avg.to_csv(snakemake.output.Normal_input_df, index=False, sep='\t')
normal_ip_avg.to_csv(snakemake.output.Normal_IP_df, index=False, sep='\t')
h2o2_input_avg.to_csv(snakemake.output.H2O2_input_df, index=False, sep='\t')
h2o2_ip_avg.to_csv(snakemake.output.H2O2_IP_df, index=False, sep='\t')
heat_input_avg.to_csv(snakemake.output.Heat_input_df, index=False, sep='\t')
heat_ip_avg.to_csv(snakemake.output.Heat_IP_df, index=False, sep='\t')
stat_input_avg.to_csv(snakemake.output.Stat_input_df, index=False, sep='\t')
stat_ip_avg.to_csv(snakemake.output.Stat_IP_df, index=False, sep='\t')
yas99_avg.to_csv(snakemake.output.yAS99_df, index=False, sep='\t')
yas113_avg.to_csv(snakemake.output.yAS113_df, index=False, sep='\t')

