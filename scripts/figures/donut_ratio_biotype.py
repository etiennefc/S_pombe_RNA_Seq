#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


ratio_df = pd.read_csv(snakemake.input.ratio_df, sep='\t')
biotype_df = pd.read_csv(snakemake.input.gene_biotype_df, sep='\t')

# Merge biotype to ratio_df
ratio_df = ratio_df.merge(biotype_df, how='left', on='gene_id')
ratio_df['gene_biotype'] = ratio_df['gene_biotype'].str.lstrip()  # remove white space before string

# Create a ratio category column for each stress ratio
for col in ['H2O2_avg_ratio_vs_normal', 'Heat_avg_ratio_vs_normal', 'Stat_avg_ratio_vs_normal']:
    ratio_df.loc[ratio_df[col] < 0.5, f'{col}_category'] = '<0.5'
    ratio_df.loc[(ratio_df[col] >= 0.5) & (ratio_df[col] <= 2), f'{col}_category'] = '[0.5; 2]'
    ratio_df.loc[ratio_df[col] > 2, f'{col}_category'] = '>2'



# Create a donut chart of the ratio category of stress ratio to normal ratio (outer donut)
# The inner donut shows the gene biotype
# First for H2O2 stress
count_datasets = []
count_attributes = []
for ratio in snakemake.params.ratio_colors.keys():  # Iterate through ratio category
    temp_df = ratio_df[ratio_df['H2O2_avg_ratio_vs_normal_category'] == ratio]
    count_datasets.append(len(temp_df))
    attributes_dict = {}
    for type in snakemake.params.biotype_colors.keys():
        attributes_dict[type] = len(temp_df[temp_df['gene_biotype'] == type])
    sorted_dict = {k: v for k, v in sorted(attributes_dict.items())}  # Sort dictionary alphabetically as in the config file
    print(ratio, sorted_dict)
    for val in sorted_dict.values():
        count_attributes.append(val)

counts_h2o2 = [count_datasets, count_attributes]

# Create a donut chart of the ratio category of stress ratio to normal ratio (outer donut)
# The inner donut shows the gene biotype
# Secondly, on heat stress
count_datasets2 = []
count_attributes2 = []
for ratio in snakemake.params.ratio_colors.keys():  # Iterate through ratio category
    temp_df = ratio_df[ratio_df['Heat_avg_ratio_vs_normal_category'] == ratio]
    count_datasets2.append(len(temp_df))
    attributes_dict = {}
    for type in snakemake.params.biotype_colors.keys():
        attributes_dict[type] = len(temp_df[temp_df['gene_biotype'] == type])
    sorted_dict = {k: v for k, v in sorted(attributes_dict.items())}  # Sort dictionary alphabetically as in the config file
    print(ratio, sorted_dict)
    for val in sorted_dict.values():
        count_attributes2.append(val)

counts_heat = [count_datasets2, count_attributes2]

# Create a donut chart of the ratio category of stress ratio to normal ratio (outer donut)
# The inner donut shows the gene biotype
# Thirdly, on stat stress
count_datasets3 = []
count_attributes3 = []
for ratio in snakemake.params.ratio_colors.keys():  # Iterate through ratio category
    temp_df = ratio_df[ratio_df['Stat_avg_ratio_vs_normal_category'] == ratio]
    count_datasets3.append(len(temp_df))
    attributes_dict = {}
    for type in snakemake.params.biotype_colors.keys():
        attributes_dict[type] = len(temp_df[temp_df['gene_biotype'] == type])
    sorted_dict = {k: v for k, v in sorted(attributes_dict.items())}  # Sort dictionary alphabetically as in the config file
    print(ratio, sorted_dict)
    for val in sorted_dict.values():
        count_attributes3.append(val)

counts_stat = [count_datasets3, count_attributes3]


print(counts_h2o2, counts_heat, counts_stat)

# Set inner_labels as a list of empty strings, and labels as outer and inner_labels
inner_labels = [None] * len(snakemake.params.ratio_colors.keys()) * len(snakemake.params.biotype_colors.keys())
labels = [list(snakemake.params.ratio_colors.keys()), inner_labels]

# Set inner colors as a repeated list of colors for each part of the inner donut (ex: same 2 inner colors repeated for each outer donut part)
inner_colors = list(snakemake.params.biotype_colors.values()) * len(snakemake.params.ratio_colors.keys())

# Set colors as outer and inner_colors
colors = [list(snakemake.params.ratio_colors.values()), inner_colors]

# Create donut chart for each stress
ft.donut_2(counts_h2o2, labels, colors, 'H2O2', list(snakemake.params.biotype_colors.keys()), list(snakemake.params.biotype_colors.values()), snakemake.output.donut_h2o2)
ft.donut_2(counts_heat, labels, colors, 'Heat', list(snakemake.params.biotype_colors.keys()), list(snakemake.params.biotype_colors.values()), snakemake.output.donut_heat)
ft.donut_2(counts_stat, labels, colors, 'Stat', list(snakemake.params.biotype_colors.keys()), list(snakemake.params.biotype_colors.values()), snakemake.output.donut_stat)
