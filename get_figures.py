#!/usr/bin/python3
from snakemake.io import expand
import os
def get_figures_path(config):
    """Return list of figures (and their path) to generate from config"""
    files = []
    files.append(os.path.join(config['figure']['donut'], 'Heat_vs_normal_biotype.svg'))
    files.append(os.path.join(config['figure']['donut'], 'H2O2_vs_normal_biotype.svg'))
    files.append(os.path.join(config['figure']['donut'], 'Stat_vs_normal_biotype.svg'))
    files.extend(expand(os.path.join(config['figure']['heatmap'], '{cond}_read_type_per_isotype.svg'), 
                        cond=['Normal_input', 'Normal_IP', 'H2O2_input', 'H2O2_IP', 'Heat_input', 
                                'Heat_IP', 'Stat_input', 'Stat_IP', 'yAS99', 'yAS113']))
    files.extend(expand(os.path.join(config['figure']['heatmap'], 'read_type_ratio_{cond}.svg'),
                        cond=['Normal', 'H2O2', 'Heat', 'Stat', 'yAS113_yAS99']))

    print(files)
    return files
