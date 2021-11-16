#!/usr/bin/python3
from snakemake.io import expand
import os
def get_figures_path(config):
    """Return list of figures (and their path) to generate from config"""
    files = []
    files.append(os.path.join(config['figure']['donut'], 'Heat_vs_normal_biotype.svg'))
    files.append(os.path.join(config['figure']['donut'], 'H2O2_vs_normal_biotype.svg'))
    files.append(os.path.join(config['figure']['donut'], 'Stat_vs_normal_biotype.svg'))

    print(files)
    return files
