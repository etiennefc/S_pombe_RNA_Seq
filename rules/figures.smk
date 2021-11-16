import os

include: "df_formatting.smk"
include: "ratio_stress_normal.smk"

rule donut_ratio_biotype:
    """ Create a donut chart per stress condition ratio (relative to the normal
        IP/input ratio) to show in the outer donut the proportion of genes that
        have a ratio below 0.5, above 2 or in between. The inner donut shows the
        gene_biotype of the IPed RNAs."""
    input:
        ratio_df = rules.ratio_tables.output.ratio_df,
        gene_biotype_df = rules.add_gene_biotype.output.gene_biotype_df
    output:
        donut_heat = os.path.join(config['figure']['donut'], 'Heat_vs_normal_biotype.svg'),
        donut_h2o2 = os.path.join(config['figure']['donut'], 'H2O2_vs_normal_biotype.svg'),
        donut_stat = os.path.join(config['figure']['donut'], 'Stat_vs_normal_biotype.svg')
    params:
        biotype_colors = config['colors']['gene_biotype'],
        ratio_colors = config['colors']['ratio']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/donut_ratio_biotype.py"


rule volcano:
    """ Create a volcano plot for each comparison (WT vs KO, Stress IP vs Normal
        IP, Stress IP vs other stress)."""
    input:
        deseq_wt_ko = rules.format_deseq_output_wt_ko.output.WT_vs_KO_clean,
        deseq_ip_stress_normal_dir = rules.format_deseq_output_stress_normal.output.output_dir,
        deseq_input_stress_normal_dir = rules.format_deseq_output_input_stress_normal.output.output_dir,
        biotype_df = rules.add_gene_biotype.output.gene_biotype_df
    output:
        volcano_output_dir = directory(os.path.join(config['figure']['volcano'], 'deseq/'))
    params:
        biotype_colors = config['colors']['gene_biotype']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/volcano.py"
