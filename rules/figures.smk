import os

include: "df_formatting.smk"
include: "ratio_stress_normal.smk"
include: "fastq_extraction.smk"

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

rule heatmap_read_type:
    """ Create one heatmap per average condition to show the number of reads per
        read type and tRNA isotype."""
    input:
        dfs = rules.average_normalized_tables.output
    output:
        heatmap = expand(os.path.join(config['figure']['heatmap'], '{cond}_read_type_per_isotype.svg'), 
			cond=['Normal_input', 'Normal_IP', 'H2O2_input', 'H2O2_IP', 'Heat_input', 
				'Heat_IP', 'Stat_input', 'Stat_IP', 'yAS99', 'yAS113'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/heatmap_read_type.py"

rule heatmap_read_type_ratio_ip_input:
    """ Create one heatmap to show the ratio of number of reads per
        read type and tRNA isotype between IPs and input average samples."""
    input:
        dfs = rules.average_normalized_tables.output
    output:
        heatmap = expand(os.path.join(config['figure']['heatmap'], 'read_type_ratio_{cond}.svg'),
                        cond=['Normal', 'H2O2', 'Heat', 'Stat', 'yAS113_yAS99']) 
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/heatmap_read_type_ratio_ip_input.py"

