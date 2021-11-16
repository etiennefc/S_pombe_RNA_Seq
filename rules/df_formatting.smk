import os

include: "DESeq2.smk"

rule add_gene_biotype:
    """ Create a gene_biotype df for all genes in a gtf."""
    input:
        gtf = config['path']['genome_gtf'],
    output:
        gene_biotype_df = config['path']['gene_biotype_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_gene_biotype.py"

rule format_deseq_output_wt_ko:
    """ Format DESeq2 output to keep only genes with non-null values for log2
        fold change and p-value."""
    input:
        WT_vs_KO = os.path.join(rules.DESeq2_sla1_WT_KO.output.results, 'knockout-wild_type.csv')
    output:
        WT_vs_KO_clean = os.path.join(config['path']['deseq_output'],
                        'sla1_WT_KO/knockout-wild_type_v2.tsv'),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/format_deseq_output_wt_ko.py"

rule format_deseq_output_stress_normal:
    """ Format DESeq2 output to keep only genes with non-null values for log2
        fold change and p-value."""
    input:
        stress_normal_deseq_dir = rules.DESeq2_stress_vs_normal_ip.output.results
    output:
        output_dir = directory(os.path.join(config['path']['deseq_output'],
                        'sla1_IP_stress_normal/'))
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/format_deseq_output_stress_normal.py"

rule format_deseq_output_input_stress_normal:
    """ Format DESeq2 output to keep only genes with non-null values for log2
        fold change and p-value."""
    input:
        stress_normal_deseq_dir = rules.DESeq2_stress_vs_normal_input.output.results
    output:
        output_dir = directory(os.path.join(config['path']['deseq_output'],
                        'sla1_input_stress_normal/'))
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/format_deseq_output_stress_normal.py"
