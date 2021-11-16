import os

rule ratio_tables:
    """ From the abundance (in TPM) merged table, create a new table of the
        ratios between IP and input. Add also the ratio of ratios to compare with
        normal samples ((Stress IP/Stress input)/(Normal IP/Normal input))."""
    input:
        tpm_df = os.path.join(config['path']['coco_merge'], "merged_tpm.tsv")
    output:
        ratio_df = os.path.join(config['path']['ratio_stress_normal'], "ratio_stress_normal.tsv")
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ratio_tables.py"
