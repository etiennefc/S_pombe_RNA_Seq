import os

rule DESeq2_sla1_WT_KO:
    """ Differential expression per gene for the different conditions (sla1 WT
        vs KO). The samples in design_wt_ko.tsv must match exactly the keys of
        the 'dataset' dictionary located in the config.json file."""
    input:
        counts = os.path.join(config["path"]["coco_merge"], "merged_counts.tsv"),
        samples = "data/design_wt_ko.tsv"
    output:
        results = directory("results/DESeq2/sla1_WT_KO"),
    log:
        "logs/DESeq2/sla1_WT_KO.log"
    conda:
        "../envs/deseq2_new.yaml"
    script:
        "../scripts/DESeq2_sla1_WT_KO.R"

rule DESeq2_stress_vs_normal_input:
    """ Differential expression per gene for the different conditions (input
        under various stresses vs input in normal conditions. The samples in
        design_stress_vs_normal_input.tsv must match exactly the keys of the
        'dataset' dictionary located in the config.json file."""
    input:
        counts = os.path.join(config["path"]["coco_merge"], "merged_counts.tsv"),
        samples = "data/design_stress_vs_normal_input.tsv"
    output:
        results = directory("results/DESeq2/sla1_input_stress_normal"),
    log:
        "logs/DESeq2/sla1_input_stress_normal.log"
    conda:
        "../envs/deseq2_new.yaml"
    script:
        "../scripts/DESeq2_sla1_input_stress_normal.R"

rule DESeq2_stress_vs_normal_ip:
    """ Differential expression per gene for the different conditions (sla1 IP
        under various stresses vs sla1 IP in normal conditions and also stress
        IP vs other stress IP). The samples in design_stress_vs_normal_ip.tsv must
        match exactly the keys of the 'dataset' dictionary located in the config.json
        file."""
    input:
        counts = os.path.join(config["path"]["coco_merge"], "merged_counts.tsv"),
        samples = "data/design_stress_vs_normal_ip.tsv"
    output:
        results = directory("results/DESeq2/sla1_IP_stress_normal"),
    log:
        "logs/DESeq2/sla1_IP_stress_normal.log"
    conda:
        "../envs/deseq2_new.yaml"
    script:
        "../scripts/DESeq2_sla1_IP_stress_normal.R"
