import os


rule download_genome:
    """Download the reference genome (fasta file) used for this analysis."""
    output:
        genome = config['path']['reference_genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
	       "gunzip {output.genome}.gz"

rule download_annotation:
    """Download the annotation (gff file) used for this analysis."""
    output:
        genome_gff = config['path']['genome_gff']
    params:
        link_genome_gff = config['download']['genome_gff']
    shell:
        "wget --quiet -O {output.genome_gff}.gz {params.link_genome_gff} && "
        "gunzip {output.genome_gff}.gz"

rule download_coco_git:
    """Download git repository of CoCo."""
    output:
        git_coco_folder = directory('git_repos/coco')
    params:
        git_coco_link = config['download']['coco_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_coco_folder} '
        '&& git clone {params.git_coco_link} {output.git_coco_folder}'

rule download_agrep:
    """ Download git repository of agrep, a bash tool to to approximate (fuzzy)
        grep searches (ex: allow N mismatches in sequence search)."""
    output:
        git_agrep_folder = directory('git_repos/agrep')
    params:
        git_agrep_link = config['download']['agrep_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_agrep_folder} && '
        'git clone {params.git_agrep_link} {output.git_agrep_folder} && '
        'cd git_repos/agrep && make'
