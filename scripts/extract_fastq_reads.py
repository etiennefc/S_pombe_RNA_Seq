#!/usr/bin/python3
import pandas as pd
import subprocess as sp

tRNA = pd.read_csv(snakemake.input.tRNA_sequences)
input_fastq = snakemake.input.fastq
output_dir = snakemake.output.fasta_dir
agrep_path = snakemake.params.agrep_path
sample_id = snakemake.wildcards.sample_ids

# Create dict of isotype_anticodon as keys and their fishing sequence(s) in a list as values
d = {}
for isotype in list(pd.unique(tRNA['Isotype'])):
    seqs = list(tRNA[tRNA['Isotype'] == isotype].sequence.values)
    d[isotype] = seqs

# For each isotype_anticodon, create a file containing all reads containing its possible corresponding fishing sequence(s) with <=2 mismatch allowed
sp.call(f'mkdir -p {output_dir}', shell=True)
temp_fastq = f'temp_{sample_id}.fastq'
sp.call(f'gunzip -c {input_fastq} > {temp_fastq}', shell=True)
for k, v in d.items():
    if len(v) == 1:  # if only one fishing sequence per isotype_anticodon
        seq = v[0]
        sp.call(f'echo {k} {seq}', shell=True)
        sp.call(f'{agrep_path} -2 {seq} {temp_fastq} | head -n -1 > {output_dir}/{k}.fa', shell=True)
        sp.call("""gawk -i inplace 'BEGIN{j=0} {j+=1; print ">read_"j"\\n"$0}' """+f"""{output_dir}/{k}.fa""", shell=True)  # add read id for all reads
    else:  # if multiple fishing sequence per isotype_anticodon
        len_v = len(v) - 1
        for i, seq in enumerate(v):
            sp.call(f'echo {k} {seq}', shell=True)
            sp.call(f'{agrep_path} -2 {seq} {temp_fastq} | head -n -1 >> {output_dir}/{k}.fa', shell=True)
            sp.call(f'wc -l {output_dir}/{k}.fa', shell=True)
            if i == len_v:  # for last iteration, add read id for all reads that match all sequences in v
                sp.call("""gawk -i inplace 'BEGIN{j=0} {j+=1; print ">read_"j"\\n"$0}' """+f"""{output_dir}/{k}.fa""", shell=True)

sp.call(f'rm {temp_fastq}', shell=True)

