#!/bin/bash

module load nixpkgs/16.09
module load gcc/5.4.0
module load kentutils/20180716

# Convert fasta file of the genome into TwoBit format (required format for Assembly hub display on UCSC genome browser)

genome_fasta_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/references/s_pombe_genome.fa
output_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/references/s_pombe_genome.2bit


faToTwoBit $genome_fasta_path $output_path

echo "Well done!"
