#!/bin/bash

module load nixpkgs/16.09
module load gcc/5.4.0
module load kentutils/20180716

# Convert gtf file into GenePred format, then bigGenePred format (required format for Assembly hub display on UCSC genome browser)

gtf_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/references/s_pombe.gtf
index_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/references
output_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/references/s_pombe.bb
chr_size=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/star_index/chrNameLength.txt

gtfToGenePred $gtf_path temp.genePred

genePredToBigGenePred temp.genePred temp.txt

LC_COLLATE=C sort -k1,1 -k2,2n temp.txt > temp2.txt
wget https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as

bedToBigBed -extraIndex=name,geneName -type=bed12+8 -tab -as=bigGenePred.as temp2.txt $chr_size $output_path

# Create search index (so we can search by gene name in the Assembly hub)
cat temp.genePred | awk '{print $1, $12, $1}' > temp_input.txt
ixIxx temp_input.txt $index_path/search_index.ix $index_path/search_index.ixx


rm bigGenePred.as
rm temp*
echo "All done!"
