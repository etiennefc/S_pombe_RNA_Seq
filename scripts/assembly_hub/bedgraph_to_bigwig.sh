#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=8000M
#SBATCH --array=[0-29]

module load nixpkgs/16.09
module load gcc/5.4.0
module load kentutils/20180716



namelist=(\
'Normal_input_1' \
'Normal_input_2' \
'Normal_input_3' \
'Normal_IP_1' \
'Normal_IP_2' \
'Normal_IP_3' \
'H2O2_input_1' \
'H2O2_input_2' \
'H2O2_input_3' \
'H2O2_IP_1' \
'H2O2_IP_2' \
'H2O2_IP_3' \
'Heat_input_1' \
'Heat_input_2' \
'Heat_input_3' \
'Heat_IP_1' \
'Heat_IP_2' \
'Heat_IP_3' \
'Stat_input_1' \
'Stat_input_2' \
'Stat_input_3' \
'Stat_IP_1' \
'Stat_IP_2' \
'Stat_IP_3' \
'yAS99_1' \
'yAS99_2' \
'yAS99_3' \
'yAS113_1' \
'yAS113_2' \
'yAS113_3' \
)

name=${namelist[$SLURM_ARRAY_TASK_ID]}
bedgraph_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/results/coco_bedgraph
output_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/results/bedgraph_to_bigwig/

mkdir -p $output_path


LC_COLLATE=C sort -k1,1 -k2,2n $bedgraph_path/${name}.bedgraph > $output_path/${name}_sorted.bedGraph &&


bedGraphToBigWig $output_path/${name}_sorted.bedGraph \
$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/data/star_index/chrNameLength.txt \
$output_path/${name}.bigwig

echo 'Well done!'
