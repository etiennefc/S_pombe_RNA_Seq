#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=32G
#SBATCH --array=[0-31]

module load nixpkgs/16.09
module load picard/2.18.9

namelist=('decoy' \
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
project_path=$SCRATCH/s_pombe_rna_seq/s_pombe_rna_seq/results/

mkdir -p $project_path/Picard/$name

BAM_FILE=$project_path/star/$name/Aligned.sortedByCoord.out.bam

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
      I=$BAM_FILE \
      O=$project_path/Picard/$name/Picard_insert_size_metrics.txt \
      H=$project_path/Picard/$name/Picard_size_histogram.pdf \
