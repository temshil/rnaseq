#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 50G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/DeletionWT_BulkRNAseq_210225
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J salmon-quant
#SBATCH --array 1-6

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" src/filenames.txt)

module load Salmon/1.4.0-IGB-gcc-8.2.0

salmon quant -i salmon/salmon_index -l A \
-1 data/${line}_1.fastq.gz \
-2 data/${line}_2.fastq.gz \
-p 12 --seqBias --gcBias --numBootstraps=30 --recoverOrphans --validateMappings -o results/${line}
