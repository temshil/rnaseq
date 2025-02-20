#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 4G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/DeletionWT_BulkRNAseq_210225
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J multiqc

module load MultiQC/1.9-IGB-gcc-8.2.0-Python-3.7.2

multiqc /home/n-z/ts16/DeletionWT_BulkRNAseq_210225
