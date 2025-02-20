#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 4G
#SBATCH --mail-type ALL
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -D /home/n-z/ts16/DeletionWT_BulkRNAseq_210225
#SBATCH -J uncom

tar -xvf /home/n-z/ts16/raw-seq-data/Ceman_RNASeq.2021223.tar.bz2
