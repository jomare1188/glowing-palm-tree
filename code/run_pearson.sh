#!/bin/sh
#SBATCH --job-name=big_pearson_cor
#SBATCH --partition=long
#SBATCH --ntasks-per-node=10
#SBATCH --mem=300G
#SBATCH --error=%j_cor.err
#SBATCH --output=%j_cor.out

Rscript pearson_cor.r
