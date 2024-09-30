#!/bin/sh
#SBATCH --job-name=many_topGO
#SBATCH --partition=long
#SBATCH --ntasks-per-node=150
#SBATCH --mem=500G
#SBATCH --error=%j_topgosorghum.err
#SBATCH --output=%j_topgosorghum.out

Rscript topGO.r
