#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --ntasks-per-node=100
#SBATCH --mem=300G
#SBATCH --job-name=big_histogram

python histogram_cor.py

