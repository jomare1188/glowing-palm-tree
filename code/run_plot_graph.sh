#!/bin/bash
#SBATCH --job-name=ok_plot
#SBATCH --partition=long
#SBATCH --ntasks-per-node=10
#SBATCH --mem=300G
#SBATCH --error=%j_plotN.err
#SBATCH --output=%j_plotN.out

Rscript plot_network.r 
