#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --ntasks-per-node=50
#SBATCH --mem=80G
#SBATCH --job-name=orthofinder

orthofinder -I 2 -t 50 -a 20 -f ../data/proteins/ -o ../results/orthofinder/
