#!/bin/sh
#SBATCH --job-name=snake_master
#SBATCH --partition=long
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --error=logs/%j_snakemaster.err
#SBATCH --output=logs/%j_snakemaster.out
#SBATCH --nodelist=n01

#conda activate snakemake
snakemake -p -s snakefile_preprocessing.py --resources load=100 --cluster "sbatch --nodelist=n01 --job-name {params.jobname} --partition {params.partition} --ntasks-per-node {threads} --mem {params.mem} --error {log} --output {log}" --jobs 3 --keep-going --conda-frontend mamba --use-conda --rerun-incomplete --rerun-triggers input
