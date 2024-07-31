#!/bin/sh
#SBATCH --job-name=snake_master
#SBATCH --partition=long
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --error=logs/%j_snakemaster.err
#SBATCH --output=logs/%j_snakemaster.out

#conda activate snakemake
snakemake -p -s Snakefile_download2.py --resources load=100 --cluster "sbatch --job-name {params.jobname} --partition {params.partition} --ntasks-per-node {threads} --mem {params.mem} --error {log} --output {log}" --jobs 10 --keep-going --conda-frontend mamba --use-conda --rerun-triggers input

