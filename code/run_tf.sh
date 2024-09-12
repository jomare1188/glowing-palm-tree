#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --array=1-20%10 
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=8gb
#SBATCH --job-name=hmmr_tf



export m=${SLURM_ARRAY_TASK_ID}
species=sugarcane
db=../data/tf_annot/db/37v_Pfam-A.hmm.gz

hmmsearch --seed 999 --cpu 5 --domtblout ../results/tf_annot/${species}/domtbl.${m}.out --cut_ga -o ./../results/tf_annot/${species}/all.${m}.out ${db} ../data/proteins/parts_sugarcane/part_${m}_sugarcane_proteins_of_longest_cds_per_OG.fasta
