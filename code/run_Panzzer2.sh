#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --array=1-50
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6gb
#SBATCH --job-name=PANNZER2

SUFFIX=$(ls /home/dmpachon/jorge/comparative_cane/data/panzzer/parts/ | head -n ${SLURM_ARRAY_TASK_ID} | tail -n1)
PATH="/home/dmpachon/jorge/comparative_cane/data/panzzer/parts/"
FILE=${PATH}${SUFFIX}
PANZZER="/home/dmpachon/TO_COPY_IN_OTHER_PLACE/tests_sugarcane_pantranscriptome_quantifications/panzzer2/SANSPANZ.3/runsanspanz.py"
SPECIE="Sorghum"
RESULTS_DIR="/home/dmpachon/jorge/comparative_cane/results/panzzer"

#/home/dmpachon/miniconda3/condabin/conda activate panzzer

source ~/.bashrc
conda activate panzzer
python ${PANZZER} -R -o ",${RESULTS_DIR}/${SUFFIX}_DE.out,${RESULTS_DIR}/${SUFFIX}_GO.out,${RESULTS_DIR}/${SUFFIX}_anno.out" -s "${SPECIE}" < ${FILE}
