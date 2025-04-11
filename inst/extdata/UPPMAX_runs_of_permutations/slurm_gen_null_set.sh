#!/bin/bash -l
#SBATCH -A naiss2023-5-517
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J null_set_permutation

i=$SLURM_ARRAY_TASK_ID

singularity exec ../../r2u_container_v1.2.sif Rscript -e "rmarkdown::render('permutation_null_matrix.Rmd', params=list(arg=${i}))"

