#!/bin/bash
#SBATCH -J run_ferman         # Job name
#SBATCH -o scripts/new-inference/logs/run_ferman_%A_%a.out   # output and error in one file
#SBATCH -p sapphire           # Partition (queue) name
#SBATCH -c 8              # Number of cores
#SBATCH --mem=8G          # Memory in GB
#SBATCH -t 12:00:00         # Runtime (hours:minutes:seconds)


my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"


singularity_command="singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image}"

$singularity_command Rscript -e "source('scripts/new-inference/replicate-ferman.R')"

