#!/bin/bash
#SBATCH -J run_sim_inference  # Job name
#SBATCH -o scripts/inference-scripts/logs/run_sim_inference_%A.out  # Output log
#SBATCH -p sapphire           # Partition (queue) name
#SBATCH -c 8                  # Number of CPU cores
#SBATCH --mem=16G             # Memory allocation (increase if needed)
#SBATCH -t 08:00:00           # Max runtime (hours:minutes:seconds)

# Define R library path and Singularity container
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18

rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

# Singularity execution command
singularity_command="singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image}"

# Run the R script
$singularity_command Rscript -e "source('/n/holylabs/LABS/pillai_lab/Users/xmeng1/scmatch2/scripts/inference-scripts/1_sim_inference_main.R')"
