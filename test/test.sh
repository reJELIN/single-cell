#!/bin/bash

########################################################################
## Single-cell RNA-seq script to launch single-cell RNA-seq pipeline
##
## using: sbatch /mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell/test/test.sh
########################################################################

## JOB PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#SBATCH --job-name=sc_lr_pipeline
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=430M
#SBATCH --partition=mediumq

source /mnt/beegfs/userdata/r_jelin/conda_install/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

module load singularity

#parameters
path_to_configfile="/mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell/test/test.yaml"
path_to_pipeline="/mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell"

#launch
snakemake --profile /mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell/profiles/slurm -s ${path_to_pipeline}/Snakefile_old --configfile ${path_to_configfile}

conda deactivate