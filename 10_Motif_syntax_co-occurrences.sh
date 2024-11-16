#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

# time limit 
#SBATCH --time=8:00:00

# queue
#SBATCH --qos=normal

# memory (MB)
#SBATCH --mem=10G

# job name
#SBATCH --job-name motifs

# job array directive
#SBATCH --array=0-112

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail

#######
# env #
#######
#source ~/.bashrc
#source /users/asebe/aelek/bin/miniconda2/etc/profile.d/conda.sh
# eval "$(conda shell.bash hook)"

###############################################
# get the fastq file based on the array index #
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################
WDIR=`dirname $PWD`
BIN=$((SLURM_ARRAY_TASK_ID+1)) 

Rscript 10_Motif_syntax_co-occurrences.R ${BIN} 
