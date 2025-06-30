#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=06:00:00

# queue
#SBATCH --qos=normal

# memory (MB)
#SBATCH --mem=80G

# job name
#SBATCH --job-name spyce-kmer

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

#######
# env #
#######
source ~/.bashrc

export PATH="/users/asebe/aelek/bin/miniconda3/bin:$PATH"
source activate spyce

#######
# job #
#######
spyce-create \
      --setup_file "/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/sPYce/spyce.embed" \
      --k 6 \
      --norm "none" \
      --n_policy "remove" \
      --verbosity 2 \
      --data_path "/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/sPYce/" \
      --fig_path "/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/plots/sPYce/" \
      --save_prefix "k6" \
      --correct_gc \
      --n_jobs 24
