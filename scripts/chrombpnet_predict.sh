#!/usr/bin/env bash
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 20G
#SBATCH --time 12:00:00
#SBATCH -A lp_big_wice_gpu
#SBATCH -p dedicated_big_gpu
#SBATCH --gpus-per-node 1
#SBATCH --cluster wice
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user anamariaelek@gmail.com
#SBATCH --output logs/%x.%j.out
#SBATCH --error logs/%x.%j.err
#SBATCH --job-name chrombpnet_predict

set -e
set -o pipefail

# environment
source /staging/leuven/stg_00002/conda/vsc36173/etc/profile.d/conda.sh
conda activate chrombpnet
module load cuDNN/8.2.1.32-CUDA-11.3.1

# hard coded file paths
wdir=/staging/leuven/stg_00002/lcb/aelek/chrombpnet_nvec/
genome=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Nvec_vc1.1_gDNA.fasta
peaks=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Peaks_cell_type_mapped_250bp.bed
chromsizes=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Nvec_vc1.1_gDNA.chrom.sizes

# function to display usage
usage() {
    echo "
    Usage: $0 model_file bigwig_file

    model_file	Bias-factorized model e.g. ${model_dir}/models/chrombpnet_nobias.h5
    bigwig_file Pseudobulk bigwig file with observed (groundtruth) signal.
    "
    exit 1
}

# parse positional arguments
if [ $# -lt 2 ]; then
    usage
fi

# cell type to train chrombpnet model for
model_file=$1
model_dir=$( dirname $model_file)
model_dir=$( dirname $model_dir)
name=$( basename $model_dir )
name=${name%%_inputlen_*}

# path to bigwig file
bigwig_file=$2

model_type=$( basename $( dirname $model_dir ) )
if [ $model_type == "bias_model" ]; then
	name=${name}_bias
fi

# function to get the current time with the new line character
function timestamp {
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# log file
logfile=${model_dir}"/logs/chrombpnet_predict_$( timestamp ).log"
touch ${logfile}

# run predictions
pred_dir=${model_dir}/predictions
if [ ! -d ${pred_dir} ]; then
    mkdir -p ${pred_dir}
fi

echo $( timestamp ): "Starting predictions for ${model_file}" | tee -a ${logfile}

chrombpnet contribs_bw \
    -m ${model_file} \
    -r ${peaks} \
    -g ${genome} \
    -c ${chromsizes} \
    -op ${pred_dir}/${name} | tee -a ${logfile}

echo $( timestamp ): "Done predictions with ${model_file}." | tee -a ${logfile}

echo $( timestamp ): "Generating bigwigs containing predictios for ${model_file} and comparing to ${bigwig_file}" | tee -a ${logfile}

chrombpnet pred_bw \
    -cmb ${model_file} \
    -r ${peaks} \
    -g ${genome} \
    -c ${chromsizes} \
    -op ${pred_dir}/${name} \
    -bw ${bigwig_file} | tee -a ${logfile}

echo $( timestamp ): "Done generating bigwigs with ${model_file}." | tee -a ${logfile}
{logfile}
echo "All done." | tee -a ${logfile}

