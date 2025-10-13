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
#SBATCH --job-name chrombpnet_train

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

# default values of params
inlen=500
outlen=250
filters=512
dilation_layers=4

# function to display usage
usage() {
    echo "Usage: $0 name bias_name [-i inlen] [-o outlen] [-f filters] [-d dilation_layers]"
    exit 1
}

# parse positional arguments
if [ $# -lt 2 ]; then
    usage
fi

# cell type to train chrombpnet model for
name=$1

# bias model to use
bias_name=$2

shift 2

# parse optional arguments
while getopts "i:o:f:d:" opt; do
    case ${opt} in
        i ) inlen=$OPTARG ;;
        o ) outlen=$OPTARG ;;
        f ) filters=$OPTARG ;;
        d ) dilation_layers=$OPTARG ;;
        \? ) usage ;;
    esac
done

# dirs
model_dir=${wdir}/chrombpnet_model/${name}_inputlen_${inlen}_outputlen_${outlen}
bias_model=${wdir}/bias_model/${bias_name}_inputlen_${inlen}_outputlen_${outlen}/models/bias.h5

# function to get the current time with the new line character
function timestamp {
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# log file
logfile=${wdir}"/scripts/logs/chrombpnet_train_$( timestamp ).log"
touch ${logfile}

echo $( timestamp ): "Parameters:
  cell type: ${name}
  bias model: ${bias_model}
  input length: ${inlen}
  output length: ${outlen}
  filters: ${filters}
  dilation layers: ${dilation_layers}
" | tee -a ${logfile}

# make nonopeaks bed if it doesn't exist
if [ ! -f "${wdir}/data/nvec_${inlen}_negatives.bed" ]; then
    echo $( timestamp ): "Making non-peaks bed file with input width ${inlen}." | tee -a ${logfile}
	chrombpnet prep nonpeaks \
	    -g ${genome} \
	    -p ${peaks} \
	    -c ${chromsizes} \
	    -st 100 \
	    -il ${inlen} \
	    -npr 4 \
	    -fl ${wdir}/data/splits.json \
	    -o ${wdir}/data/nvec_${inlen} | tee -a ${logfile}
fi

# check if bias model exists
if [ ! -f "${bias_model}" ]; then
    echo $( timestamp ): "Bias model ${bias_model} doesn't exist! Exiting." | tee -a ${logfile}
    exit 1
fi

# model training
if [ ! -d  ${model_dir} ]; then
    echo $( timestamp ): "Training the model for ${name}." | tee -a ${logfile}

    chrombpnet pipeline \
        -ibam /staging/leuven/stg_00002/lcb/aelek/data_nvec/GroupBam/${name}.bam \
        -d "ATAC" \
        -g ${genome} \
        -c ${chromsizes} \
        -p ${peaks} \
        -n ${wdir}/data/nvec_${inlen}_negatives.bed \
        -fl ${wdir}/data/splits.json \
        -b ${bias_model} \
        --inputlen ${inlen} \
        --outputlen ${outlen} \
        -o ${model_dir} \
        --filters ${filters} \
        --n-dilation-layers ${dilation_layers} | tee -a ${logfile}

    echo $( timestamp ): "Trained bias-factorized chrombpnet model for ${name}." | tee -a ${logfile}
    logname=$( basename ${logfile} )
    mv ${logfile} ${model_dir}/logs/${logname}
else
    echo $( timestamp ): "Model directory ${model_dir} already exists! Exiting." | tee -a ${logfile}
    exit 1
fi
echo "All done!" | tee -a ${logfile}

