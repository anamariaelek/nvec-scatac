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
#SBATCH --job-name chrombpnet_bias

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
bias_threshold=0.5
inlen=500
outlen=250
filters=512
dilation_layers=4

# function to display usage
usage() {
    echo "Usage: $0 bias_name [-i inlen] [-o outlen] [-f filters] [-d dilation_layers]"
    exit 1
}

# parse positional arguments
if [ $# -lt 1 ]; then
    usage
fi

# cell type to train bias model for
bias_name=$1

shift 1

# parse optional arguments
while getopts "b:i:o:f:d:" opt; do
    case ${opt} in
        b ) bias_threshold=$OPTARG ;;
        i ) inlen=$OPTARG ;;
        o ) outlen=$OPTARG ;;
        f ) filters=$OPTARG ;;
        d ) dilation_layers=$OPTARG ;;
        \? ) usage ;;
    esac
done

# directory to svae model to
model_dir=${wdir}/bias_model/${bias_name}_inputlen_${inlen}_outputlen_${outlen}

# function to get the current time with the new line character
function timestamp {
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# log file
logfile=${wdir}"/scripts/logs/chrombpnet_bias_$( timestamp ).log"
touch ${logfile}

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

# bias model training
if [ ! -d  ${model_dir} ]; then
    echo $( timestamp ): "Training the bias model for ${bias_name}." | tee -a ${logfile}

    chrombpnet bias pipeline \
        -ibam /staging/leuven/stg_00002/lcb/aelek/data_nvec/GroupBam/${bias_name}.bam \
        -d "ATAC" \
        -g ${genome} \
        -c ${chromsizes} \
        -p ${peaks} \
        -n ${wdir}/data/nvec_${inlen}_negatives.bed \
        -fl ${wdir}/data/splits.json \
        -b ${bias_threshold} \
        --inputlen ${inlen} \
        --outputlen ${outlen} \
        -o ${wdir}/bias_model/${bias_name}_inputlen_${inlen}_outputlen_${outlen} \
        --filters ${filters} \
        --n-dilation-layers ${dilation_layers} | tee -a ${logfile}

    echo $( timestamp ): "Trained bias model for ${bias_name}." | tee -a ${logfile}
    logname=$( basename ${logfile} )
    mv ${logfile} ${model_dir}/logs/${logname}
else
    echo $( timestamp ): "Bias model directory ${model_dir} already exists! Exiting." | tee -a ${logfile}
    exit 1
fi
echo "All done." | tee -a ${logfile}

