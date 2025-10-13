#!/usr/bin/env bash
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 50G
#SBATCH --time 12:00:00
#SBATCH -A lp_big_wice_gpu
#SBATCH -p dedicated_big_gpu
#SBATCH --gpus-per-node 1
#SBATCH --cluster wice
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user anamariaelek@gmail.com
#SBATCH --output logs/%x.%j.out
#SBATCH --error logs/%x.%j.err
#SBATCH --job-name pipeline.500.250.512.4

set -e
set -o pipefail

source /staging/leuven/stg_00002/conda/vsc36173/etc/profile.d/conda.sh
conda activate chrombpnet
module load cuDNN/8.2.1.32-CUDA-11.3.1

# file paths
wdir=/staging/leuven/stg_00002/lcb/aelek/chrombpnet_nvec/
genome=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Nvec_vc1.1_gDNA.fasta
peaks=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Peaks_cell_type_mapped_250bp.bed
chromsizes=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Nvec_vc1.1_gDNA.chrom.sizes

# params
inlen=500
outlen=250
filters=512
dilation_layers=4

# bias model to train
bias_name=neuron_Pou4_FoxL2_2
bias_model=${wdir}/bias_model/${bias_name}_inputlen_${inlen}_outputlen_${outlen}/models/bias.h5

# cell type to train chrombpnet model for
name=neuron_Pou4_FoxL2_2
model_dir=${wdir}/chrombpnet_model/${name}_inputlen_${inlen}_outputlen_${outlen}

# function to get the current time with the new line character
function timestamp {
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# make nonopeaks bed if it doesn't exist
if [ ! -f "${wdir}/data/nvec_${inlen}_negatives.bed" ]; then
    echo $( timestamp ): "Making non-peaks bed file with input width ${inlen}."
	chrombpnet prep nonpeaks \
	    -g ${genome} \
	    -p ${peaks} \
	    -c ${chromsizes} \
	    -st 100 \
	    -il ${inlen} \
	    -npr 4 \
	    -fl ${wdir}/data/splits.json \
	    -o ${wdir}/data/nvec_${inlen}
fi

# bias model training if doesn't already exist
if [ ! -f "${bias_model}" ]; then
    echo $( timestamp ): "Training bias model for ${bias_name}."
	chrombpnet bias pipeline \
	    -ibam /staging/leuven/stg_00002/lcb/aelek/data_nvec/GroupBam/${bias_name}.bam \
	    -d "ATAC" \
	    -g ${genome} \
	    -c ${chromsizes} \
	    -p ${peaks} \
	    -n ${wdir}/data/nvec_${inlen}_negatives.bed \
	    -fl ${wdir}/data/splits.json \
	    -b 0.5 \
	    -o ${wdir}/bias_model/${bias_name}_inputlen_${inlen}_outputlen_${outlen} \
	    --inputlen ${inlen} \
	    --outputlen ${outlen} \
	    --filters ${filters} \
	    --n-dilation-layers ${dilation_layers}
	wait
	echo $( timestamp ): "Trained bias model for ${bias_name}. Moving on."
fi

# model training
if [ ! -d  ${model_dir} ]; then
    echo $( timestamp ): "Training the model for ${name}."
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
        --n-dilation-layers ${dilation_layers}

    echo $( timestamp ): "Trained bias-factorized chrombpnet model."
fi