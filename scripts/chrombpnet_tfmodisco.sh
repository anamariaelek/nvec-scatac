#!/usr/bin/env bash
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem 20G
#SBATCH --time 8:00:00
#SBATCH -A lp_cbd_stae
#SBATCH -p bigmem
#SBATCH --cluster wice
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user anamariaelek@gmail.com
#SBATCH --output logs/%x.%j.out
#SBATCH --error logs/%x.%j.err
#SBATCH --job-name chrombpnet_tfmodisco

set -e
set -o pipefail

# environment
source /staging/leuven/stg_00002/conda/vsc36173/etc/profile.d/conda.sh
conda activate chrombpnet


# hard coded file paths
wdir=/staging/leuven/stg_00002/lcb/aelek/chrombpnet_nvec/
genome=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Nvec_vc1.1_gDNA.fasta
peaks=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Peaks_cell_type_mapped_250bp.bed
chromsizes=/staging/leuven/stg_00002/lcb/aelek/data_nvec/Nvec_vc1.1_gDNA.chrom.sizes

# path to chrombpnet code
cbp="/staging/leuven/stg_00002/lcb/aelek/software/chrombpnet/chrombpnet/"

# default params
meme_motifs=${cbp}"/data/motifs.meme.txt"
max_seqlets=1000000
window=250
sliding_window=20

# function to display usage
usage() {
    echo "
    Usage: $0 name out_prefix [-s <one-hot encoded sequences>] [-a <SHAP hypothetical scores>] [-i <h5 file with prediction scores>] [-m <motifs file in meme format>] [-n <max seqlets>] [-w <window>]

    -s one-hot encoded sequences
	-a SHAP hypothetical scores
	-i h5 file with prediction scores, e.g.
		${model_dir}/predictions/${name}.profile_scores.h5 
		${model_dir}/predictions/${name}.counts_scores.h5
	-m motifs file in meme format, if not supplied, default file at
		chrombpnet/data/motifs.meme.txt is used
	-n max seqlets, see tfmodisco -h
	-w window, see tfmodisco -h

	Either -s and -a, or -i must be provided.

    "
    exit 1
}

# parse positional arguments
name=$1
out_prefix=$2

if [ $# -lt 2 ]; then
    usage
fi

shift 2

# parse optional arguments
while getopts "s:a:i:m:n:w:z:" opt; do
    case ${opt} in
	s ) ohe=$OPTARG ;;
	a ) shap=$OPTARG ;;
	i ) h5=$OPTARG ;;
        m ) meme_motifs=$OPTARG ;;
        n ) max_seqlets=$OPTARG ;;
        w ) window=$OPTARG ;;
	z ) sliding_window=$OPTARG ;;
        \? ) usage ;;
    esac
done

# check for mandatory arguments
if { [ -z "$ohe" ] || [ -z "$shap" ]; } && [ -z "$h5" ]; then
    echo "Error: Either -s and -a, or -i must be provided."
    usage
fi

if [ -n "$ohe" ] && [ -z "$shap" ]; then
    echo "Error: -a is required when -s is provided."
    usage
fi

if [ -z "$ohe" ] && [ -n "$shap" ]; then
    echo "Error: -s is required when -a is provided."
    usage
fi

# function to get the current time with the new line character
function timestamp {
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# tfmodisco output
motifs_dir=$( dirname $out_prefix )
if [ ! -d "${motifs_dir}" ]; then
    mkdir ${motifs_dir}
fi
motifs_out=${out_prefix}.modisco_results.h5
report_out=${out_prefix}.modisco_report

# log file
logfile=${motifs_dir}"/logs/chrombpnet_tfmodisco_$( timestamp ).log"
if [ ! -d "${motifs_dir}/logs" ]; then
    mkdir ${motifs_dir}/logs
fi
touch ${logfile}

# run TFModisco
echo $( timestamp ): "Running tfmodisco-lite for "${name} | tee -a ${logfile}
echo "Params:"
echo "    -w" ${window}
echo "    -z" ${sliding_window}
echo "    -n" ${max_seqlets}

if [ -n "${h5}" ]; then
	echo "Input: "${h5}
	modisco motifs -i ${h5} -n ${max_seqlets} -o ${motifs_out} -w ${window} -z ${sliding_window} | tee -a ${logfile}
fi

if [ -n "$ohe" ] && [ -n "$shap" ]; then
	echo "Input: "${ohe}" & "${shap}
	modisco motifs -s ${ohe} -a ${shap} -n ${max_seqlets} -o ${motifs_out} -w ${window} -z ${sliding_window} | tee -a ${logfile}
fi

# reports
echo $( timestamp ): "Generating reports for "${name} | tee -a ${logfile}
modisco report -i ${motifs_out} -o ${report_out} -m ${meme_motifs} | tee -a ${logfile}

echo $( timestamp ): "Done tfmodisco-lite for "${name} | tee -a ${logfile}
echo "All done." | tee -a ${logfile}
