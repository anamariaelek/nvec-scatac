#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=6:00:00

# queue
#SBATCH --qos=normal

# memory (MB)
#SBATCH --mem=80G

# job name
#SBATCH --job-name streme

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

##########
# params #
##########
beddir=${1}
name=${2}
genome=${3}
bedfg=${beddir}"/"${name}.bed
bedbg=${beddir}"/"${name}_bg.bed
fastafg=${bedfg%%bed}fasta
fastabg=${bedbg%%bed}fasta
outdir=${beddir}"/"${name}
mkdir -p ${outdir}
log=${outdir}/${name}_streme.log

#######################
# Run streme analysis #
#######################
echo "# Creating fasta file from peaks ${bedfg}" |& tee -a ${log} 
bedtools getfasta -fi ${genome} -bed ${bedfg} -fo ${fastafg}
echo "Saved ${fastafg}" |& tee -a ${log} 


echo "# Creating fasta file from peaks ${bedbg}" |& tee -a ${log} 
bedtools getfasta -fi ${genome} -bed ${bedbg} -fo ${fastabg}
echo "Saved ${fastabg}" |& tee -a ${log} 


echo "# Streme motif enrichment" |& tee -a ${log}
streme_outdir="${outdir}/streme"
mkdir -p ${streme_outdir}
streme --oc ${streme_outdir} --dna --minw 6 --maxw 18 --p ${fastafg} --n ${fastabg} |& tee -a ${log} # --nmotifs 9
echo "Streme finished. Result saved in ${streme_outdir}." |& tee -a ${log} 

#########################
# Compare to archetypes #
#########################
motifs="/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Archetypes/motif-archetypes-PPM-PCC-0.8-IC0.5-5bp-pwms.meme"
streme_motifs="${streme_outdir}/streme.xml"

echo "# Running Tomtom to compare motifs" |& tee -a ${log}
tomtom_outdir="${outdir}/tomtom"
rm -rf ${tomtom_outdir}
tomtom -o ${tomtom_outdir} -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 0.5 ${streme_motifs} ${motifs} |& tee -a ${log}

echo "Tomtom finished. Results saved in ${tomtom_outdir}" |& tee -a ${log}

#############################
# Extract top match per motif #
#############################
top_matches="${tomtom_outdir}/top_matches.txt"
awk 'NR==1 || !seen[$1]++' ${tomtom_outdir}/tomtom.tsv > ${top_matches}
echo "Saved top matches in ${top_matches}" |& tee -a ${log}
