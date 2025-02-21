# Train gkmSVM models on cell type level differential peaks
# using LS-GKM + gkmexplain implementation
# https://github.com/kundajelab/lsgkm
# note: it does not precompute kernel matrix 
# and cannot be used for scoring kmers
#
# input args: 
#   1. cell type
#   2. gkm input directory (it should contain ${cell_type}.fg.fasta and ${cell_type}.bg.fasta files)
#
# output:
# 
#
# example usage:
# nohup bash gkmSVM_train.sh cnidocyte gkm_input_dir > gkmSVM_train_cnidocyte.out 2>&1 &

function train_gkSVM {
    pfx=$1
    dir=$2
    posfile=${dir}/${pfx}.fg.fasta
    negfile=${dir}/${pfx}.bg.fasta
    echo "Training model for" ${pfx} "in" ${dir}
    gkmtrain -m 4000 -T 4 ${posfile} ${negfile} ${pfx}
    echo "Training finished for" ${pfx} "in" ${dir}
    echo "CV for" ${pfx} "in" ${dir}
    gkmtrain -m 4000 -T 4 -x 5 ${posfile} ${negfile} ${pfx}
    echo "CV finished for" ${pfx} "in" ${dir}
}

ct=${1}
gkmdir=${2}
cd ${gkmdir}
train_gkSVM ${ct} ${gkmdir}
