# Train gkmSVM models on cell type level peaks,
# score kmers using these models, and derive PWMs
# https://www.beerlab.org/gkmsvm/
# note: slower but it can be used to score kmers
# 
# example usage:
# nohup bash gkmSVM_train_kmers.sh cnidocyte > gkmSVM_train_kmers_cnidocyte.out 2>&1 &

function train_gkSVM {
    pth=$1
    pfx=$2
    dir=$3
    k=$4
    posfile=${dir}/${pfx}.fg.fasta
    negfile=${dir}/${pfx}.bg.fasta
    kernel=${dir}/${pfx}.kernel.out
    svmtrain=${dir}/${pfx}
    nrkmers=${dir}/${pfx}.nr${k}mers.fa
    nrscore=${dir}/${pfx}.nr${k}mers_scores.txt
    pwm=${dir}/${pfx}.pwm
    
    # training
    echo "Building kernel for" ${pfx} "in" ${dir}
    gkmsvm_kernel -d 3 ${posfile} ${negfile} ${kernel}
    echo "Training model for" ${pfx} "in" ${dir}
    gkmsvm_train ${kernel} ${posfile} ${negfile} ${svmtrain}
    
    # alternatively, training with CV
    # echo "Training model for" ${pfx} "in" ${dir}
    # python ${pth}/scripts/cksvm_train.py ${kernel} ${positives_fa} ${negatives_fa} ${svmtrain}
    # Rscript ${pth}/scripts/rocprcurve.R ${svmtrain}_cvpred.out ${svmtrain}_rocpr.pdf
    
    # PWMs
    echo "Generating "${k}"mers for" ${pfx} "in" ${nrkmers}
    python2 ${pth}/scripts/nrkmers.py ${k} ${nrkmers}
    echo "Classifying "${k}"mers for" ${pfx} "in" ${nrkmers}
    gkmsvm_classify -d 3 ${nrkmers} ${svmtrain}_svseq.fa ${svmtrain}_svalpha.out ${nrscore}
    echo "Generating PWM for" ${pfx} "from" ${nrscore}
    python2 ${pth}/scripts/svmw_emalign.py -n 30 -c 3 ${nrscore} ${k} ${pwm}
    echo "Done "${pfx}
}

# conda activate gkmsvm
gkmsvm_path=/home/anamaria/bin/gkmsvm
gkmsvm_dir=/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/gkmSVM/
cd ${gkmsvm_dir}
ct=${1}
train_gkSVM ${gkmsvm_path} ${ct} ${gkmsvm_dir} 10
