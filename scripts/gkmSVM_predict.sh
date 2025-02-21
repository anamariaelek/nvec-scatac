
function predict_gkSVM {
    model_file=$1
    test_seqfile=$2
    output_file=$3
    echo "Prediction for" ${outprefix} "using" ${model_file}
    gkmpredict ${test_seqfile} ${model_file} ${output_file}
}

moddir=${1}
modelprefix=${2}
fastafile=${3}
outprefix=${fastafile%%.fasta}
predict_gkSVM \
    ${moddir}/${modelprefix}.model.txt \
    ${moddir}/${fastafile} \
    ${moddir}/${modelprefix}.${outprefix}.gkmpredict.txt