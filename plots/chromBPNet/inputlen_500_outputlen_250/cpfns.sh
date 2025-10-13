for i in chrombpnet_model/*_inputlen_500_outputlen_250/evaluation/chrombpnet_only_peaks.counts_pearsonr.png
do
    dr=$(basename $(dirname $(dirname $i)))
    fn=${dr}_counts_pearsonr.png
    echo $fn
    cp $i plots/${fn}
done