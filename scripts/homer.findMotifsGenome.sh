
genome="/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/genome/Nvec_vc1.1_gDNA.fasta"
beddir="/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Motifs/homer/muscle/"
motifs="/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Archetypes/motif-archetypes-PPM-PCC-0.8-IC0.5-8bp-pwms.homer"
bedfg=${beddir}"/"peaks_${1}.bed
bedbg=${beddir}"/"peaks_${1}_bg.bed
outdir=${beddir}"/"peaks_${1}
echo ""
echo "Starting HOMER analysis for" $1
echo "using the intervals in" $bedfg
echo "Output will be saved to" $outdir
findMotifsGenome.pl $bedfg $genome $outdir -size 250 -len 6,8,10,12 -bg $bedbg #-nomotif #-mknown $motifs 
echo ""
echo "Finished de novo analysis for" $1
echo ""