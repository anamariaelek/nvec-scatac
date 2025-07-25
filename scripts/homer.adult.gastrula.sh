findMotifs(){
	genome="/users/asebe/aelek/proj/scATAC_nvec_v2/genome/Nvec_vc1.1_gDNA.fasta"
	beddir="/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Motifs/homer/adult_gastrula"
	#beddir="/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Motifs/homer/adult_gastrula_vs_bg"
	#beddir="/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Motifs/homer/adult_gastrula_vs_bg_archs"
	motifs="/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/Archetypes/motif-archetypes-PPM-PCC-0.8-IC0.5-5bp-pwms.homer"
	bedfg=${beddir}"/"peaks_${1}.bed
	bedbg=${beddir}"/"peaks_${1}_bg.bed
	outdir=${beddir}"/"peaks_${1}
	echo ""
	echo "Starting HOMER analysis for" $1
	echo "using the intervals in" $bedfg
	echo "Output will be saved to" $outdir
	findMotifsGenome.pl $bedfg $genome $outdir -size 250 -len 6,8,10,12 
	#findMotifsGenome.pl $bedfg $genome $outdir -size 250 -len 6,8,10,12 -bg $bedbg 
	#findMotifsGenome.pl $bedfg $genome $outdir -size 250 -len 6,8,10,12 -bg $bedbg -mknown $motifs -nomotif
	echo ""
	echo "Finished analysis for module" $1
}

for bed in gastrula_genes_gastrula gastrula_genes_adult gastrula_genes_both adult_genes_gastrula adult_genes_adult adult_genes_both both_genes_gastrula both_genes_adult both_genes_both 
do
	findMotifs "$bed" &
done
wait
echo "Done."

