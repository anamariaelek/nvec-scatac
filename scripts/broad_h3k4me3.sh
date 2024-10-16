# directories
scatac_dir=$( dirname $(pwd) )
chip_dir="/home/anamaria/cluster/aelek/proj/ChIPseq_nvec_v2"
plot_dir=${scatac_dir}"/plots/Peaks"
resu_dir=${scatac_dir}"/results/Peaks"

# inputs
genes=${scatac_dir}"/genome/Nvec_v4_merged_annotation_sort.bed"
bw_h3k4me3=${chip_dir}"/H3K4me3/alignments/H3K4me3.bw"
bw_h3k4me3_gastrula=${chip_dir}"/H3K4me3_gastrula/H3K4me3_gastrula_bwa.bw"

# # # # # # # # # # # #
#  CLUSTER AND PLOT   #
# # # # # # # # # # # #

# outputs
matrix=${resu_dir}"/deeptools.matrix.H3K4me3.gz"
heatmap=${plot_dir}"/deeptools.heatmap.H3K4me3.pdf"
k=7

# adult
computeMatrix scale-regions \
  -R ${genes} \
  -S ${bw_h3k4me3} \
  -b 0 -a 0 \
  --regionBodyLength 5000 \
  --skipZeros \
  -o ${matrix%%.gz}.adult.gz \
  --outFileNameMatrix ${matrix%%.gz}.adult.tab \
  --outFileSortedRegions ${matrix%%.gz}.adult.bed \
  -p 32

plotHeatmap -m ${matrix%%.gz}.adult.gz \
  -out ${heatmap%%.pdf}.k${k}.adult.pdf \
  --outFileSortedRegions ${heatmap%%.pdf}.k${k}.adult.bed \
  --colorMap Blues \
  --zMin 0 --zMax 30 \
  --kmeans ${k} \
  --legendLocation lower-center \
  --dpi 100 &

# gastrula
computeMatrix scale-regions \
  -R ${genes} \
  -S ${bw_h3k4me3_gastrula} \
  -b 0 -a 0 \
  --regionBodyLength 5000 \
  --skipZeros \
  -o ${matrix%%.gz}.gastrula.gz \
  --outFileNameMatrix ${matrix%%.gz}.gastrula.tab \
  --outFileSortedRegions ${matrix%%.gz}.gastrula.bed \
  -p 32

plotHeatmap -m ${matrix%%.gz}.gastrula.gz \
  -out ${heatmap%%.pdf}.k${k}.gastrula.pdf \
  --outFileSortedRegions ${heatmap%%.pdf}.k${k}.gastrula.bed \
  --colorMap Blues \
  --zMin 0 --zMax 30 \
  --kmeans ${k} \
  --legendLocation lower-center \
  --dpi 100 &

# # # # # # # # # # # # # # # # # # # # # # # # 
#  PLOT REFERENCE POINT FOR SELECTED CLUSTERS #
# # # # # # # # # # # # # # # # # # # # # # # # 

# adult
clust_broad="cluster_2"
genes_broad=${plot_dir}"/Nvec_v4_merged_annotation_sort_broad_adult.bed"
grep ${clust_broad} ${heatmap%%.pdf}.k${k}.adult.bed > ${genes_broad}

clust_control="cluster_4"
genes_control=${plot_dir}"/Nvec_v4_merged_annotation_sort_control_adult.bed"
grep ${clust_control} ${heatmap%%.pdf}.k${k}.adult.bed > ${genes_control}

computeMatrix reference-point \
  -R ${genes_broad} ${genes_control} \
  -S ${bw_h3k4me3} \
  -b 500 -a 3000 \
  --referencePoint TSS \
  --skipZeros --missingDataAsZero \
  -o ${matrix%%.gz}.adult.broad.gz \
  --outFileNameMatrix ${matrix%%.gz}.adult.tab \
  --outFileSortedRegions ${matrix%%.gz}.adult.broad.bed \
  -p 32

plotHeatmap -m ${matrix%%.gz}.adult.broad.gz \
  -out ${heatmap%%.pdf}.adult.pdf \
  --outFileSortedRegions ${heatmap%%.pdf}.adult.broad.bed \
  --colorMap Blues \
  --sortUsing region_length \
  --linesAtTickMarks \
  --zMin 0 --zMax 30 \
  --legendLocation lower-center \
  --dpi 100 &
  

# gastrula
clust_broad="cluster_1"
genes_broad=${plot_dir}"/Nvec_v4_merged_annotation_sort_broad_gastrula.bed"
grep ${clust_broad} ${heatmap%%.pdf}.k${k}.gastrula.bed > ${genes_broad}

clust_control="cluster_4"
genes_control=${plot_dir}"/Nvec_v4_merged_annotation_sort_control_gastrula.bed"
grep ${clust_control} ${heatmap%%.pdf}.k${k}.gastrula.bed > ${genes_control}

computeMatrix reference-point \
  -R ${genes_broad} ${genes_control} \
  -S ${bw_h3k4me3} \
  -b 500 -a 3000 \
  --referencePoint TSS \
  --skipZeros --missingDataAsZero \
  -o ${matrix%%.gz}.gastrula.broad.gz \
  --outFileNameMatrix ${matrix%%.gz}.gastrula.tab \
  --outFileSortedRegions ${matrix%%.gz}.gastrula.broad.bed \
  -p 32

plotHeatmap -m ${matrix%%.gz}.gastrula.broad.gz \
  -out ${heatmap%%.pdf}.gastrula.pdf \
  --outFileSortedRegions ${heatmap%%.pdf}.gastrula.broad.bed \
  --colorMap Blues \
  --sortUsing region_length \
  --linesAtTickMarks \
  --zMin 0 --zMax 30 \
  --legendLocation lower-center \
  --dpi 100 &
