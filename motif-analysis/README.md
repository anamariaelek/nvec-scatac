# Motif analysis pipeline

This is our internal pipeline for motif analysis.

## Environment

```bash
# create environment from scratch:
conda create -n mots
conda activate mots
conda install -c conda-forge r-base=4.1 r-doparallel r-plyr r-igraph r-circlize bioconductor-universalmotif bioconductor-pwmenrich bioconductor-genomicranges  r-wgcna bioconductor-rtracklayer  bioconductor-complexheatmap r-dbscan
# save
# conda env export --name mots > environment.yaml

# create from yaml file:
conda env create -n mots --file environment.yaml
conda activate mots
```

## How-to

1) Load libraries

```R
# load library
source("mta_downstream_functions.R")
# require("universalmotif")
# require("PWMEnrich")
# require("GenomicRanges")
# require("IRanges")
# require("rtracklayer")
# require("plyr")
# require("doParallel")
# require("WGCNA")
# require("igraph")
# require("circlize")
# require("ComplexHeatmap")
# require("data.table")
```

2) Find marker regulatory regions associated with specific cell types:

```R
# We'll take the regualtory regions from peaks and summits from an ATAC experiment
pka_fn = "test_data/macs_peaks.narrowPeak"
sum_fn = "test_data/macs_peaks.summits.bed"
# genome data
gtf_fn = "test_data/amphi.gtf"
gen_fn = "test_data/amphi.fasta"
gix_fn = "test_data/amphi.fasta.fai"
# promoters (this is optional, but it can improve the assignment of peaks to genes)
pro_fn = "test_data/promoters.bed"

# this function will take a list of ATAC peaks and match them to specific genes, based 
# on their relative positions in the genome, with the option to use promoter coordinates 
# to refine the assignment (e.g. H3K4me3 peaks).
hits_pka_gen = mta_match_peaks_to_genes(
  gff_object = gtf_fn, 
  peak_object = pka_fn, 
  index_object = gix_fn, 
  feature_to_match = "transcript",
  max_tss_dist = 20000,
  promoter_object = pro_fn)

# load a table with fold changes of your genes in different cell types
ct_footprint = read.table("test_data/footprints_per_ct.tsv", sep = "\t", header = TRUE)

# this function will take the peak-gene associations and find marker-specific peaks,
# which will be stored in BED files (and also in the `marker_peaks` object)
marker_peaks = mta_marker_peaks_per_ct(
  fp = ct_footprint,
  gene_peak_table = hits_pka_gen,
  peaks_object = sum_fn,
  index_object = gix_fn,
  fixed_peak_size = 250,
  threshold = 1.2,
  n_top_genes = 500,
  cell_types = c("endoderm","neural","mesodermal_somites","mesodermal_other"),
  save_bed = TRUE,
  bed_prefix = "test_data/marker_peaks")

# check it out
head(marker_peaks[["endoderm"]][["fg"]])
```

3) Create a HOMER library for your species using the foreground/background peaks identified above.

```bash
# just an example command for one cell type
mkdir -p test_data/homer_motifs
findMotifsGenome.pl \
  test_data/marker_peaks.endoderm.fg.bed \
  test_data/amphi.fasta \
  test_data/homer_motifs/motifs_endoderm \
  -size 250 -len 8,10,12,14 -mis 4 \
  -bg test_data/marker_peaks.endoderm.bg.bed \
  -mknown <path to known homer motif library> \
  -p 2 -noweight
```

4) Curate a HOMER motif library:

```R
#### Load motifs ####
# load a file with HOMER motifs
mot_fn = "test_data/motif_library.homer.txt"
mot = mta_read_homer_mod(mot_fn, pval_thr = 1e-3, ignore_SeqBias = TRUE, sanitise_names = TRUE, prefix = "motif")
# you can concatenate `mot` objects to load various libraries in R, e.g.:
# mot = c(mot1, mot2)

#### Motif filtering ####
# remove low-quality motifs by IC block filtering
keep_ixs = mta_filter_by_ic_block(motifs = mot, ic_thr = 0.5, len_uniblock = 4, len_multiblock = 3)
mot_f = mot [ keep_ixs ]

#### Motif merging ####
# merge motifs by similarity
mot_m = mta_merge_motifs_by_similarity(
  motifs = mot_f,
  matrix_class = "ICM", 
  threshold = 0.9, 
  dist_method = "WPCC", 
  clus_method = "hclust",
  min_overlap = 6,
  do_heatmap = TRUE)

# take a look at the merging heatmap
pdf("test_data/merged_motifs.pdf", height = 5, width = 5)
print(mot_m$matrix_hm)
dev.off()

# You can also merge motifs by similarity of their genome-wide motif energy tracks, using
# the function `mta_merge_motifs_by_gw_score_correlation` (much slower though).

# get the best representative motif for each cluster: in this case, we'll keep the motif with
# the highest IC score among those with at >=10 foreground sites (`min_nsites`). If no such motif
# is available, just keep the one with the highest IC score.
mot_m_b = mta_merge_motifs_find_best(motifs = mot_f, clusters_df = mot_m$clusters, criterion = "max_ic_min_fg", min_nsites = 10)

# get a curated motif library in `universalmotif` format:
mot_m = mot_m_b$motifs
```

5) Scan your genome with a curated motif library:

```R
#### Define regions of interest ####
# define regions of 250 around each peak summit that will be scanned for motif presence
# load genome index
gix = read.table(gix_fn, stringsAsFactors = FALSE)
gix = gix[,1:2]
colnames(gix) = c("chr","length")
# load peak summits
pka = mta_bed_to_granges(sum_fn)
seqlevels(pka)  = gix$chr
seqlengths(pka) = gix$length
# expand summits to width = 250 bp 
pka = pka + 125
pka = GenomicRanges::trim(pka)
pka = pka [ which(width(pka) == 251) ]

#### Define motif-specific reporting thresholds & perform genome-wide scan ####
# This function will:
# 1) obtain genome-wide quantiles of motif energy to use as empirical detection thresholds
# 2) scan the genome with these quantiles
mot_m_l = mta_convert_umot_to_monalisa(mot_m)
gw_scan_out = mta_gw_motif_score_monalisa(
  motifs = head(mot_m_l,1),  # list of motifs to scan in monalisa format
  genome_object = gen_fn,
  index_object = gix_fn,
  # regulatory regions to scan:
  given_gr = pka,          
  # width of genome bins used for empirical motif threshold scoring:
  bin_width = 250,         
  # list of score quantiles to keep for each motif in the list
  score_quantiles = c(0.00, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98, 0.99, 1.00),
  # This is the REPORTING quantile that'll be used if you perform a genome-wide scan. 
  # If you know which threshold you want to use, you can choose this value as a reporting 
  # threshold. If you don't, you can choose something lower and filter-out poor alignments
  # later on:
  score_quantile_thr = 0.9,
  # genome fraction that is scanned identify quantile-based motif-specific thresholds:
  subsample_fraction = 0.05,
  do_gw_scan = TRUE,
  nthreads = 4)
  
# Alternatively, you can use this function to obtain a similar output with p-values,
# obtained from comparing the alignment scores in your regions of interest to a background 
# distrubution in GC-matched random genomic bins. Output is a GRanges object with best hits.
mot_m_l = mta_convert_umot_to_monalisa(mot_m)
gw_scan_out_gr = mta_gw_motif_score_monalisa_fdr(
  motifs = mot_m_l,  # list of motifs to scan in monalisa format
  genome_object = gen_fn,
  index_object = gix_fn,
  # regulatory regions to scan:
  fg_r = pka,
  # background regions (set to NULL to use the rest of the genome):
  bg_r = NULL,
  bg_not_in_fg = TRUE,
  # subsample the background genome by this fraction (to speed up things)
  # (don't set it too low or you won't have enough background bins to get meaningful pvalues)
  subsample_fraction = 0.5,
  # width of genome bins used for background scoring (should match that of your foreground regions):
  bin_width = 250,
  # use this number of GC content intervals to match your foreground regions to the appropriate background
  # each bin will contain the same number of foreground regions, approximately.
  gcc_intervals = 20,
  # number of threads
  nthreads = 4)
  
# Alternatively, you can just scan the motifs directly using the 80% of the IC score as the 
# threshold (a common approach which we don't love, but it's MUCH faster because you can skip 
# the scoring step).
gw_scan_out_gr = mta_gw_scan_to_granges(
  mot_pwm_list = mot_m,
  genome_object = gen_fn,
  given_gr = pka,
  fractional_threshold = 0.8,
  nthreads_mot = 20,
  nthreads_chr = 1)

# output: a GRanges object with all motif presence coordinates above `score_quantile_thr` threshold
# (same as `gw_scan_out$gw_scan`).
```

4) Find cell type-specific motif enrichments

```R
# this should be run for each cell type separately
enri_mot = data.frame()
for (cell_type in c("endoderm","neural","mesodermal_somites","mesodermal_other")) {
  enri_mot_i = mta_motif_enrichment_test(
    sites_object = gw_scan_out$gw_scan,
    # we can reuse the foreground and background BED files obtained in the first step
    fg_object = sprintf("test_data/marker_peaks.%s.fg.bed", cell_type),
    bg_object = sprintf("test_data/marker_peaks.%s.bg.bed", cell_type),
    label = cell_type,
    thresholds_vector = gw_scan_out$score_quantiles[,"q0.98"],
    nthreads = 2, 
    pval_adjust = "fdr")
  # concatenate data for each cell type
  enri_mot = rbind(enri_mot, enri_mot_i)
}

# output: a table with motifs and their relative enrichments in each cell type
head(enri_mot)
# table of fold changes and table of p-values
enri_mot_m_fc = reshape2::dcast(enri_mot, motif ~ label, value.var="fc", drop = TRUE)
enri_mot_m_pv = reshape2::dcast(enri_mot, motif ~ label, value.var="padj", drop = TRUE)
# schit, check det hääär
```

Profit!

```R
 .              +   .                .   . .     .  .
                   .                    .       .     *
  .       *                        . . . .  .   .  + .
            "You Are Here"            .   .  +  . . .
.                 |             .  .   .    .    . .
                  |           .     .     . +.    +  .
                 \|/            .       .   . .
        . .       V          .    * . . .  .  +   .
           +      .           .   .      +
                            .       . +  .+. .
  .                      .     . + .  . .     .      .
           .      .    .     . .   . . .        ! /
      *             .    . .  +    .  .       - O -
          .     .    .  +   . .  *  .       . / |
               . + .  .  .  .. +  .
.      .  .  .  *   .  *  . +..  .            *
 .      .   . .   .   .   . .  +   .    .            +
```
