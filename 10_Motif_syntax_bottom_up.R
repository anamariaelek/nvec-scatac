# log time
message("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ")
message(
  "#                        ",
  Sys.time(),
  "                        #"
)
message("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ")

# load packages
library(GenomicRanges)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggrepel)
source("motif-analysis/mta_downstream_functions.R")

# set ggplot2 theme
theme_py <- theme_light() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA),
  text = element_text(size=20),
  strip.placement = "outside", 
  strip.text = element_text(size=20, color="black"),
  strip.background = element_rect(fill="white")
)
theme_set(theme_py)

# directories
adult_dir <- "ArchRProj_Nvec_TSS4_frag200"
gastr_dir <- "ArchRProj_Nvec_gastrula"
ann_dir <- "annotation"
pks_dir <- "results/Peaks"
mta_dir <- "results/Motifs"
arc_dir <- "results/Archetypes"
grn_dir <- "results/GRN"
syn_dir <- "results/Syntax"
dir.create(syn_dir, showWarnings = FALSE)
fig_dir <- "plots/Syntax"
dir.create(fig_dir, showWarnings = FALSE)

# cell types
ct_cols <- c(
  "cnidocyte"                  = "#ff42ff",
  "cnidocyte_gastrula"         = "#f7abf7",
  "ecto_pharynx"               = "#5bc0e8",
  "ectoderm"                   = "#51a0be",
  "ecto_aboral"                = "#045170",
  "EMS"                        = "#bdf5bd",
  "EMS_ecto_boundary"          = "#93dbce",
  "gastro_circular_muscle_1"   = "#85c90e",
  "gastro_circular_muscle_2"   = "#73b009",
  "gastro_parietal_muscle"     = "#8ceb10",
  "gastro_IRF1_2"              = "#c1eb05",
  "gastro_somatic_gonad"       = "#bde314",
  "muscle_tentacle_retractor"  = "#ffd700",
  "muscle_mesentery_retractor" = "#f0e229",
  "digestive_filaments_1"      = "#e33d3d",
  "digestive_filaments_2"      = "#d10606",
  "digestive_filaments_3"      = "#ad0303",
  "epidermis_1"                = "#04ccd4",
  "epidermis_2"                = "#16bacc",
  "precursors_PGC"             = "#bebebe",
  "precursors_endoNPC"         = "#8a8686",
  "precursors_NPC"             = "#636363",
  "NPC_1"                      = "#808d91",
  "NPC_2"                      = "#758d92",
  "neuron_GATA_Islet_1"        = "#0c82f7",
  "neuron_GATA_Islet_2"        = "#1175f0",
  "neuron_Pou4_FoxL2_1"        = "#101cde",
  "neuron_Pou4_FoxL2_2"        = "#0b16bf",
  "neuron_Pou4_FoxL2_3"        = "#2e39dd",
  "neuronal_gastrula"          = "#063cb9",
  "gland"                      = "#ff6f08",
  "gland_mucin"                = "#ff8f12"
)
cell_types <- names(ct_cols)
adult_cell_types <- c(
  "cnidocyte",
  "gastro_circular_muscle_1", 
  "gastro_circular_muscle_2",
  "gastro_parietal_muscle",
  "gastro_IRF1_2",
  "gastro_somatic_gonad",
  "muscle_mesentery_retractor",
  "muscle_tentacle_retractor",
  "digestive_filaments_1",
  "digestive_filaments_2",
  "digestive_filaments_3",
  "epidermis_1",
  "epidermis_2",
  "precursors_PGC",
  "precursors_endoNPC",
  "precursors_NPC",
  "neuron_GATA_Islet_1",
  "neuron_GATA_Islet_2",
  "neuron_Pou4_FoxL2_1",
  "neuron_Pou4_FoxL2_2",
  "neuron_Pou4_FoxL2_3",
  "gland"
)
gastr_cell_types <- c(setdiff(cell_types, adult_cell_types))

# parse arguments

# cell type is the first argument passed to this script
ct <- as.character(commandArgs(trailingOnly = TRUE)[1])

# overlap fraction is the second argument passed to this script
reciprocal_overlap <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# quantile for motif scores is the third argument passed to this script
q <- as.numeric(commandArgs(trailingOnly = TRUE)[3])

# which motifs to subset is the fourth argument passed to this script
# "enrich" OR "enrich-assign" OR "active"
subset_motifs <- as.character(commandArgs(trailingOnly = TRUE)[4])

# subdirectory where to save the results
res_dir <- sprintf("q%s-hits-%s-motifs-ovl-%s-all-pks", q, subset_motifs, reciprocal_overlap)
dir.create(file.path(syn_dir, res_dir), showWarnings = FALSE)

# starting analysis
message("")

# load motif scores
message(sprintf("%s | Loading motif scores above %s quantile", Sys.time(), q))
arc_id <- "PPM-PCC-0.8-IC0.5-5bp"
mta_dt <- rbindlist(lapply(c(
  # archetypes
  file.path(arc_dir, sprintf("motif-scores-archetypes-%s-mona-q%s.tsv.gz", arc_id, q)),
  # experimental motifs
  file.path(mta_dir, sprintf("motif-scores-mona-q%s.tsv.gz", q))
), fread))

# # in silico ChIP score filtered peaks
# ics_fn <- file.path(syn_dir, "motif-hits-ics-filtered.tsv.gz")
# # if filtered peaks file doesn't exist, create it
# if (!file.exists(ics_fn)) {
#   message("Filtering motif hits by in silico ChIP score")
#   lvl <- "metacell"
#   id <- "genes_exp_FC2_acc_FC4_spearman"
#   thr <- 0.1
#   ics_dt <- fread(file.path(
#     grn_dir, lvl, sprintf(
#       "insilico-chip-binding-score-%s.tsv.gz", id
#     )
#   ))
#   ics_dt <- unique(ics_dt[in_silico_chip_score > thr, .(seqnames, start, end, peak, motif)])
#   fwrite(ics_dt, file.path(syn_dir, "motif-hits-ics-filtered.tsv.gz"), sep = "\t")
# } else {
#   message("Loading filtered motif hits: ", ics_fn)
#   ics_dt <- fread(ics_fn, sep = "\t")
# }
# 
# # subset motifs hits in peaks
# mta_dt <- merge.data.table(mta_dt, ics_dt[, .(peak, motif)], by = c("peak", "motif"))
# mta_dt <- unique(mta_dt)

# enriched motifs
mta_en <- fread(file.path(mta_dir, sprintf("motif-enrichment-mona-q%s-FC-1-padj-0.001.tsv", q)))

# assigned motifs
mta_as <- fread(file.path(arc_dir, "motif-assignment-archetypes-PPM-PCC-0.8-IC0.5-5bp.tsv"))

# active TFs motifs
grn_tfs <- fread(file.path(
  grn_dir,
  "networks",
  "grn_tfs_info_expression_fc_0.4_accessibility_access_0.4_chromvar_4.tsv"
))
grn_tfs <- unique(grn_tfs[cell_type == ct, .(motif, gene, zscore, cell_type)])

# subset motifs
if (subset_motifs == "enrich") {
  # enriched motifs
  mta_dt <- mta_dt[motif %in% mta_en$archetype_name]  
} else if (subset_motifs == "enrich|assign") {
  # enriched and assigned motifs
  mta_dt <- mta_dt[motif %in% c(mta_en$archetype_name, mta_as$archetype_name)]  
} else if (subset_motifs == "active") {
  # motifs assigned to active TFs
  mta_dt <- mta_dt[motif %in% grn_tfs$motif]
}
message(sprintf("%s | Subset %s motifs", Sys.time(), subset_motifs))

# enrichment values
mta_en <- rbindlist(lapply(c(
  file.path(mta_dir, sprintf("motif-enrichment-cell-type-mona-q%s.tsv", q)),
  file.path(arc_dir, sprintf("motif-enrichment-cell-type-archetypes-PPM-PCC-0.8-IC0.5-5bp-mona-q%s.tsv", q))
), fread))
mta_en <- unique(mta_en[,.(cell_type, motif, fc, pval, padj)])

# all peaks
peaks <- fread(file.path(pks_dir, "Peaks_cell_type_mapped.bed"))
setnames(peaks, c("seqnames", "start", "end", "peak", "score", "strand"))
peaks <- unique(peaks)
peaks_gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)

# cell type peaks
message(sprintf("%s | Reducingp peaks for %s", Sys.time(), ct))
stg_dir <- ifelse(ct %in% adult_cell_types, adult_dir, gastr_dir)
#pks_fn <- file.path(
#  stg_dir, "ArchRProj", "PeakDifferential", "cell_type_filtered",
#  sprintf("Peaks-%s-vs-others.tsv", ct)
#)
pks_fn <- file.path(
  stg_dir, "ArchRProj", "Peaks", "cell_type_filtered",
  sprintf("Peaks-%s.tsv", ct)
)
pks_dt <- fread(pks_fn, select = 1:3)
setnames(pks_dt, c("seqnames", "start", "end"))
pks_gr <- makeGRangesFromDataFrame(pks_dt)

# overlap peaks (to map per-stage peak ids to final peak set ids)
pks_ovl <- findOverlaps(query = peaks_gr, subject = pks_gr)
pks_dat <- peaks_gr[queryHits(pks_ovl)]

# peaks in cell type
mta_ct_dt <- mta_dt[peak %in% pks_dat$peak]

# motif enrichment in cell type
mta_ct_en <- mta_en[cell_type == ct]#[padj < 0.05]

# combine
mta_ct <- merge.data.table(mta_ct_dt, mta_ct_en, by = "motif")

# make ranges
mta_gr <- makeGRangesFromDataFrame(mta_ct, keep.extra.columns = TRUE)

# reduce motif hits
red_res_ls <- lapply(unique(seqnames(mta_gr)), function(x) {
  message(sprintf("%s | Reducing hits for %s", Sys.time(), x))
  mta_reduce_motif_hits(
    mta_gr, seqnames = x,
    reciprocal_overlap = reciprocal_overlap,
    order_col = "fc",
    order_decrease = TRUE
  )
})
red_res <- sapply(
  c("select_hits", "reduce_hits", "inputs_hits"), function(x) {
    ls <- lapply(red_res_ls, function(y) y[[x]])
    do.call("c", ls)
}, USE.NAMES = TRUE, simplify = FALSE)
message(sprintf("%s | Done reducing hits for %s", Sys.time(), ct))

# save results
out_fn <- file.path(syn_dir, res_dir, sprintf("motif-hits-reduced-%s.rds", ct))
saveRDS(red_res, out_fn)
message("Saved to: ", out_fn)

# log time
message("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ")
message(
  "#                        ",
  Sys.time(),
  "                        #"
)
message("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ")
