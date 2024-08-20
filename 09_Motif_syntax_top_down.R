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
ann_dir <- "annotation"
pks_dir <- "Results/Peaks"
mta_dir <- "Results/Motifs"
arc_dir <- "Results/Archetypes"
syn_dir <- "Results/Syntax"
dir.create(syn_dir, showWarnings = FALSE)
fig_dir <- "Plots/Syntax"
dir.create(fig_dir, showWarnings = FALSE)

# chromosome is the first argument passed to this script
chr <- as.character(commandArgs(trailingOnly = TRUE)[1])

# overlap fraction is the second argument passed to this script
reciprocal_overlap <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# load motif scores
message("Loading motif scores")
arc_id <- "PPM-PCC-0.8-IC0.5-5bp"
q <- 0.95
mta_dt <- rbindlist(lapply(c(
  # archetypes
  file.path(arc_dir, sprintf("motif-scores-archetypes-%s-mona-q%s.tsv.gz", arc_id, q)),
  # experimental motifs
  file.path(mta_dir, sprintf("motif-scores-mona-q%s.tsv.gz", q))
), fread))

# subset enriched motifs
mta_en <- fread(file.path(mta_dir, "motif-enrichment-mona-q0.98-FC-1-padj-0.001.tsv"))
mta_dt <- mta_dt[motif %in% mta_en$archetype_name]

# calculate relative motif scores (independent of the motif length)
mta_dt[, relative_motif_score := motif_score / max(motif_score)]
mta_dt <- unique(mta_dt[seqnames %in% chr])

# make ranges
mta_gr <- makeGRangesFromDataFrame(mta_dt, keep.extra.columns = TRUE)
rm(mta_dt)

# reduce motif hits
red_res <- mta_reduce_motif_hits(
  mta_gr, 
  seqnames = chr, 
  reciprocal_overlap = reciprocal_overlap,
  order_col = "relative_motif_score"
)
message("Done reducing hits for chromosome ", chr)

# save results
out_fn <- file.path(syn_dir, sprintf("motif-hits-reduced-%s.rds", chr))
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
