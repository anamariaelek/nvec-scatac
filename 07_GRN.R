
# global settings
options(stringsAsFactors = FALSE)

# load packages
library(RColorBrewer)
library(BSgenome.jaNemVect1.1.DToL.Assembly)
library(GenomicRanges)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(openxlsx)

# set ggplot2 theme
theme_blank <- theme_minimal() + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  text = element_text(size=20)
)
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


adult_dir <- "ArchRProj_Nvec_TSS4_frag200"
gastr_dir <- "ArchRProj_Nvec_gastrula"
pks_dir <- "Results/Peaks"
mta_dir <- "Results/Motifs"
arc_dir <- "Results/Archetypes"
map_dir <- "Results/Metacells"
grn_dir <- "Results/GRN"
dir.create(grn_dir, showWarnings = FALSE)
dir.create(file.path(grn_dir, "metacell"), showWarnings = FALSE)
dir.create(file.path(grn_dir, "cell_type"), showWarnings = FALSE)
dir.create(file.path(grn_dir, "networks"), showWarnings = FALSE)
fig_dir <- "Plots/GRN"
dir.create(fig_dir, showWarnings = FALSE)
ann_dir <- "annotation"


# gene annotation
gnan <- fread(file.path(
  ann_dir, "Nematostella_DToL_FINAL.tsv"
))

# TF annotation
tfan <- fread(file.path(
  ann_dir, "Nematostella_DToL_TFs_FINAL.tsv"
))

# golden markers
gold <- fread(file.path(
  ann_dir, "golden-marks-231124.tsv"
), header = FALSE)
setnames(gold, c("common_name", "gene", "remark"))


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
bct_cols <- toupper(c(
  "cnidocyte"                 = "#FF42FF",
  "ecto"                      = "#51a0be",
  "EMS"                       = "#bdf5bd",
  "gastro_circular_muscle"    = "#73b009",
  "gastro_parietal_muscle"    = "#8ceb10",
  "gastro"                    = "#85c90e",  
  "muscle"                    = "#FFD700",
  "digestive_filaments"       = "#e33d3d",
  "precursors"                = "#BEBEBE",
  "NPC"                       = "#808d91",
  "epidermis"                 = "#04ccd4",
  "neuron_GATA_Islet"         = "#1175f0",
  "neuron_Pou4_FoxL2"         = "#101cde",
  "neuronal"                  = "#063cb9",
  "gland"                     = "#ff6f08"
))
bct_maps <- setDT(cbind.data.frame(
  cell_type = cell_types,
  broad_cell_type = str_extract(cell_types, paste(names(bct_cols), collapse = "|"))
))



# how peak-TF correlations were calculated
lvl <- "metacell"

# mapping between metacells
id <- "genes_exp_FC2_acc_FC4_spearman"

# in silico ChIP threshold
thrs <- 0.1

# expression quantification: fc or umifrac
exp <- "fc"

# accessibility quantification: access or score
acc <- "access"

# quantile thresholds for expression and accessibility
thr_q_exp <- 0.4
thr_q_acc <- 0.4

# chromvar threshold for TF motif activity
chr_thrs <- 4

# this is a global expression threshold for all cell types;
# if not NULL, it overrides expression quantification and quantile threshold!
thr_fc_exp <- NULL
if (!is.null(thr_fc_exp)) exp <- "fc"

# file name
if (is.null(thr_fc_exp)) {
  fn <- sprintf(
    "expression_%s_%s_accessibility_%s_%s_chromvar_%s", 
    exp, thr_q_exp, acc, thr_q_acc, chr_thrs
  )
} else {
  fn <- sprintf(
    "expression_%s_%s_accessibility_%s_%s_chromvar_%s", 
    exp, thr_fc_exp, acc, thr_q_acc, chr_thrs
  )
}



# global GRN data
grn_dt <- readRDS(file.path(
    grn_dir, lvl, sprintf("insilico-chip-grn-%s.rds", id)
))
grn_dt <- unique(grn_dt[in_silico_chip_score > thrs])

# global network targets
grn_glob <- unique(grn_dt[, .(gene, target_gene)])
setnames(grn_glob, "gene", "tf_gene")


# gene developmental stage data
devel_dt <- fread(file.path(
    grn_dir, "development_data_for_grn.tsv"
))
stopifnot(all(devel_dt$cell_type %in% cell_types))
devel_dt[, cell_type := factor(cell_type, levels = cell_types)]


# chromvar data 
exp_chromvar_dt <- fread(
    file.path(
        grn_dir, lvl,
        sprintf("gene_expression_%s_chromVAR_%s.tsv.gz", exp, id)
    )
)
exp_chromvar_dt <- exp_chromvar_dt[insilico_ChIP_threshold == thrs]
stopifnot(all(exp_chromvar_dt$cell_type %in% cell_types))
exp_chromvar_dt[, cell_type := factor(cell_type, levels = cell_types)]

# aggregate zscores per cell type
exp_chromvar_dt[, zscore := mean(zscore), .(gene, motif, cell_type)]
exp_chromvar_dt[, expression := mean(expression), .(gene, motif, cell_type)]
exp_chromvar_dt[, c("metacell", "expression", "insilico_ChIP_threshold") := NULL]
exp_chromvar_dt <- unique(exp_chromvar_dt)




# expression and accessibility data for base network
gns_exp_acc_list <- readRDS(
    file.path(
        grn_dir, "cell_type", "grn_tf_targets_expression_accessibility.rds"
    )
)
gns_exp_acc_dt <- rbindlist(gns_exp_acc_list)

# combine withh cromVAR scores
gns_exp_acc_act_dt <- merge.data.table(
    gns_exp_acc_dt, exp_chromvar_dt,
    by = intersect(colnames(exp_chromvar_dt), colnames(gns_exp_acc_dt)),
    all.x = TRUE, sort = FALSE
)
setcolorder(gns_exp_acc_act_dt, c(
  colnames(gns_exp_acc_dt)[1 : grep("target", colnames(gns_exp_acc_dt))[1]-1],
  "zscore"
))



# per cell type thresholds of expression and accessibility
data_thrs <- fread(file.path(
    grn_dir, "cell_type", sprintf(
        "qthreshold_%s_expression_%s_accessibility_%s.tsv",
        id, thr_q_exp, thr_q_acc
    )
))
data_thrs[cell_type %in% adult_cell_types, stage := "adult"]
data_thrs[cell_type %in% gastr_cell_types, stage := "gastrula"]
stopifnot(all(data_thrs$cell_type %in% cell_types))
data_thrs[, cell_type := factor(cell_type, levels = cell_types)]
if (exp == "umifrac") {
    data_thrs[, exp_thrs := thr_exp_umifrac]
} else if (exp == "fc") {
    data_thrs[, exp_thrs := thr_exp_fc]
}
if (!is.null(thr_fc_exp)) (
  data_thrs[, exp_thrs := thr_fc_exp]
)
if (acc == "score") {
    data_thrs[, acc_thrs := thr_gene_score]
} else if (acc == "access") {
    data_thrs[, acc_thrs := thr_accessibility]
}

# combine
grn_base_dt <- merge.data.table(
  gns_exp_acc_act_dt, data_thrs[, .(cell_type, stage, exp_thrs, acc_thrs)], 
  by = c("cell_type", "stage")
)




# determine which TFs are active in which cell type
# based on mean expression and z score
if (exp == "fc") {
  grn_filt_dt <- grn_base_dt[expression_fc > exp_thrs & zscore > chr_thrs]
} else if (exp == "umifrac") {
  grn_filt_dt <- grn_base_dt[expression_umifrac > exp_thrs & zscore > chr_thrs]
}



# filter target genes per cell type
# based on mean expression and accessibility/gene score
if (exp == "fc") {
  grn_filt_dt <- grn_filt_dt[target_expression_fc > exp_thrs]
} else if (exp == "umifrac") {
  grn_filt_dt <- grn_filt_dt[target_expression_umifrac > exp_thrs]
}
if (acc == "access") {
  grn_filt_dt <- grn_filt_dt[target_accessibility > acc_thrs]
} else if (exp == "score") {
  grn_filt_dt <- grn_filt_dt[target_gene_score > acc_thrs]
}


# auto regulating TFs per cell type
grn_auto_dt <- unique(grn_filt_dt[, .(gene, target_gene, cell_type)])
grn_auto_dt <- rbindlist(lapply(cell_types, function(ct) {
  auto_reg_gns <- unique(grn_auto_dt[cell_type == ct][gene == target_gene]$gene)
  data.table(gene = auto_reg_gns)[, cell_type := ct][, target_self := TRUE]
}))
grn_filt_dt <- merge.data.table(
  grn_filt_dt, grn_auto_dt, by = c("gene", "cell_type"), 
  all.x = TRUE, sort = FALSE
)
grn_filt_dt[is.na(target_self), target_self := FALSE]


# TF info table
tfs_info_dt <- unique(grn_filt_dt[, .(
  cell_type, stage, gene, gene_name, common_name, og, pfam, 
  expression_fc, expression_umifrac, gene_score, zscore, motif, id, target_self
)])
stopifnot(nrow(tfs_info_dt[,.N,.(cell_type,gene)][N>1])==0)

# save
fwrite(
  tfs_info_dt,
  file.path(
    grn_dir, "networks", sprintf("grn_tfs_info_%s.tsv", fn)
  ),
  sep = "\t"
)
require(openxlsx)
wb <- createWorkbook()
for (ct in cell_types) {
  sheetName <- substr(ct, 1, 30)
  grn_ct <- tfs_info_dt[cell_type == ct]
  addWorksheet(wb, sheetName = sheetName)
  writeData(wb, sheet = sheetName, grn_ct)
}
saveWorkbook(
  wb,
  file.path(grn_dir, "networks", sprintf("grn_tfs_info_%s.xlsx", fn)),
  overwrite = TRUE
)




# any TFs that have no target genes
tfs_info_dt[, stage_cell_type_gene := paste(stage, cell_type, gene)]
grn_filt_dt[, stage_cell_type_gene := paste(stage, cell_type, gene)]
grn_fill_dt <- rbindlist(list(
  grn_filt_dt,
  tfs_info_dt[! stage_cell_type_gene %in% grn_filt_dt$stage_cell_type_gene]
), use.names = TRUE, fill = TRUE)
grn_fill_dt[, exp_thrs := .SD[!is.na(exp_thrs)]$exp_thrs[1], .(cell_type, stage)]
grn_fill_dt[, acc_thrs := .SD[!is.na(acc_thrs)]$acc_thrs[1], .(cell_type, stage)]
grn_fill_dt[, cell_type := factor(cell_type, levels = cell_types)]
grn_fill_dt[, stage_cell_type_gene := NULL]
tfs_info_dt[, stage_cell_type_gene := NULL]
setorder(grn_fill_dt, cell_type, stage)

# target TF annotations
grn_fill_dt[, target_TF := target_gene %in% tfan$gene]
grn_fill_dt[, target_active_TF := target_gene %in% .SD$gene, .(cell_type, stage)]

# save
fwrite(
  grn_fill_dt,
  file.path(
    grn_dir, "networks", sprintf("grn_peaks_%s.tsv", fn)
  ),
  sep = "\t"
)
require(openxlsx)
wb <- createWorkbook()
for (ct in cell_types) {
  sheetName <- substr(ct, 1, 30)
  grn_ct <- grn_fill_dt[cell_type == ct]
  addWorksheet(wb, sheetName = sheetName)
  writeData(wb, sheet = sheetName, grn_ct)
}
saveWorkbook(
  wb,
  file.path(grn_dir, "networks", sprintf("grn_peaks_%s.xlsx", fn)),
  overwrite = TRUE
)


# aggregate peaks per gene
grn_gene_dt <- copy(grn_fill_dt)
grn_gene_dt[, c("target_peak", "seqnames", "start", "end") := NULL]
grn_gene_dt[, target_accessibility := mean(target_accessibility), .(gene, target_gene, cell_type, stage)]
grn_gene_dt[, peak_tf_correlation_score := mean(peak_tf_correlation_score), .(gene, target_gene, cell_type, stage)]
grn_gene_dt[, in_silico_chip_score := mean(in_silico_chip_score), .(gene, target_gene, cell_type, stage)]
grn_gene_dt <- unique(grn_gene_dt)

# save
fwrite(
  grn_gene_dt,
  file.path(
    grn_dir, "networks", sprintf("grn_genes_%s.tsv", fn)
  ),
  sep = "\t"
)
require(openxlsx)
wb <- createWorkbook()
for (ct in cell_types) {
  sheetName <- substr(ct, 1, 30)
  grn_ct <- grn_gene_dt[cell_type == ct]
  addWorksheet(wb, sheetName = sheetName)
  writeData(wb, sheet = sheetName, grn_ct)
}
saveWorkbook(
  wb,
  file.path(grn_dir, "networks", sprintf("grn_genes_%s.xlsx", fn)),
  overwrite = TRUE
)


grn_tfs_dt <- copy(grn_gene_dt)
grn_tfs_dt <- grn_tfs_dt[target_TF == TRUE]

# any TFs that have no target genes
tfs_info_dt[, stage_cell_type_gene := paste(stage, cell_type, gene)]
grn_tfs_dt[, stage_cell_type_gene := paste(stage, cell_type, gene)]
grn_tfs_dt <- rbindlist(list(
  grn_tfs_dt,
  tfs_info_dt[! stage_cell_type_gene %in% grn_tfs_dt$stage_cell_type_gene]
), use.names = TRUE, fill = TRUE)
grn_tfs_dt[, exp_thrs := .SD[!is.na(exp_thrs)]$exp_thrs[1], .(cell_type, stage)]
grn_tfs_dt[, acc_thrs := .SD[!is.na(acc_thrs)]$acc_thrs[1], .(cell_type, stage)]
grn_tfs_dt[, cell_type := factor(cell_type, levels = cell_types)]
grn_tfs_dt[, stage_cell_type_gene := NULL]
tfs_info_dt[, stage_cell_type_gene := NULL]
setorder(grn_tfs_dt, cell_type, stage)

# save
fwrite(
  grn_tfs_dt,
  file.path(
    grn_dir, "networks", sprintf("grn_tfs_%s.tsv", fn)
  ),
  sep = "\t"
)
require(openxlsx)
wb <- createWorkbook()
for (ct in cell_types) {
  sheetName <- substr(ct, 1, 30)
  grn_ct <- grn_tfs_dt[cell_type == ct]
  addWorksheet(wb, sheetName = sheetName)
  writeData(wb, sheet = sheetName, grn_ct)
}
saveWorkbook(
  wb,
  file.path(grn_dir, "networks", sprintf("grn_tfs_%s.xlsx", fn)),
  overwrite = TRUE
)
