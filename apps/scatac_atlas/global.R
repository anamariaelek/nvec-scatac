# packages
options(repos = BiocManager::repositories())
library(ape)
library(BiocManager)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(devtools)
library(dplyr)
library(DT)
library(eulerr)
library(ggiraph)
require(ggtree)
library(ggplot2) 
#library(ggpubr)
library(ggrepel)
library(gplots)
library(kableExtra)
library(Matrix)
library(patchwork)
library(RColorBrewer)
library(R.utils)
library(reshape2)
library(scales)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(stringr)
require(treeio)
require(tidytree)
library(viridis)
library(XML)
library(zoo)

# Load functions --------------------------------------------------------------
source("functions.R")
source("modules.R")


# Input data ------------------------------------------------------------------

# read config file
conf <- yaml::yaml.load_file("config.yaml", eval.expr=TRUE)
INPUT_DIR <- conf[['default']]$data_dir
dir.create(INPUT_DIR, showWarnings=FALSE)
overwrite <- FALSE

gene_annotation <-  conf[['default']]$gene_annot_file
tfs_annotation <-  conf[['default']]$tf_annot_file

GENE_ANNT <- fread_gene_annotation(
  file = file.path(INPUT_DIR, gene_annotation),
  cols=1:4,
  colnames=c("gene_id","gene name","PFAM domain","orthogroup")
)
for (i in 2:ncol(GENE_ANNT))
  GENE_ANNT[is.na(GENE_ANNT[[i]]),colnames(GENE_ANNT)[i]:="-"]
GENE_ANNT[search_id=='',search_id:=gene_id]
TF_ANNT <- fread(file.path(INPUT_DIR,tfs_annotation), header=TRUE)
setnames(TF_ANNT,c("gene_id","gene name","PFAM domain","orthogroup"))

# graphical params
heatmap_colors <- c("white","gray99","orange","orangered2","#520c52")
proj_colors <- c("gray95","lightyellow","khaki1","orange","orangered2","#520c52")

# experiments
samples <- c(
  # whole organism
  "25m_PFA","50m_PFA",
  "3_FixSor", "9_F6k", "10_F6k", "11_F5k", "12_F", "13_F", 
  "16_FC", 
  "20_Fn", "21_Fn", "22_Fn", "23_Fn", 
  "24_Fn", 
  "25_Fn", "26_Fn", "27_Fn",
  # elav
  "18_Elav", "19_Elav",
  "28_Elav", "29_Elav", "30_Elav",
  # multiome
  "Multiome_07563AAD",
  # gastrula
  "G26_fx", "G26_fs"
)

# assign colors to samples
orig_cols <- c(
  # whole organism
  hsv(h=c(0.08,0.1), s=1, v=1),
  hsv(h=c(0.14,0.16,0.18,seq(0.2,0.45,length.out=3)), s=1, v=1),
  hsv(h=0.5, s=1, v=1),
  hsv(h=seq(0.55,0.7,length.out=4), s=1, v=1),
  hsv(h=0.75, s=1, v=1),
  hsv(h=seq(0.8,0.85,length.out = 3), s=1, v=1),
  # elav
  hsv(h=c(0.94,0.96), s=1, v=1),
  hsv(h=c(0.98, 0.99, 1), s=1, v=1),
  # multiome
  "#4D4D4D",
  # gastrula
  "#ccff00", "#cc00ff"
)
orig_cols <- structure(orig_cols, names = samples)

# colors
ct_cols <- c(
  cnidocyte                  = "#ff42ff",
  cnidocyte_gastrula         = "#f7abf7",
  ecto_pharynx               = "#5bc0e8",
  ectoderm                   = "#51a0be",
  ecto_aboral                = "#045170",
  EMS                        = "#bdf5bd",
  EMS_ecto_boundary          = "#93dbce",
  gastro_circular_muscle_1   = "#85c90e",
  gastro_circular_muscle_2   = "#73b009",
  gastro_parietal_muscle     = "#8ceb10",
  gastro_IRF1_2              = "#c1eb05",
  gastro_somatic_gonad       = "#bde314",
  muscle_tentacle_retractor  = "#ffd700",
  muscle_mesentery_retractor = "#f0e229",
  digestive_filaments_1      = "#e33d3d",
  digestive_filaments_2      = "#d10606",
  digestive_filaments_3      = "#ad0303",
  epidermis_1                = "#04ccd4",
  epidermis_2                = "#16bacc",
  precursors_PGC             = "#bebebe",
  precursors_endoNPC         = "#8a8686",
  precursors_NPC             = "#636363",
  NPC_1                      = "#808d91",
  NPC_2                      = "#758d92",
  neuron_GATA_Islet_1        = "#0c82f7",
  neuron_GATA_Islet_2        = "#1175f0",
  neuron_Pou4_FoxL2_1        = "#101cde",
  neuron_Pou4_FoxL2_2        = "#0b16bf",
  neuron_Pou4_FoxL2_3        = "#2e39dd",
  neuronal                   = "#063cb9",
  gland                      = "#ff6f08",
  gland_mucin                = "#ff8f12"
)
cell_types <- names(ct_cols)
cell_type_rename <- c(
  "cnidocyte_precursors" = "NPC_2",
  "^muscle_1" = "muscle_tentacle_retractor",
  "^muscle_2" = "muscle_mesentery_retractor",
  "adult_muscle_1" = "adult_muscle_tentacle_retractor",
  "adult_muscle_2" = "adult_muscle_mesentery_retractor",
  "gastro_unknown_1" = "gastro_IRF1_2",
  "gastro_unknown_2" = "gastro_somatic_gonad",
  "gastrula_gland" = "gastrula_gland_mucin",
  "ectoderm_embryonic_oral" = "ecto_pharynx",
  "ectoderm_embryonic$" = "ectoderm",
  "ectoderm_embryonic_aboral" = "ecto_aboral",
  "mesendoderm_embryonic" = "EMS",
  "mesendoderm_ectoderm" = "EMS_ecto_boundary",
  "neuron_GATA_Islet_3" = "neuron_Pou4_FoxL2_3",
  "neuronal_gastrula" = "neuronal",
  "NPC$" = "NPC_1",
  "precursors_1" = "precursors_PGC",
  "precursors_2" = "precursors_endoNPC",
  "precursors_3" = "precursors_NPC"
)