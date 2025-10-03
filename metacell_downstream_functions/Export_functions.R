require(metacell)
require(Seurat)
require(data.table)
require(stringr)
require(ggplot2)
require(anndata)


#' Export Metacell solution as Seurat object
#' 
#' @param mat_id character, id of single cell matrix object in scdb
#' @param mc_id character, id of metacell object in scdb
#' @param keep_single_cells logical, keep single cell resolution? Default `TRUE`, if `FALSE`, 
#'   use metacells in place of single cells in seurat object
#' @param mc_ann character, path to file containing metacell annotation, this should be a tsv file
#'   with the three columns containing metacell, cell type and color annotaions
#' @param gset_id id of markers gene set in scdb (used for building Metacell solution)
#'   if this is not specified, `nfeatures` will be selected by Seurat
#' @param nfeatures integer, number of variable features to select (using Seurat), 
#'   this doesn't have to be specified if `gset_id` is used (default: NULL)
#' @param dims dimensions of reduction to use as input, passed to `Seurat::FindNeighbors` 
#'   and `Seurat::RunUMAP`
#' @param resolution resolution parameter for `Seurat::FindClusters`
#' @param k number of neighbors for `Seurat::FindNeighbors`
#' @param seurat_rds character, path to file where to save the Seurat object (default: NULL)
#' @param out_dir character, where to save the object and plots (default: working directory)
#' @param verbose logical, print detailed log messages (default: TRUE)
#' @param ... Other arguments passed to Seurat methods
#' 
#' @return Seurat object of class `Assay`. Also saves 2D projection plots to `out_dir`.
#' 
mc_export_seurat <- function(
  mat_id, 
  mc_id, 
  mc_ann, 
  gset_id, 
  keep_single_cells=FALSE, 
  cluster=TRUE, 
  nfeatues=NULL, 
  dims=1:10, 
  resolution=0.5,  
  k=30, 
  seurat_rds=NULL, 
  out_dir=".", 
  verbose=TRUE) {
  
  # TBA check that scdb is initialized...
  # TBA create object w/o cell type annotations
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # seurat object
  if (verbose) message("Loading matrix data")
  mat <- scdb_mat(mat_id)
  if (verbose) message("Gene names must not contain '_', replacing them with '-'.")
  rownames(mat@mat) <- stringr::str_replace_all(rownames(mat@mat),"_","-")
  if (verbose) message("Loading metacell data")
  mc <- scdb_mc(mc_id)
  
  if (keep_single_cells==TRUE) {
    mc_cells <- intersect(colnames(mat@mat),names(mc@mc))
    if (verbose) message("Creating Seurat object with ",length(mc_cells)," cells.")
    seurat_obj <- CreateSeuratObject(counts=mat@mat[,mc_cells], project=mc_id)
  } else {
    if (verbose) message("Creating Seurat object with ",dim(mc@mc_fp)[2]," metacells.")
    mc_counts <- sca_mc_gene_counts(mc, mat)
    seurat_obj <- CreateSeuratObject(counts=mc_counts, project=mc_id)
  }
  
  # metadata
  seurat_obj@meta.data$orig.ident=mc_id
  if (!is.null(mc_ann)) {
    if (verbose) message("Loading metacell annotations")
    ann <- fread(mc_ann,header=TRUE)
    setnames(ann,colnames(ann)[1:3],c("metacell","cell_type","color"))
    ann[, metacell:=as.character(metacell)]
    cell_asgn <- data.table(cell=names(mc@mc),metacell=as.character(mc@mc))
    cell_ann <- merge.data.table(cell_asgn,ann,by="metacell",all.x=TRUE)
    cell_ann[is.na(cell_type),c("cell_type","color"):=list("unassigned","gray")]
    
    if (keep_single_cells==TRUE) {
      meta_dt <- data.table(seurat_obj@meta.data,keep.rownames="cell")
      meta_dt_ann <- merge.data.table(meta_dt,cell_ann,by="cell",all.x=TRUE)
      meta_dt_ann[is.na(cell_type),c("cell_type","color"):=list("unassigned","gray")]
    } else {
      meta_dt <- data.table(seurat_obj@meta.data,keep.rownames="metacell")
      meta_dt[, metacell:=as.character(metacell)]
      meta_dt_ann <- merge.data.table(meta_dt,ann,by="metacell",all.x=TRUE)
    }
    meta <- copy(meta_dt_ann); class(meta) <- "data.frame"; rownames(meta) <- meta_dt_ann[[1]]
    seurat_obj@meta.data <- meta
    Idents(seurat_obj) <- structure(meta$orig.ident, names=rownames(meta))
    Idents(seurat_obj) <- Idents(seurat_obj)[colnames(seurat_obj)]
    
  } else {
    stop("You have to specify metacell annotation file!")
  }
  
  
  # if not using metacells, need to do normalization and clustering
  if (keep_single_cells==TRUE) {
    # normalize
    if (verbose) message("Normalizing data.")
    seurat_obj <- NormalizeData(seurat_obj)

    # scaling
    if (verbose) message("Scaling data.")
    seurat_obj <- ScaleData(seurat_obj)
    
    # variable features
    if (!is.null(gset_id)) {
      marker_gset <- scdb_gset(gset_id)
      marker_genes <- names(marker_gset@gene_set)
      if (verbose) message("Gene names must not contain '_', replacing them with '-'.")
      marker_genes <- stringr::str_replace_all(marker_genes,"_","-")
      variable_genes <- intersect(marker_genes,rownames(mat@mat))
      if (verbose) message("Supplied ",length(marker_genes)," marker genes; using ",length(variable_genes)," variable genes.")
      VariableFeatures(seurat_obj) <- marker_genes
    } else if (!is.null(nfeatures)) {
      seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=3000)
      if (verbose) message("Selecting ",nfeatures," variable genes.")
    } else {
      stop("You have to specify either gset_id or nfeatures!")
    }
    
    # clustering
    if (cluster==TRUE) {
      if (verbose) message("Running PCA.")
      seurat_obj <- RunPCA(seurat_obj, features=VariableFeatures(object=seurat_obj), verbose=verbose)
      if (verbose) message("Finding clusters.")
      seurat_obj <- FindNeighbors(seurat_obj, dims=dims, verbose=verbose)
      seurat_obj <- FindClusters(seurat_obj, resolution=resolution, k.param=k, verbose=verbose)
      if (verbose) message("Running UMAP.")
      seurat_obj <- RunUMAP(seurat_obj, dims=dims, verbose=verbose)
      if (verbose) message("Plotting 2d projection.")
      gp_col_dt <- unique(meta_dt_ann[,.(cell_type,color)])
      if (verbose) print(gp_col_dt)
      gp_col <- gp_col_dt$color
      names(gp_col) <- gp_col_dt$cell_type
      plots <- lapply(c("pca","umap"), function(x) 
        DimPlot(seurat_obj, reduction=x, group.by="cell_type", combine=TRUE) + 
          theme(legend.position="top") + labs(title = x) +
          scale_color_manual(values=gp_col) +
          guides(color=guide_legend(nrow=3, byrow=TRUE, override.aes=list(size=3))) 
      )
      if (verbose) message("Saving plots.")
      seurat_pdf <- file.path(out_dir,sprintf("seurat_%s.pdf",mc_id))
      pdf(seurat_pdf,height=6,width=6)
      print(plots)
      dev.off()
    }
    
  } else {
    
    # scaling
    if (verbose) message("Scaling data.")
    seurat_obj <- ScaleData(seurat_obj)
    
  }
  
  # saving outputs
  if (!is.null(seurat_rds)) {
    #seurat_rds <- file.path(out_dir,sprintf("seurat_%s.RDS",mc_id))
    if (verbose) message("Saving Seurat object to ",seurat_rds)
    saveRDS(seurat_obj, file=seurat_rds)
    if (verbose) message("Done.")
  }
  
  return(seurat_obj)
  
}


#' Export Metacell solution as `anndata` object
#' 
#' @param mat_id character, id of single cell matrix object in scdb
#' @param mc_id character, id of metacell object in scdb
#' @param gset_id id of markers gene set in scdb (used for building Metacell solution)
#' @param mc2d_id id of metacell 2d projection object
#' @param keep_single_cells logical, keep single cell resolution? Default `TRUE`, if `FALSE`, 
#'   use metacells in place of single cells in seurat object
#' @param mc_ann character, path to file containing metacell annotation, this should be a tsv file
#'   with the three columns containing metacell, cell type and color annotaions
#' @param gene_ann character, path to file containig gene annotations, this should be a tsv file 
#'   with gene ids in the first column
#' @param output_file character, where tosave resulting `h5ad` file (if NULL (default), no file is saved)
#' @param verbose logical, print detailed log messages (default: TRUE)
#'   
mc_export_ann <- function(
  mat_id, 
  mc_id, 
  gset_id = NULL, 
  gstat_id = NULL,
  mc2d_id = NULL,
  mc_ann = NULL, 
  gene_ann = NULL, 
  output_file = NULL, 
  keep_single_cells = TRUE, 
  verbose = TRUE) {
  
  if (verbose) message("Loading metacell objects")
  mat <- scdb_mat(mat_id)
  mc <- scdb_mc(mc_id)
  cells <- names(mc@mc)
  
  # raw counts
  X <- t(mat@mat[,cells])
  
  # cell metadata
  meta_dt <- as.data.table(mat@cell_metadata[cells,], keep.rownames="cell")
  if (!is.null(mc_ann)) {
    # load annotation
    if (verbose) message("Loading metacell annotations")
    ann <- fread(mc_ann, header=TRUE)
    setnames(ann,colnames(ann)[1:3],c("metacell","cell_type","color"))
    ann[, metacell:=as.character(metacell)]
    cell_asgn <- data.table(cell=names(mc@mc),metacell=as.character(mc@mc))
    cell_ann <- merge.data.table(cell_asgn,ann,by="metacell",all.x=TRUE)
    cell_ann[is.na(cell_type),c("cell_type","color"):=list("unassigned","gray")]
    # account for potential names colisions
    for (i in intersect(colnames(meta_dt),colnames(cell_ann)))
      if (i!="cell") setnames(meta_dt,i,paste0("dataset_",i))
    meta_dt_ann <- merge.data.table(cell_ann,meta_dt,by="cell",all.x=TRUE)
    meta_dt_ann[is.na(cell_type),c("cell_type","color"):=list("unassigned","gray")]
    # change colors to hex codes
    meta_dt_ann[!grepl("^#",color), color:=gplots::col2hex(color)]
  } else {
    meta_dt_ann <- meta_dt[,metacell:=as.character(mc@mc[cell])]
  }
  obs_df <- copy(meta_dt_ann); class(obs_df) <- "data.frame"; rownames(obs_df) <- meta_dt_ann$cell; obs_df <- obs_df[cells,]
  
  
  # dim red
  if (!is.null(mc2d_id)) {
    if (verbose) message("Loading dimensionaluty reduction")
    mc2d <- scdb_mc2d(mc2d_id)
    mc2d_sc <- matrix(c(
      mc2d@sc_x[cells],
      mc2d@sc_y[cells]
    ), ncol = 2, byrow = FALSE)
    colnames(mc2d_sc) <- c("x","y")
    rownames(mc2d_sc) <- cells
    obsm <- list(
      X_draw_graph_mc = mc2d_sc
    )
  } else {
    obsm <- NULL
  }

  
  # gene meta data
  genes <- colnames(X)
  gene_dt <- data.table("gene" = genes)
  if (!is.null(gene_ann)) {
    if (verbose) message("Loading gene annotations")
    gann_dt <- fread(gene_ann, header=FALSE)
    setnames(gann_dt,colnames(gann_dt)[1],"gene")
    gene_dt <- merge.data.table(gene_dt, gann_dt, all.x=TRUE)
  }
  gene_dt[,variable_gene:="FALSE"][gene %in% rownames(mc@mc_fp), variable_gene:=TRUE]
  # markers
  if (!is.null(gset_id)) {
    marker_gset <- scdb_gset(gset_id)
    marker_genes <- names(marker_gset@gene_set)
    variable_genes <- intersect(marker_genes,rownames(mat@mat))
    if (verbose) message("Supplied ",length(marker_genes)," marker genes; using ",length(variable_genes)," variable genes.")
    gene_dt[,marker_gene:="FALSE"][gene %in% variable_genes, marker_gene:=TRUE]
  }
  # gene stats
  if (!is.null(gstat_id)) {
    gstat <- scdb_gstat(gstat_id)
    gstat <- setDT(gstat)
    setnames(gstat,colnames(gstat)[1], "gene")
    gene_dt <- merge.data.table(gene_dt, gstat, by="gene", all.x=TRUE)
  }
  var_df <- copy(gene_dt); class(var_df) <- "data.frame"; rownames(var_df) <- gene_dt$gene; var_df <- var_df[genes,]
  
  # anndata object
  ad <- AnnData(X = X, obs = obs_df, var = var_df, obsm = obsm)
  
  # save 
  if (!is.null(output_file)) {
    write_h5ad(ad, output_file)
    if (verbose) message("Saved to ", output_file)
  }

  return(ad)  
}

#' Import Metacell from `anndata` object
#' 
#' @param ad_path character, path to h5ad file
#' @param mat_id character, id of single cell matrix object in scdb
#' @param mc_id character, id of metacell object in scdb
#' 
mc_import_ann <- function(
  ad_path,
  mat_id, 
  mc_id,
  mc_column = "SEACell",
  mc_cols = NULL,
  mc2d_id = NULL,
  graph_id = NULL,
  K = 100,
  obs_2d_key = 'X_umap'
) {
  
  # load anndata object
  ad <- read_h5ad(ad_path)
  
  # make metacell matrix object
  cell_meta <- ad$obs
  gene_meta <- ad$var
  umis <- Matrix::t(ad$X)
  mat <- scm_new_matrix(umis[1:nrow(umis),1:ncol(umis)], cell_meta)
  scdb_add_mat(mat_id, mat)
  message("new matrix added to scdb: ", mat_id)
  
  # get mutacell assignments
  mcs <- as.integer(cell_meta[,mc_column])
  mc_assignment <- structure(mcs, names=rownames(cell_meta))
  
  # make metacell object
  scdb_add_mc(mc_id, tgMCCov(mc = mc_assignment, outliers = character(length = 0), scmat = mat))
  mc <- scdb_mc(mc_id)
  if (!is.null(mc_cols)) {
    if (!is.null(names(mc_cols))) {
      mc@colors <- mc_cols[as.character(1:max(mc@mc))]
    } else {
      mc@colors <- structure(mc_cols, names=1:max(mc@mc))
    }
  } else {
    mc@colors <- structure(rep("white", max(mc@mc)), names=1:max(mc@mc))
  }
  scdb_add_mc(mc_id, mc)
  message("new metacell added to scdb: ", mc_id)
  
  # single cell reduced dim representation
  # redim <- ad$obsm[[obs_2d_key]]
  # sc_x = redim[,1]
  # sc_y = redim[,2]
  
  # build graph
  if (is.null(graph_id)) graph_id = sprintf("%s_graphk%s",mat_id,K)
  d_mat = ad$obsp$connectivities
  colnames(d_mat) <- ad$obs_names
  mcell_add_cgraph_from_distmat(d_mat, graph_id, K, balance = FALSE, k_expand = 0.5*K)
  message("new metacell graph object to scdb: ", graph_id)
  
  # build 2d object
  mc2d_id <- sprintf("%s_2dproj",mc_id)
  mcell_mc2d_force_knn(mc2d_id=mc2d_id, mc_id=mc_id, graph_id=graph_id)
  message("new metacell 2d object to scdb: ", mc2d_id)

}
