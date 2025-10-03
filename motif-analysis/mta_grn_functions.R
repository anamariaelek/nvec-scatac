# libraries
require(universalmotif)
require(GenomicRanges)
require(IRanges)
require(rtracklayer)
require(plyr)
require(tidyr)
require(doParallel)
require(igraph)
require(data.table)

#' Prepare base GRN from motif scanning results
#' (e.g. to be used with CellOracle)
#' 
#' @param gw_scan GRanges object output from `mta_gw_motif_score()`
#' @param peaks_to_genes data.frame with two columns: peak and gene name
#' @param motifs_to_genes data.frame with two columns: motif and TF gene name
#' @param output_file character, tab-separated text file to which the results 
#' will be saved  (recommended to use "tsv.gz" filename extension, so that output 
#' file is saved in compressed format)
#' @param binarize logical, whether to binarize motif scores (default: `FALSE`)
#' For CellOracle set to `TRUE`.
#' @param motifs_per_gene character, determines how to treat multiple motifs 
#' assigned to a single TF gene, one of `c("top", "all")`.
#' If set to `top` (default), the best-scoring motif (across all peaks) is selected.
#' If set to `all`, scores for all motifs assigned to a TF will be summed for each peak - 
#' this only makes sense if you also set `binarize` to `TRUE`.
#' @param nthreads integer, check `parallel::detectCores()`
#' @param chunk_size number of genes to process in single chunk for parallelization.
#' @param return_shape character, one of `c("long","wide")`. See details below.
#' For CellOracle, set this parameter to `wide`.
#' 
#' @details 
#' The function can be run in parallel, but bear in mind that this will require
#' much more RAM, as memory is shared between cores. For example, if each core 
#' needs to hold in memory an initial base GRN for 3000 target genes assigned to 
#' 7k peaks (this expands to over over 1 million rows!), for 700 TF motifs, 
#' each core then needs more than 10GB of RAM for this! So, sometimes it's maybe 
#' better to run the function without paralellization (`nthreads=1`) and let it 
#' run sequentially on chunked input data.
#' 
#' @return data.table with the base GRN.
#' If `return_shape` is set to `long` (default), the resulting data.table has 
#' the following columns: peak_id, gene_sort_name, TF, and binding_score.
#' If `return_shape` is set to `wide`, first two columns of the resulting data.table 
#' are peak_id and gene_short_name, followed by column for each TF, 
#' with binding_score values in rows. 
#' Values in peak_id are composed of peak number, chromosome, motif hit start and 
#' end coordinates, separated by "_".
#' 
mta_prepare_base_grn <- function(
  gw_scan,
  peaks_to_genes,
  motifs_to_genes,
  output_file = NULL,
  binarize = FALSE,
  motifs_per_gene = c("sum", "all", "top")[1],
  nthreads = 1,
  chunk_size = 3000,
  return_shape = c("long", "wide")[1]
) {
  
  require("plyr")
  require("doParallel")
  
  # register cores
  if (nthreads>1) {
    doparallel = TRUE
    registerDoParallel(cores = nthreads)
  } else {
    doparallel = FALSE
  }
    
  
  # inputs
  peaks_to_genes <- unique(as.data.table(peaks_to_genes))
  setnames(peaks_to_genes, c("peak","gene"))
  motifs_to_genes <- unique(as.data.table(motifs_to_genes))
  setnames(motifs_to_genes, c("motif","gene"))
  
  # motif scanning
  if ("GRanges" %in% class(gw_scan)) {
    gw_scan_dt <- setDT(as.data.frame(gw_scan))
    setnames(gw_scan_dt, c("bin_name","name"), c("peak","motif"))
  } else if ("data.frame" %in% class(gw_scan)) {
    gw_scan_dt <- copy(gw_scan)
  }
  gw_scan_dt[,peak_id:=sprintf("%s_%s_%s_%s",peak,seqnames,start,end)]
  gw_scan_dt <- gw_scan_dt[,.(peak,peak_id,motif,score)]
  
  # assignments of peaks to genes
  motif_data <- merge.data.table(peaks_to_genes, gw_scan_dt, allow.cartesian=TRUE)
  motif_data[,empty:=!(sum(score)>0),.(peak_id)]
  motif_data_filtered <- motif_data[empty==FALSE] 
  motif_data_filtered <- motif_data_filtered[, .(peak_id, gene, motif, score)]
  motif_data_filtered <- motif_data_filtered[order(-score), .SD[[1]], .(peak_id, gene, motif)]
  
  # function to build base grn
  fun_grn_building <- function(i) {
    
    # base GRN is a data frame with peaks assigned to genes in rows, and motifs in columns
    message(sprintf("Constructing base GRN from inputs (%s/%s)",i,length(genes_split)))
    genes_subset <- genes_split[[i]]
    grn <- dcast.data.table(
      unique(motif_data_filtered[gene %in% genes_subset]),
      peak_id + gene ~ motif,
      value.var = "score"
    )
    setnames(grn, "gene", "gene_short_name") # harcoded in CellOracle
    message(sprintf(
      "Initial base GRN dims: %s rows, %s columns, %s target genes, %s peaks, %s motifs (%s/%s)", 
      nrow(grn), ncol(grn),
      length(unique(grn$gene_short_name)),
      length(unique(grn$peak)),
      ncol(grn)-2,
      i, length(genes_split)
    ))
    
    # replace NAs with 0s
    for (j in seq_len(ncol(grn)))
      set(grn,which(is.na(grn[[j]])),j,0)
    
    # map motifs to TFs
    message(sprintf("Matching motifs to TFs (%s/%s)",i,length(genes_split)))
    motifs <- colnames(grn)[-c(1:2)]
    missing_motifs <- motifs[!(motifs %in% motifs_to_genes$motif)]
    if (length(missing_motifs)>0) {
      warning(length(missing_motifs), " motif(s) not found in the motifs_to_genes table!")
      motifs_to_genes <- motifs_to_genes[!(motif %in% missing_motifs)]
    } 
    missing_motifs <- unique(motifs_to_genes[!(motif %in% motifs)]$motif)
    if (length(missing_motifs)>0) {
      warning(length(missing_motifs), " motif(s) not found in the motifs scanning results!")
      motifs_to_genes <- motifs_to_genes[!(motif %in% missing_motifs)]
    } 
    motifs_to_genes[,motif:=factor(motif,levels=motifs)]
    setorder(motifs_to_genes,motif)
    
    # replicate table entries for motifs assigned to multiple genes
    motifcols <- as.character(motifs_to_genes$motif)
    motifgenes <- as.character(motifs_to_genes$gene)
    grncols <- c(colnames(grn)[1:2], motifcols)
    grngenes <- c(colnames(grn)[1:2], motifgenes)
    grn_ext <- grn[,..grncols]
    colnames(grn_ext) <- grngenes
    
    # remove motifs that are not mapped to genes 
    unassign_mot <- unique(motifs_to_genes[gene==""]$motif)
    tot_mot <- unique(motifs)
    if (length(unassign_mot)>0) {
      message(sprintf(
        "Removing motifs that are not assigned to genes (%s/%s; %.2f%%) (%s/%s)",
        length(unassign_mot), length(tot_mot), length(unassign_mot)/length(tot_mot)*100, i, length(genes_split)
      ))
      grn_ext[,which(colnames(grn_ext)==""):=NULL]
    }
    
    # binding info for motifs assigned to same gene
    if (motifs_per_gene == "sum") {
      message(sprintf("Summing binding for all motifs for each TFs (%s/%s)",i,length(genes_split)))
    } else if (motifs_per_gene == "all") {
      message(sprintf("Keeping binding for all motifs for each TFs (%s/%s)",i,length(genes_split)))
    } else if (motifs_per_gene == "top") {
      message(sprintf("Selecting motif with highest binding score for each TFs (%s/%s)",i,length(genes_split)))
    }
    genecols <- unique(colnames(grn_ext)[-c(1,2)])
    # mcounts <- sort(table(genecols)) # motifs per gene
    sum_list <- lapply(genecols, function(x) {
      colids <- unique(which(colnames(grn_ext)==x))
      if (motifs_per_gene == "sum" & binarize == TRUE) {
        dt <- data.table(rowSums(grn_ext[,..colids]))
        colnames(dt) <- paste(colnames(dt), grncols[colids], sep = "__")
      } else if (motifs_per_gene == "sum" & binarize == FALSE) {
        stop("When summing info from multiple motifs per gene, you have to set binarize=TRUE")
      } else if (motifs_per_gene == "top") {
        colid <- colids[which.max(colSums(grn_ext[,..colids]))]
        dt <- grn_ext[,..colid]
        colnames(dt) <- paste(colnames(dt), grncols[colid], sep = "__")
        dt
      } else if (motifs_per_gene == "all") {
        dt <- grn_ext[,..colids]
        colnames(dt) <- paste(colnames(dt), grncols[colids], sep = "__")
        dt
      }
    })
    sum_mat <- do.call(cbind, sum_list)
    grn_sum <- as.data.table(as.data.frame(
      sum_mat
    ))
    grn_sum[, c("peak_id", "gene_short_name") := list(
      grn_ext$peak_id, grn_ext$gene_short_name
    )]
    setcolorder(grn_sum, c("peak_id", "gene_short_name"))
    message(sprintf(
      "Final base GRN dims: %s rows; %s columns, %s target genes, %s peaks, %s motifs (%s/%s)", 
      nrow(grn_sum), ncol(grn_sum),
      length(unique(grn_sum$gene_short_name)),
      length(unique(grn_sum$peak)),
      ncol(grn_sum)-2,
      i, length(genes_split)
    ))
    
    # binarize
    if (binarize == TRUE) {
      for (j in 3:ncol(grn_sum))
        set(grn_sum,which(grn_sum[[j]]>0),j,1)
    }
    
    return(grn_sum)
  }
  
  # chunk work because the resulting GRN is too large to fit in memory
  genes <- unique(motif_data_filtered$gene)
  ids <- 1:ceiling(length(genes)/chunk_size)
  genes_split <- tapply(genes, rep(ids,each=chunk_size)[1:length(genes)], function(x) x)
  base_grn_list = plyr::alply(.data = ids, .margins = 1, .fun = fun_grn_building, .parallel = doparallel)
  
  # combine results
  base_grn <- rbindlist(base_grn_list, use.names = TRUE, fill = TRUE) 
  for (j in seq_len(ncol(base_grn)))
    set(base_grn,which(is.na(base_grn[[j]])),j,0)
  
  # format of output data.table
  if (return_shape=="long") {
    tryCatch({
      # sometimes this fails because molten dt would exceed the vector length limit of 2^31-1
      # (https://stackoverflow.com/a/48676389) so we set 0s to NAs and ignore them to have less rows
      for (j in seq_len(ncol(base_grn)))
        set(base_grn,which(base_grn[[j]]==0),j,NA)
      base_grn <- melt.data.table(base_grn, id.vars = c("peak_id","gene_short_name"), variable.name = "TF", value.name = "binding_score", na.rm = TRUE)
      #base_grn <- base_grn[binding_score>0]
    },
    error = function(e) {
      warning("Could not convert to long format (", e, "), returning wide.")
      base_grn
    })
  } else if (return_shape=="wide") {
    # nothng to do
  } else {
    warning("Unrecognized format (it should be either wide or long); returning output data in the wide format.")
  }
  
  # save
  if (!is.null(output_file)) {
    fwrite(base_grn, output_file, sep="\t", col.names=TRUE)
    message(sprintf("Done, saved output to %s", output_file))
  } else {
    message("Done")
  }
  
  return(base_grn)
}

#' Calculate binding energy for all pairs in base GRN
#'
#' @param base_grn data.frame output by `mta_prepare_base_grn()`
#' @param output_file character, tab-separated text file to which the results 
#' will be saved 
#' @param quant_norm logical, quantile normalize motif scores (default: TRUE) - 
#' this puts all score to the same range and removes bias towards longer motifs 
#' (because most motifs scores are not independent of the motif length)
#' @param aggregate_function function to aggregate multiple peaks scores for one TF,
#' should be either `sum` or `max` or `mean`
#' @param peaks_stats optionally, a data.frame with three columns containing peak, 
#' it's variability and distance weights (default: NULL, no weighting is performed)
#'
#' @details
#' For each transcription factor $TF$ - target gene $TG$ pair in the base GRN,
#' calculates TF binding energy $B_{TF,TG}$ as a sum across target gene's
#' regulatory regions (peaks) of TF motif binding scores $S_{TF,peak}$, optionally
#' weighted by peak distance from TSS and peak variability.
#'
#' $B_{TF,TG} = \sum_{peaks} w_{dist} \space w_{var} \space motif_{peak}$$
#'
#' @return data.table with three columns: TF, gene_short_name, and binding_energy
#'
mta_grn_binding_energy <- function(
  base_grn,
  output_file = NULL,
  quant_norm = TRUE,
  aggregate_function = sum,
  peaks_stats = NULL
) {
  
  # base grn input
  base_grn <- as.data.table(base_grn)
  # check format
  if (!all(c("TF","binding_score") %in% colnames(base_grn))) {
    base_grn <- melt.data.table(base_grn, id.vars = c("peak_id","gene_short_name"), variable.name = "TF", value.name = "binding_score")
    # filter out 0s
    base_grn <- base_grn[binding_score>0]
  }
  
  if (quant_norm==TRUE) {
    base_grn[, binding_score := ecdf(.SD$binding_score)(binding_score), TF]
  }
    
  # calculate binding energy
  if (!is.null(peaks_stats)) {
    # extract only peak name from peak_id (also contains coordinates)
    base_grn <- copy(base_grn)[,peak:=str_extract(peak_id,"peak\\d+")]
    base_grn <- merge.data.table(base_grn, peaks_stats, by="peak", all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
    peaks_stats <- unique(as.data.table(peaks_stats))
    setnames(peaks_stats, c("peak",weigth_cols))
    weight_cols = colnames(peaks_stats)[-1]
    for (i in weight_cols) {
      base_grn[,binding_score:=get(i)*binding_score]
    }
  }
  binding_en <- base_grn[, .(
    binding_energy = do.call(
      aggregate_function, list(binding_score)
    )
  ), .(TF, gene_short_name)]
  
  # save
  if (!is.null(output_file))
    fwrite(binding_en, output_file, sep="\t", col.names=TRUE)
  
  return(binding_en)
}

#' Build a linear regularized model to predict gene expression from  TF binding 
#' and TF expression
#' 
#' @param binding data.frame with three columns: tf, gene, and binding energy for tg-gene pair
#' @param expression metrix of expression values, genes in rows, cells in columns
#' @param output_file character, tab-separated text file to which the results 
#' will be saved  (recommended to use "tsv.gz" filename extension, so that output 
#' file is saved in compressed format)
#' 
#' @return data.table with the folowing columns: gene, TF, model coefficient,  
#' coefficient variance and p value
#' 
mta_grn_model <- function(
  binding,
  expression,
  scale_cell_size = TRUE,
  quant_norm_gene = FALSE,
  output_file
) {
  
  require(reticulate)
  
  # binding data
  binding <- as.data.table(binding)
  setnames(binding, c("TF", "gene_short_name", "binding_energy"))

  # scale binding data
  binding[,binding_energy_scaled:=scale(binding_energy, center = FALSE)]

  # binding data in matrix format (genes x TFs)
  binding_scaled_dt <- dcast.data.table(binding, gene_short_name~TF, value.var="binding_energy_scaled")
  binding_scaled <- as.matrix(binding_scaled_dt[,-1])
  rownames(binding_scaled) <- binding_scaled_dt[,1][[1]]
  binding_scaled[is.na(binding_scaled)] <- 0

  # scale expression data
  if (scale_cell_size==TRUE) {
    expression_scaled <- Matrix::t(expression)/Matrix::colSums(expression)
  } else {
    expression_scaled <- Matrix::t(expression)
  }
  if (quant_norm_gene==TRUE) {
    expression_scaled_norm <- apply(expression_scaled, 2, function(x) ecdf(x)(x))
    rownames(expression_scaled_norm) <- rownames(expression_scaled)
  } else {
    expression_scaled_norm <- expression_scaled
  }
  
  # tf genes as features
  tf_features <- intersect(colnames(binding_scaled),colnames(expression_scaled_norm))
  tf_expression_scaled_norm <- expression_scaled_norm[,tf_features]
  binding_scaled <- binding_scaled[,tf_features]
  
  # genes as response variable
  expressed_genes <- names(which(Matrix::colSums(expression_scaled_norm)>0))
  genes <- intersect(rownames(binding_scaled),expressed_genes)
  
  # check genes
  if (all(all.equal(colnames(binding_scaled), colnames(tf_expression_scaled_norm))!=TRUE))
    stop("Genes are not the same or not in the same order!")
  
  np <- reticulate::import("numpy")
  # sc <- reticulate::import("scipy.stats")
  sk <- reticulate::import("sklearn.linear_model")
  
  # models for each gene
  .grn_model <- function(gene) {
    # message(gene)
    # binding data for potential TFs
    train_binding <- binding_scaled[gene,]
    # binding x expression
    x <- Matrix::t(Matrix::t(tf_expression_scaled_norm) * train_binding)
    # remove TFs which are all 0s
    coef_names <- colnames(x)[Matrix::colSums(x)>0]
    # input for model
    x <- as.matrix(x[,coef_names,drop=FALSE])
    y <- as.matrix(expression_scaled_norm[rownames(x),gene,drop=FALSE])
    reg <- sk$BayesianRidge()
    reg <- reg$fit(x, np$ravel(y))
    coefs <- reg$coef_
    coef_vars <- diag(reg$sigma_)
    # p values
    # p = sc$norm$cdf(x=0, loc=abs(coefs), scale=np$sqrt(coef_vars))*2
    p = 2*(1-pnorm(q=abs(coefs), mean=0, sd=sqrt(coef_vars)))
    q = p.adjust(p = p, method = "BH")
    data.table(TF=coef_names, coef=coefs, coef_var=coef_vars, pval=p, qval=q)
  }
  
  coef_list <- lapply(genes, function(x) 
    tryCatch(.grn_model(x), error=function(e) {
      message(e)
      NULL
    })
  )
  names(coef_list) <- genes
  coef_dt <- rbindlist(coef_list, idcol = "gene")
  
  # save
  if (!is.null(output_file))
    fwrite(coef_dt, output_file, sep="\t", col.names=TRUE)
  
  return(coef_dt)
}


#' Calculate network metrics and cartography
#' 
#' @param coef_dt data.frame with at least the following columns: gene and TF 
#' in the first two columns, and weights in another column
#' @param weight_column integer, which column of `coef_dt` contains weights
#' @param output_file character, tab-separated text file to which the results 
#' will be saved 
#' 
#' @return data.frame
#' 
mta_network_scores <- function(coef_dt, weight_column=3, output_file=NULL) {

  message("Graph from data")
  g <- graph.data.frame(coef_dt, directed = TRUE)
  E(g)$weight <- coef_dt[[weight_column]]
  
  message("Nodes")
  res <- igraph::degree(g, mode = "all")
  res <- as.data.table(res, keep.rownames = "gene")
  setnames(res, c("gene","degree_all"))
  
  # degrees
  res[, degree_in := igraph::degree(g, mode = "in")]
  res[, degree_out := igraph::degree(g, mode = "out")]
  
  # clustering coefficient
  message("Clustering coef")
  tryCatch(
    res[, clustering_coefficient := transitivity(g,type="local",isolates="zero")],
    error = function(e) {
      warning("Could not calculate clustering coefficient", e)
    }
  )
  # clustering coefficient for weighted network
  tryCatch(
    res[, clustering_coefficient_weighted := transitivity(g,vids=V(g),type="weighted",isolates="zero")],
    error = function(e) {
      warning("Could not calculate clustering coefficient for weighted network", e)
    }
  )
  
  # degree centrality
  message("Centrality")
  res[, degree_centrality_all := degree_all / (vcount(g) - 1)]
  res[, degree_centrality_in := degree_in / (vcount(g) - 1)]
  res[, degree_centrality_out := degree_out / (vcount(g) - 1)]
  # betweenness centrality
  res[, betweenness_centrality := betweenness(g)]
  # closeness centrality
  res[, closeness_centrality := closeness(g)]
  # eigenvector centrality
  res[, eigenvector_centrality := evcent(g)$vector]
  # page rank
  res[, page_rank := page.rank(g)$vector]
  
  res[, assortative_coefficient := assortativity.degree(g)]
  res[, average_path_length := average.path.length(g)]
  
  # community detection (non-overlapping)
  message("Random walk community detection")
  data <- walktrap.community(g, modularity=TRUE)
  res[, community_random_walk := data$membership]
  
  get_scores_cartography <- function(g){
    adj <- get.adjacency(g,sparse=FALSE)
    score <- rnetcarto::netcarto(adj)[[1]]
    rownames(score) <- score$name
    return(score)
  }
  
  message("Cartography scores")
  result <- get_scores_cartography(g)
  result <- as.data.table(result, keep.rownames="gene")[,name:=NULL]
  
  result_dt <- merge.data.table(res, result, by="gene", all.x=TRUE, sort=FALSE)
  
  if (!is.null(output_file)) {
    fwrite(result_dt, output_file, sep="\t", col.names = TRUE)
    message("Saved to ", output_file)
  }
    
  message("Done")
  return(result_dt)
}

#' Plot graph
#' @param expression named numeric vector
#' @param gene_annotation data.frame with gene IDs in the first column, and any annotations in the following columns'
#'   If there is a column named "og", it is assumed these are ortogroups, and node labels will be extracted from this column.
#'   Otherwise, labels will be taken from `label_column`
#' @param label_column character, name of the column in gene_annotation to be used for node labels
mta_graph <- function(
    df, nodes,
    expression = NULL, 
    min_fc_scale = 0, max_fc_scale = 2, 
    min_fc_filt_gene = -Inf, max_fc_filt_gene = Inf,
    min_fc_filt_tf = -Inf, max_fc_filt_tf = Inf,
    gene_annotation = NULL, tf_annotation = NULL, 
    label_column = "og", show_label = TRUE, 
    label_TFs_only = TRUE, label_TFs_all=TRUE, 
    label_expression_thrs = -Inf, fill_labels_with_gene_ids = FALSE,
    layout = c('fr', 'stress', 'lgl', 'graphopt')[1], 
    fill = c('TF', 'expression')[1], 
    output_file = NULL, width = 10, height = 10, return_data = FALSE
) {
  # graph data.table
  dt <- copy(df)[,1:3]
  setnames(dt, c("TF","gene","weight"))
  
  # graph nodes
  nd <- copy(nodes)
  setnames(nd, colnames(nd)[1], "gene")
  
  # intersect node genes
  node_genes <- unique(intersect(nd$gene, c(dt$TF, dt$gene)))
  dt <- dt[TF %in% node_genes & gene %in% node_genes]
  nd <- nd[gene %in% node_genes]
  
  # gene annots
  if (!is.null(gene_annotation)) {
    ga <- copy(gene_annotation)
    setnames(ga, colnames(ga)[1], "gene")
  } else {
    ga <- NULL
  }
  if (!is.null(tf_annotation)) {
    ta <- copy(tf_annotation)
    setnames(ta, colnames(ta)[1], "gene")
    tfs <- ta$gene
    if (!is.null(gene_annotation)) {
      ga <- rbindlist(list(
        ga[!gene %in% ta$gene],
        ta
      ))
    } else {
      ga <- ta
    }
  } else {
    tfs <- dt$TF
  }
  nd[,TF:=FALSE][gene %in% tfs, TF:=TRUE]
  
  # add expression
  if (!is.null(expression) & length(expression)>0) {
    nd[,mcfp:=expression[gene]]
    # filter genes by expression
    nd <- nd[(TF==FALSE & (mcfp>min_fc_filt_gene & mcfp<max_fc_filt_gene)) | TF==TRUE]
    # filter TFs by expression
    nd <- nd[(TF==TRUE & (mcfp>min_fc_filt_tf & mcfp<max_fc_filt_tf)) | TF==FALSE]
    # scale expression
    nd[,mcfp:=pmax(min_fc,pmin(mcfp,max_fc))]
    # update graph data frame and nodes
    dt <- dt[TF %in% nd$gene & gene %in% nd$gene]
    nd <- nd[gene %in% c(dt$gene,dt$TF)]
  }
  
  # labels
  nd[,label:=gene][,show_label:=show_label]
  if (!is.null(expression) & length(expression)>0)
    nd[mcfp<label_expression_thrs, show_label:=FALSE]
  
  # add gene annotaitons
  if (!is.null(ga)) {
    nd <- merge.data.table(nd, ga, by="gene", all.x = TRUE, sort = FALSE)
    
    # construct labels to show
    if (label_column == "og") {
      nd[,label:=substr(str_remove(og,"[A-z0-9_-]*\\.HG\\d+\\.\\d+(?=:):"),1,25)]
      nd[grepl("ribosom",label,ignore.case=TRUE),show_label:=FALSE]
      nd[grepl("Heat shock",label,ignore.case=TRUE),show_label:=FALSE]
      nd[grepl("cDNA",label,ignore.case=TRUE),show_label:=FALSE]
    } else {
      setnames(nd,label_column,"label")
      nd[,label:=substr(label,1,25)]
    }
    
  } else {
    nd[,label:=substr(gene,1,25)]
  }
  
  # show TF labels?
  if (label_TFs_only==TRUE) {
    nd[TF==FALSE & show_label==TRUE, show_label:=FALSE]
  }
  if (label_TFs_all==TRUE) {
    nd[TF==TRUE, show_label:=TRUE]
  }
  
  # replace missing labels with gene ids
  nd[show_label==TRUE & is.na(label), label:=gene]
  nd[show_label==TRUE & label=="NA", label:=gene]
  if (fill_labels_with_gene_ids) {
    nd[show_label==TRUE & label=="", label:=gene]
  }
  
  # build graph
  net <- graph_from_data_frame(d = dt, vertices = nd, directed = TRUE)
  net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE) 
  
  # degree as node size and add labels
  V(net)$size <- igraph::degree(net)
  V(net)$label <- nd$label
  
  # TFs vs non-TFs in graph
  nodecols <- structure(c("goldenrod","firebrick1"), names=c(FALSE, TRUE)) 
  nodeshapes <- structure(c(21,22), names=c(FALSE, TRUE)) 
  heatmap_colors <- c("#ffe8bd","orange","orangered2","#520c52")
  col_fun <- colorRampPalette(colors = heatmap_colors)
  
  # graph
  gg <- ggraph(net, layout = layout) +
    geom_edge_link0(aes(edge_width = weight), edge_colour = "grey66")
  if (fill == "TF") {
    gg <- gg + 
      geom_node_point(shape=21, aes(size = size, fill = TF)) +
      scale_fill_manual(values=nodecols)
  } else if (fill == "expression") {
    gg <- gg + 
      geom_node_point(aes(size = size, shape = TF, fill = mcfp)) +
      scale_fill_gradientn(colours = heatmap_colors) +
      scale_shape_manual(values = nodeshapes)
  }
  if (show_label == TRUE) {
    gg <- gg + 
      geom_node_text(aes(filter = show_label==TRUE, label = label), repel = TRUE, fontface="bold", color="black")
  }
  gg <- gg + 
    scale_edge_width(range = c(0.05,1)) +
    scale_size(range = c(2,8)) +
    theme_graph(base_family = "Helvetica") +
    theme(legend.position = "right") +
    labs(caption=layout)
  
  # save plot
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    print(gg)
    dev.off()
  }
  
  # data
  data_list <- list(dt= dt, nodes = nd)
    
  # return
  if (return_data==TRUE) {
    data_list
  } else {
    gg
  }
}
