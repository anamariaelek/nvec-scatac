require(tgconfig)
require(metacell)
require(ComplexHeatmap)

#' Plot color bar
#'
plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=NULL)
{
  if (!is.null(fig_fn)) {
    .plot_start(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

  if (is.null(show_vals_ind)) {
    show_vals_ind = rep(T, length(cols))
  }
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

  if (!is.null(fig_fn)) {
    dev.off()
  }
}

#' Compute metacell gene footprint directly from tables (no need to use proper mc objects)
#'
#' @param mc a metacell assignment object (equivalent to `mc@mc`)
#' @param us umi matrix (equivalent to `mat@mat`)
#' @param norm_by_mc_size normalize by mean total umis and then multiply by median mean mc size (or at least 1000). This means umis per 1000 molecules.
#' @param min_total_umi consider genes with at least min_total_umi total umis
#'
mc_compute_fp_noobj = function (mc, us, norm_by_mc_size = TRUE, min_total_umi = 10) {
	us = as.matrix(us)
	us = us [ , names(mc) ]
	f_g_cov = rowSums(us) > min_total_umi
	# if (0) {
	# 	mc_cores = get_param("mc_cores")
	# 	doMC::registerDoMC(mc_cores)
	# 	all_gs = rownames(us[f_g_cov, ])
	# 	n_g = length(all_gs)
	# 	g_splts = split(all_gs, 1 + floor(mc_cores * (1:n_g) / (n_g + 
	# 														  	1)))
	# 	fnc = function(gs) {
	# 		.row_stats_by_factor(us[gs, ], mc, function(y) {
	# 			exp(rowMeans(log(1 + y))) - 1
	# 		})
	# 	}
	# 	clust_geomean = do.call(rbind, mclapply(g_splts, fnc, 
	# 											mc.cores = mc_cores))
	# }
	clust_geomean = t(tgs_matrix_tapply(us[f_g_cov, ], mc, 
										function(y) {
											exp(mean(log(1 + y))) - 1
										}))
	rownames(clust_geomean) = rownames(us)[f_g_cov]
	if (norm_by_mc_size) {
		mc_meansize = tapply(colSums(us), mc, mean)
		ideal_cell_size = pmin(1000, median(mc_meansize))
		g_fp = t(ideal_cell_size * t(clust_geomean) / as.vector(mc_meansize))
	}
	else {
		g_fp = clust_geomean
	}
	fp_reg = 0.1
	g_fp_n = (fp_reg + g_fp) / apply(fp_reg + g_fp, 1, median)
	return(g_fp_n)
}


#' Plot 2d projection with more options
#' 
mcell_mc2d_plot_mas = function(
  mc2d_id, legend_pos="topleft", fn_suf="mc", 
  show_mc=TRUE, show_mcid=TRUE, colors=NULL,
  plot_edges=TRUE, min_edge_l=0, edge_w = 1, short_edge_w=0, 
  cell_outline=F, sc_colors=NULL, sc_cex=1, sc_alpha=1, mcp_2d_id_cex=NULL,
  filt_mc=NULL, plot_format="png"
) {
  mcp_2d_height = get_param("mcell_mc2d_height",package="metacell")
  mcp_2d_width = get_param("mcell_mc2d_width",package="metacell")
  mcp_2d_plot_key = get_param("mcell_mc2d_plot_key",package="metacell")
  mcp_2d_cex = get_param("mcell_mc2d_cex",package="metacell")
  if (is.null(mcp_2d_id_cex)) mcp_2d_id_cex = mcp_2d_cex
  mcp_2d_legend_cex = get_param("mcell_mc2d_legend_cex",package="metacell")
  
  mc2d = scdb_mc2d(mc2d_id)
  if(is.null(mc2d)) {
    stop("missing mc2d when trying to plot, id ", mc2d_id)
  }
  mc = scdb_mc(mc2d@mc_id)
  if(is.null(mc)) {
    stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
  }
  if(!is.null(filt_mc)) {
    f_sc = filt_mc[mc@mc[names(mc2d@sc_x)]]
    mc2d@sc_x[!f_sc] = NA
    mc2d@sc_y[!f_sc] = NA
    mc2d@mc_x[!filt_mc] = NA
    mc2d@mc_y[!filt_mc] = NA
  }
  fig_nm = scfigs_fn(fn_suf, ifelse(plot_edges, "2d_graph_proj", "2d_proj"), ext=plot_format)
  #.plot_start(fig_nm, w=mcp_2d_width, h = mcp_2d_height)
  if (plot_format=="png") {
    png(fig_nm, width = mcp_2d_width, height = mcp_2d_height)
  } else if (plot_format=="pdf") {
    pdf(fig_nm, width = mcp_2d_width/72, height = mcp_2d_height/72)
  } else {
    stop("plot_format should be either 'png' or 'pdf'")
  }
  par(bg=NA)
  if(is.null(colors)) {
    cols = mc@colors
  } else {
    cols = colors
  }
  cols[is.na(cols)] = "gray"
  if (is.null(sc_colors)) {
    sc_colors=cols[mc@mc[names(mc2d@sc_x)]]
  } else {
    if (!is.null(names(sc_colors)))
      sc_colors <- sc_colors[names(mc2d@sc_x)]
  }
  if(cell_outline) {
    raster::plot(mc2d@sc_x, mc2d@sc_y, pch=21, bg=alpha(sc_colors,sc_alpha), cex=sc_cex, lwd=0.5,
                 axes=FALSE, frame.plot=TRUE, xlab="", ylab="")
  } else {
    raster::plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=alpha(sc_colors,sc_alpha), cex=sc_cex,
                 axes=FALSE, frame.plot=TRUE, xlab="", ylab="")
  }
  if(show_mc) {
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    if (plot_edges) {
      dx = mc2d@mc_x[fr]-mc2d@mc_x[to]
      dy = mc2d@mc_y[fr]-mc2d@mc_y[to]
      f = sqrt(dx*dx+dy*dy) > min_edge_l
      segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to], 
               lwd=ifelse(f, edge_w, short_edge_w))
    }  
    points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex, col="black", pch=21, bg=cols)
    if(show_mcid) {
      text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_id_cex)
    }  
  }
  if(nrow(mc@color_key)!=0 & mcp_2d_plot_key) {
    key = mc@color_key[ mc@color_key$color %in% mc@colors, ]
    #		if(nrow(key!=0)) {
    if(!is.null(key) & is.vector(key) & nrow(key) != 0) {
      #group	gene	color	priority	T_fold
      gmark = tapply(key$gene, key$group, paste, collapse=", ")
      gcol = unique(data.frame(col=key$color, group=key$group))
      rownames(gcol) = gcol$group
      if(is.vector(gmark)) {
        gmark = gmark[order(names(gmark))]
      }
      if(legend_pos == "panel") {
        dev.off()
        fig_nm = scfigs_fn(mc2d_id, "2d_proj_legend")
        pdf(fig_nm, width = 600/72, height=(length(gmark)*40+400)/72)
        plot.new()
        legend_pos = "topleft"
      }
      legend(legend_pos,
             legend=gsub("_", " ", paste0(names(gmark), ": ", gmark)),
             pch=19, cex=mcp_2d_legend_cex,
             col=as.character(gcol[names(gmark), 'col']), bty='n')
    }
  }
  
  dev.off()
}

#' Plot 2d projection coloured by values
#'
mcell_mc2d_plot_values <- function (
    mc2d_id, mc_values=NULL, sc_values=NULL, name="", show_mc_ids = F, show_legend = T, neto_points = F, 
    max_v = NA, min_v = NA, color_cells = F,
    zero_sc_v = 0, one_sc_v = 1, two_sc_v = 2, filt_mc = NULL, 
    plot_format="png", raster=TRUE, verbose=FALSE) 
{
    height = get_param("mcell_mc2d_gene_height",package="metacell")
    width = get_param("mcell_mc2d_gene_width",package="metacell")
    mc_cex = get_param("mcell_mc2d_gene_mc_cex",package="metacell")
    sc_cex = get_param("mcell_mc2d_gene_cell_cex",package="metacell")
    colspec = get_param("mcell_mc2d_gene_shades",package="metacell")
    if (is.na(max_v)) {
        max_v = get_param("mcell_mc2d_gene_max_v",package="metacell")
        min_v = -max_v
    }
    mc2d = scdb_mc2d(mc2d_id)
    if (is.null(mc2d)) {
        stop("missing mc2d when trying to plot, id ", mc2d_id)
    }
    mc = scdb_mc(mc2d@mc_id)
    if (is.null(mc)) {
        stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
    }
    if (!is.null(filt_mc)) {
        f_sc = filt_mc[mc@mc[names(mc2d@sc_x)]]
        mc2d@sc_x[!f_sc] = NA
        mc2d@sc_y[!f_sc] = NA
        mc2d@mc_x[!filt_mc] = NA
        mc2d@mc_y[!filt_mc] = NA
    }
    if (is.null(mc_values) & !is.null(sc_values)) {
       mc_values <- unlist(tapply(names(mc@mc),mc@mc, FUN=function(x) median(sc_values[x]),simplify=FALSE),use.names=FALSE) 
       if (verbose) message("mc_values: ", paste(mc_values,collaps=" "))
    } else if (is.null(mc_values) & is.null(sc_values)) {
      stop("You need to specify either mc_values or sc_values, or both!")
    }
    x = pmin(pmax(log2(mc_values), min_v), max_v) - 
        min_v
    shades = colorRampPalette(colspec)(100 * (max_v - min_v) + 1)
    mc_cols = shades[round(100 * x) + 1]
    if (verbose) message("mc_cols: ", paste(mc_cols,collaps=" "))
    fig_nm = scfigs_fn(mc2d_id, sub("\\/", "", name), sprintf("%s/%s.values", 
        .scfigs_base, mc2d_id))
    if (plot_format=="png") {
      png(paste0(fig_nm,".png"), w = width * ifelse(show_legend & !neto_points, 1.25, 1), h = height, res=200)
    } else if (plot_format=="pdf") {
      if(raster==TRUE){
        rasterpdf::raster_pdf(paste0(fig_nm,".pdf"), w = width * ifelse(show_legend & !neto_points, 1.25, 1) / 72, h = height / 72)
      } else {
        pdf(paste0(fig_nm,".pdf"), w = width * ifelse(show_legend & !neto_points, 1.25, 1) / 72, h = height / 72)
      }
    } else {
      stop("plot_format should be either 'png' or 'pdf'")
    }
    if (show_legend & !neto_points) {
        layout(matrix(c(1, 1:3), nrow = 2, ncol = 2), widths = c(4, 
            1))
    }
    if (neto_points) {
        par(mar = c(1, 1, 1, 1))
    }
    else {
        par(mar = c(4, 4, 4, 1))
    }
    sc_cols = "gray80"
    if (color_cells & !is.null(sc_values)) {
        cnms = intersect(names(mc2d@sc_x), names(sc_values))
        sc_umi = rep(NA, length(mc2d@sc_x))
        names(sc_umi) = names(mc2d@sc_x)
        sc_umi[cnms] = sc_values[cnms]
        sc_umi[is.na(sc_umi)] = 0
        base_shade = 1 + floor(length(shades) * max_v/(max_v - 
            min_v))
        l_shade = length(shades) - base_shade - 1
        collow = shades[base_shade + floor(l_shade/4)]
        colmid = shades[base_shade + floor(l_shade/2)]
        colhigh = shades[base_shade + floor(3 * l_shade/4)]
        sc_cols = ifelse(sc_umi <= zero_sc_v, "gray80", ifelse(sc_umi <= 
            one_sc_v, collow, ifelse(sc_umi <= two_sc_v, colmid, 
            colhigh)))
       if (verbose) message("sc_cols: ", paste(unique(sc_cols),collaps=" "))
    }
    plot(mc2d@sc_x, mc2d@sc_y, pch = 19, cex = sc_cex, col = sc_cols, 
        xlab = "", ylab = "", main = ifelse(neto_points, "", 
            name), cex.main = mc_cex, bty = ifelse(neto_points, 
            "n", "o"), xaxt = ifelse(neto_points, "n", "s"), 
        yaxt = ifelse(neto_points, "n", "s"))
    points(mc2d@mc_x, mc2d@mc_y, pch = 21, bg = mc_cols, cex = mc_cex)
    if (show_mc_ids) {
        text(mc2d@mc_x, mc2d@mc_y, seq_along(mc2d@mc_y), cex = mc_cex * 
            0.5)
    }
    if (show_legend & !neto_points) {
        par(mar = c(4, 1, 4, 1))
        plot_color_bar(seq(min_v, max_v, l = length(shades)), 
            shades, show_vals_ind = c(1, 100 * max_v + 1, 200 * 
                max_v + 1))
    }
    dev.off()
}

#' Calculate confusion matrix - needed for plotting function
#' 
mc_compute_norm_confu_matrix <- function(
  mc_id, graph_id, max_deg=NULL
){
  cgraph <- scdb_cgraph(graph_id)
  if(is.null(max_deg))	max_deg <- nrow(cgraph@edges)
  confu <- mcell_mc_confusion_mat(mc_id, graph_id, max_deg,ignore_mismatch=T)
  r_confu <- rowSums(confu)
  c_confu <- colSums(confu) 
  norm <- r_confu %*% t(c_confu)
  confu_n <- confu/norm
  confu_n [ is.na(confu_n) ] = 0
  confu_nodiag <- confu_n
  diag(confu_nodiag) <- 0
  confu_n <- pmin(confu_n, max(confu_nodiag, na.rm = TRUE))
  if (nrow(confu_n) > 3 ) {
    confu_n <- pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n)))
  }
  return(confu_n)
}

#' Cluster confusion matrix - needed for plotting function
#' 
mc_confusion_clustering <- function(
  confu_n, clust_method="average"
){
  epsilon <- quantile(confu_n[confu_n != 0], 0.02)
  hc <- hclust(as.dist(-log10(epsilon + confu_n)), clust_method)
  hc$height <- round(hc$height, 6) 
  return(hc)
}

#' Plot confusion matrix and save as raster pdf
#'
mcell_mc_plot_confusion_mat <- function (
    mc_id, graph_id, coc_id = NULL, 
    use_orig_order = FALSE, mc_order = NULL, 
    log_scale = FALSE, label_size = 6,
    fig_fn = NULL, w=7, h=7, plot_format = "pdf"
){
  mc = scdb_mc(mc_id)
  if (is.null(mc)) {
    stop("undefined meta cell object ", mc_id)
  }
  if (is.null(fig_fn)) {
    fig_fn = scfigs_fn(mc_id, sprintf("graph%s_confusion", graph_id), ext = plot_format)
  }
  if (!is.null(graph_id)) {
    if (!is.null(coc_id)) {
      stop("cannot specify both a graph and coclust graph when plotting confusion")
    }
    cgraph = scdb_cgraph(graph_id)
    if (is.null(cgraph)) {
      stop("undefined cgraph object when trying to plot confusion, id ",
           graph_id)
    }
    max_deg = nrow(cgraph@edges)
    confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg,
                                   ignore_mismatch = T)
  } else if (!is.null(coc_id)) {
    coc = scdb_coclust(coc_id)
    if (is.null(coc)) {
      stop("undefined coclust object when trying to plot confusion, id ",
           coc_id)
    }
    max_deg = median(table(mc@mc))/2
    confu = mcell_mc_coclust_confusion_mat(mc_id, coc_id = coc_id,
                                           K = max_deg, ignore_mismatch = T, alpha = 2)
  }
  r_confu = rowSums(confu)
  c_confu = colSums(confu)
  norm = r_confu %*% t(c_confu)
  confu_n = confu/norm
  
  # colors
  if (!use_orig_order) {
    if (is.null(mc_order)) {
      hc = mcell_mc_hclust_confu(mc_id, graph_id = NULL, confu)
      mc_order = hc$order
    }
    confu = confu[mc_order, mc_order]
    confu_n = confu_n[mc_order, mc_order]
    colnames(confu_n) = (1:ncol(confu_n))[mc_order]
    rownames(confu_n) = (1:ncol(confu_n))[mc_order]
  }
  
  # log scale?
  if (log_scale) {
    confu <- log2(1 + confu)
  }
  else {
    confu_nodiag = confu_n
    diag(confu_nodiag) = 0
    confu_n = pmin(confu_n, max(confu_nodiag))
    confu_n = pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n), na.rm = TRUE))
    confu <- confu_n
  }
  
  # heatmap
  mc_cols <- structure(mc@colors, names = seq(1,max(mc@mc)))
  ha_row <- HeatmapAnnotation(
    which = "row",
    metacell = rownames(confu_n),
    col = list("metacell" = mc_cols),
    show_legend = c("metacell" = FALSE),
    show_annotation_name = FALSE,
    border = TRUE
  )
  ha_column <- HeatmapAnnotation(
    which = "column",
    metacell = colnames(confu_n),
    col = list("metacell" = mc_cols),
    show_legend = c("metacell" = FALSE),
    show_annotation_name = FALSE,
    border = TRUE
  )
  col_fun = colorRamp2(
    seq(from = 0, to = max(confu_n), length.out = 6), 
    c("white", "pink", "red", "black", "brown", "orange")
  )
  hm <- Heatmap(
    confu_n, name = "co-clustering",  
    col = col_fun, border = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = label_size),
    column_names_gp = gpar(fontsize = label_size),
    left_annotation = ha_row,
    top_annotation = ha_column
  )
  
  if (!is.null(fig_fn)) pdf(fig_fn, width = w, height = h)
  draw(hm)
  if (!is.null(fig_fn)) dev.off()
}

#' Override pre-set parameters from config file
#'
#' Modified version of tgconfig::set_param() that only overrired parameters known to the specified package 
#' (i.e. doesn't throw an error if there are unknown parameters but ignores them)
override_params <- function (config_file, package = NULL) 
{
  package <- package %||% guess_package(parent.frame(n = 2))
  for (conf_file in config_file) {
    conf <- yaml::yaml.load_file(config_file, eval.expr = TRUE)
    params <- names(conf)
    for (i in 1:length(conf)) {
      tryCatch(
        set_param(params[i], conf[[params[i]]], package = package), error = function(e) NULL
      )
    }
  }
}

#' Convenience function to import kallisto outputs as metacell count objects
#' 
#' @param mat_nm name of the new metacell matrix object (will be saved to the initialised scdb)
#' @param matrix_fn path to kallisto sparse expression matrix (typically .mtx format)
#' @param genes_fn path to list of gene names (used as colnames)
#' @param cells_fn path to a list of cell barcodes/names (used as rownames)
#' 
mcell_import_scmat_kallisto = function(mat_nm, matrix_fn, genes_fn, cells_fn, transpose = FALSE, header_cells = FALSE, header_genes = FALSE) {
  
  require("Matrix")
  require("metacell")
  
  # timestamp
  stamp = runif(1, min = 1e6, max = 9e6)
  
  # load sparse matrix
  mat = Matrix::readMM(file = matrix_fn)
  if (transpose) {
    mat = t(mat)
  }
  colnames(mat) = read.table(cells_fn, header = header_cells)[,1]
  rownames(mat) = read.table(genes_fn, header = header_genes)[,1]
  
  # temporary matrix
  # write.table(mat, file = sprintf("tmp.%s", stamp), sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  saveRDS(mat, sprintf("tmp.%s.RDS", stamp))
  
  # load as metacell object
  # metacell::mcell_import_scmat_tsv(
  #   mat_nm = mat_nm,
  #   dset_nm = mat_nm,
  #   fn = sprintf("tmp.%s", stamp))
  mcell_import_scmat_10x_custom(mat_nm, sprintf("tmp.%s.RDS", stamp))
  
  # remove temporary matrix
  # file.remove(sprintf("tmp.%s", stamp))
  file.remove(sprintf("tmp.%s.RDS", stamp))
  
}

# from https://github.com/tanaylab/metacell/blob/e7176750d5eb7465a23491671f301a055c0b6ee6/R/scmat_10x.r#L67 not modified
# imports custom count matrix from 10x run. Matrix file must be in RDS format. Either sparse or regular matrix format are supported.
mcell_import_scmat_10x_custom = function(mat_nm, matrix_fn = NULL, force = FALSE)
{
  if(!scdb_is_valid()) {
    stop("MCERR - scdb is not initialized, cannot import")
  }
  if(!force & !is.null(scdb_mat(mat_nm))) {
    warning(sprintf("%s already exist in scdb - set force=TRUE to overwrite",mat_nm))
    return(TRUE)
  }
  scdb_add_mat(mat_nm, scmat_read_scmat_10x_custom(matrix_fn, dataset_id=mat_nm))
  return(TRUE)
}
 
#' Modified metacell function to make sure that rownames and colnames are imported properly
#' 
#' @param genes name of the column in `rowData(sce)` that contains gene names
#' @param cells name of the column in `colData(sce)` that contains cell names
#' 
mcell_import_sce_to_mat <- function(sce, counts_slot = "counts", cells = NULL, genes = NULL) {
  umis <- tryCatch(
    expr = SummarizedExperiment::assay(sce, counts_slot),
    error = function(e) {
      stop(paste0("No data in provided assay - ", counts_slot))
  })
  if (!is.null(cells)) 
    colnames(umis) <- colData(sce)[,cells]
  if (!is.null(genes)) 
    rownames(umis) <- rowData(sce)[,genes]
  err <- vector(length=0)
  if (is.null(rownames(umis)))
    err <- "Problem with gene names"
  if (is.null(colnames(umis)))
    err <- c(err,"Problem with cell names")
  if (length(err)>0)
    stop(err)
  md <- as.data.frame(SummarizedExperiment::colData(sce))
  if (is.null(rownames(md)))
    rownames(md) <- colnames(umis)
  scm_new_matrix(umis, md)
}

#' Plot gene/feature statistics (modified function)
#'
#' @param gstat_id id of gene state object
#' @param gset_id optional gene set, will be used to mark select genes if specified
#' @param output_file,width,height,res control plot (default is current device)
#'
mcell_plot_gstats_mod = function(gstat_id, gset_id = NULL, output_file = NULL, max_vm = 4, width = 12, height = 4, res = NA) {
	
  gstat = scdb_gstat(gstat_id)
	if (is.null(gstat)) {
		stop("MC-ERR non existing gstats in mcell_plot_gstats, id ", gstat_id)
	}
	if (!is.null(gset_id)) {
		gset = scdb_gset(gset_id)
		if (is.null(gset)) {
			stop("MC-ERR non existing gset id ", gset_id, " when calling plot gstat")
		}
		marks = names(gset@gene_set)

		set_cols = rep(RColorBrewer::brewer.pal(n=max(3, min(length(gset@set_names), 9)) , "Set1"), times=length(gset@set_names))[1:length(gset@set_names)]
		names(set_cols) = gset@set_names

		cols = unlist(set_cols[gset@gene_set])

	}
  
  # plots
  plotting_function(output_file, width, height, res, EXP = {
    
    # set par
    mfrow_old <- par()$mfrow
    par(mfrow=c(1,3))
    # varmin plot
    vm = pmin(gstat$ds_log_varmean, ifelse(is.null(max_vm), max(gstat$ds_log_varmean), max_vm))
    names(vm) = rownames(gstat)
    plot(log2(gstat$ds_mean), gstat$ds_log_varmean, cex=0.8, pch=1,
      ylim = c(min(vm), max(vm) + 0.25),
        xlab = "log(downsampled mean)", ylab="log2(var/mean downsampled)", 
        col = alpha("gray",0.3),
        main = "varmean")
    if (!is.null(gset_id)) {
      points(log2(gstat[marks,"ds_mean"]),
          vm[marks],
          cex=0.8, pch=19, col=alpha(cols,0.4))
      legend("topleft", legend=names(set_cols), pch=19, cex=0.8, col=set_cols, bty="n")
    }

    # sz correlation plot
    plot(log2(gstat$ds_mean), gstat$sz_cor, cex=0.8, pch=1,
        xlab = "log(downsampled mean)", ylab="sz correlation",
        col = alpha("gray",0.3),
        main = "szcor")
    if (!is.null(gset_id)) {
      points(log2(gstat[marks,"ds_mean"]),
          gstat[marks, "sz_cor"],
          cex=0.8, pch=19, col=alpha(cols,0.4))
      legend("topleft", legend=names(set_cols), pch=19, cex=0.8, col=set_cols, bty="n")
    }

    # top3 plot
    plot(log2(gstat$ds_mean), log2(1 + gstat$ds_top3), cex=0.8, pch=1,
        xlab = "log(downsampled mean)", ylab="log third highest umi (downsamp)",
        col = alpha("gray",0.3),
        main = "top3")
    if (!is.null(gset_id)) {
      points(log2(gstat[marks,"ds_mean"]),
          log2(1 + gstat[marks, "ds_top3"]),
          cex=0.8, pch=19, col=alpha(cols,0.4))
      legend("topleft", legend=names(set_cols), pch=19, cex=0.8, col=set_cols, bty="n")
    }
    # restore par
    par(mfrow=mfrow_old)
  })
}

#' Plot total cell UMIs distribution histogram (modified function)
#'
#' @param mat_id id of the mat to use
#' @param min_umis_cutoff If not null, will mark this cutoff by a line. 
#'   If NA, auto calculate by finding the first bin in the histogram with more cells then its predecessor. 
#' @param bin_for_cutoff bin size for histogram used to auto find min_umis_cutoff
#' @param breaks breaks in histogram, see `?histogram`
#' @param output_file,width,height,res control plot (default is current device)
#'
mcell_plot_umis_per_cell_mod <- function(mat_id, min_umis_cutoff=NA, bin_for_cutoff=50, breaks=NULL, output_file = NULL, width = 4, height = 4, res = NA) {
  mat = scdb_mat(mat_id)
  if (is.null(mat)) {
    stop(sprintf("MCERR: mat with id %s not found", mat_id))
  }
  uc = Matrix::colSums(mat@mat)
  plotting_function(output_file, width, height, res, EXP = {
    if (is.na(min_umis_cutoff)) {
      h = hist(uc, breaks = seq(0, max(uc) + bin_for_cutoff, 
                                by = bin_for_cutoff), plot = FALSE)
      min_umis_cutoff = min(which(diff(h$counts) > 0)) * bin_for_cutoff
    }
    if (!is.null(breaks)) {
      h = hist(log2(uc), 200, xlab = "total cell UMIs (log2)", main = mat_id, breaks=breaks)
    } else {
      h = hist(log2(uc), 200, xlab = "total cell UMIs (log2)", main = mat_id)
    }
    if (!is.null(min_umis_cutoff)) {
      abline(v = log2(min_umis_cutoff), col = "red", lty = 2)
      text(log2(min_umis_cutoff), max(h$count)/2, min_umis_cutoff, 
           pos = 2, col = "red")
    }
    box()
  })
  return(min_umis_cutoff)
}


#' Merge multiple single cell matrix object.
#' 
#' Merge multiple single cell matrix object. 
#' This is meacell function `sca_merge_mats()` modified to explicitly remove ignored genes and cells from `tgScMat` object 
#' (wheras the original function would sometimes fail at setting these slots to `NULL`).
#' Return the merged matrix, with merged meta data and issues an error if there are overlapping cell names between the two matrices.
#' 
#' In case genes sets differs between the matrices, the union is used, with zeros (not NAs!) filling up the missing genes in the respective matrix.
#'
#' @param ... `tgScMat` objects to merge.
#'   Each parameter can be either a single `tgScMat` or a list of `tgScMat` (that will be merged)
#'
#'
sca_merge_mats = function(...)
{
  scmats = list(...)
  scmats = do.call(c, scmats)
  
  if (length(scmats) == 0) {
    return(NULL)
  }
  
  if (!all(sapply(scmats, function(scmat) {"tgScMat" %in% class(scmat)}))) {
    stop("Trying to merge a non-tgScMat object using scm_merge_mats() - if you want to merge from the scdb, call mcell_merge_mats")
  }
    
  res = scmats[[1]]
  
  tryCatch(scm_ignore_cells(res, NULL), error=function(e) {
    mat@ignore_cells <- vector("character",0L)
    mat@ignore_cmat <- new("dgCMatrix")
  })
  tryCatch(scm_ignore_genes(res, NULL), error=function(e) {
    mat@ignore_genes <- vector("character",0L)
    mat@ignore_gmat <- new("dgCMatrix")
  })

  mds = lapply(scmats, function(scmat) {scmat@cell_metadata})
  md_cols = lapply(mds, colnames)
  md_cols_all = Reduce(function(a, b) {c(a, setdiff(b, a))}, md_cols)
  mds = lapply(mds, function(md) {
    missing = setdiff(md_cols_all, colnames(md))
    if (length(missing) > 0) {
      md[, missing] = NA
    }
    return(md[, md_cols_all])
  })
  res@cell_metadata = do.call(rbind, mds)
  
  mats = lapply(scmats, function(scmat) {scmat@mat})
  genes = lapply(scmats, function(scmat) {scmat@genes})
  genes_all = Reduce(function(a, b) {c(a, setdiff(b, a))}, genes)
  if (any(sapply(genes, length) != length(genes_all))) {
    warning("Merged tgScMats have different gene sets. Missing genes will be set to 0.")
  }
  mats = lapply(mats, function(mat) {
    missing <- setdiff(genes_all, rownames(mat))
    if (length(missing) > 0) {
      missing = Matrix::sparseMatrix(i = c(), j = c(), x = 0,
                                     dims = c(length(missing), ncol(mat)),
                                     dimnames = list(missing, colnames(mat)))
      mat = rbind(mat, missing)
    }
    return(mat[genes_all, ])
  })
  res@mat = do.call(cbind, mats)
  
  res@genes = rownames(res@mat)
  res@cells = colnames(res@mat)
  res@ncells = ncol(res@mat)
  res@ngenes = nrow(res@mat)
  res@ignore_genes = vector(l=0)
  res@ignore_cells = vector(l=0)
  res@ignore_cmat =  as(matrix(0, nrow=nrow(res@mat), ncol=0), 'dgCMatrix')
  res@ignore_gmat =  as(matrix(0, nrow=0, ncol=ncol(res@mat)), 'dgCMatrix')
  res@ignore_gcmat =  as(matrix(0, nrow=0, ncol=0), 'dgCMatrix')
  
  return(res)
}


#' Generate a new matrix object with a given ignore cell list
#'
#' @param new_mat_id id of matrix in scdb
#' @param mat_id existing matrix
#' @param ig_cells cells names to ignore
#' @param reverse set this to true to focus on cells instead of ignore them. False by default
#'
#' @export
mcell_mat_ignore_cells_modified = function(new_mat_id, mat_id, ig_cells, reverse=F)
{
	if(is.null(scdb_mat(mat_id))) {
		stop("mcell_mat_ignore_cells called with missing mat id, ", mat_id)
	}
	mat = metacell::scdb_mat(mat_id)

	new_mat = mcell_ignore_cells_modified(mat, ig_cells, reverse)
	metacell::scdb_add_mat(new_mat_id, new_mat)
}

#' Set ignored (i.e. blacklisted) cells
#'
#' Given a list of cells to ignore, this will cancel any previous policy for blacklisting and remove the given cells to the ignore_mat. Downstream algorithm usually ignore these cells altogether, for any purpose including normalization. However, ignored cells can be accessed and analyzed seperately for validation/tests or when they represent some relevant biology (e.g. cell cycle)
#'
#' @param scmat the matrix object
#' @param ig_cell a list of cell names to ignore
#' @param reverse false by default, if this is true the set of cells to ingore is the complement of the given list
#'
#' @export
mcell_ignore_cells_modified = function(scmat, ig_cells, reverse=FALSE) {
  if(is.null(ig_cells) | length(ig_cells) == 0) {
    ig_cells = vector(l=0)
  }
  if(length(scmat@ignore_cells) > 0) {
    scmat@mat = cbind(scmat@mat, scmat@ignore_cmat)
    if(length(scmat@ignore_genes) > 0) {
      scmat@ignore_gmat = cbind(scmat@ignore_gmat, scmat@ignore_gcmat)
    }
  }
  if(!reverse) {
    good_cells = setdiff(colnames(scmat@mat), ig_cells)
    scmat@ignore_cells = intersect(colnames(scmat@mat), ig_cells)
  } else {
    scmat@ignore_cells = setdiff(colnames(scmat@mat), ig_cells)
    good_cells = intersect(colnames(scmat@mat), ig_cells)
    if(length(good_cells) != length(ig_cells)) {
      stop("Some cells to focus on (ignore_cell, reverse=T), are missing from the current matrix, check your list. len(intersect) = ", length(good_cells), " len(ig_cells) = ", length(ig_cells))
    }
  }
  scmat@cells = good_cells
  scmat@ncells = length(good_cells)
  scmat@ignore_cmat = as(scmat@mat[,scmat@ignore_cells], "sparseMatrix")
  if(length(scmat@ignore_genes) > 0) {
    scmat@ignore_gcmat = as(scmat@ignore_gmat[,scmat@ignore_cells], "sparseMatrix")
    scmat@ignore_gmat = as(scmat@ignore_gmat[,good_cells], "sparseMatrix")
  }
  scmat@mat = scmat@mat[,good_cells]
  return(scmat)
}
