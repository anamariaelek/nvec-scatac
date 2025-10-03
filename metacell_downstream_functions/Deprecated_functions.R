#' DEPRECATED
#' Plot heatmap of gene expression
#'
#' Plots heatmap of gene expression fold change for metacells (and optionally,
#' for single cells).
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param black_list character, blacklisted genes
#' @param output_file output file name without extension, all output files
#'    will have this prefix; if NULL (default), will not save heatmap to file
#' @param height numeric, height of the heatmap image saved to output file
#' @param width numeric, width of the heatmap image saved to output file
#' @param res numeric, resolution of the output file
#' @param plot_format character, format of the output image file, either
#'    "png" (default) or "pdf"
#' @param print_heatmap logical, whether to draw heatmap to console
#' @param sub_list_mc character, subset of metacells to plot (default: NULL)
#' @param plot_sc logical, whether to also plot a heatmap of single cells
#'    (default: FALSE)
#' @param plot_sc_width_cex numeric, factor by which to scale single cell plot
#'    width (default: 2)
#' @param plot_sc_height_cex numeric, factor by which to scale single cell plot
#'    height (default: 1)
#' @param gene_list character, list of genes to plot, if NULL (default) ...
#' @param order_genes logical, whether to cluster genes (default: TRUE)
#' @param highlight_genes logical or character, whether to highlight top expressed
#'    and annotted genes on the heatmap; if character, a set of genes which to
#'    highlight on the heatmap
#' @param gene_annot_file character, file path to the file containig gene annotations,
#'    it should have three tab separated columns containing gene ID, pfam architecture
#'    and any additional annotation in the last column
#' @param tfs_annot_file character, file path to the file containig TFs annotations,
#'    minimum required is one column with the gene IDs; if this file is provided, TF
#'    genes annotations will be highlighted in red on the heatmap
#' @param annot_header logical, gene annotation file has column names?
#' @param show_gene_names logical, show gene names as heatmap row names
#'    (default: FALSE)
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param gene_chr_limit numeric, limit gene annotations to given number of charaters
#' @param clust_ord character, metacells in the order in which they should be
#'    plotted; if cluster order is not specified (default: NULL), it is
#'    determined by hierarchical clustering
#' @param clust_col character, colours to asssign to clusters, should be either
#'    a named vecor where names are names of metacells, or an unamed vector
#'    ordered in the same order as specified by clust_ord; if NULL (default),
#'    cluster colour annotation bar is not printed
#' @param clust_bars numeric, optional values to be ploted as bars annotation on top of
#'    heatmap columns (default:  NULL)
#' @param clust_anno_size unit, height of the column annotation bar (default: `unit(15,"mm")`)
#' @param show_mc_names logical, show metacell names as heatmap column names
#'    (default: TRUE) ~NOT IMPLEMENTED~
#' @param mc_font_size numeric, size of the metacell names (default: 4)
#' @param mc_label_cex numeric, scalling of the metacell names on sc plot (default: 2)
#' @param per_clust_genes integer, how many genes per cluster to aim to show in the heatmap
#'    (default: 20)
#' @param gene_min_fold numeric, minimum fold change for a gene to be considered for plotting
#'    (default: 2)
#' @param smoothen integer, width of the window used to calculate rolling mean
#'    of expression in single cells (default: 5)
#' @param max_expression_fc,max_expression_fc_sc numeric, max expression value
#'    to scale metacell and single cell heatmap coloring to (default: 5)
#' @param transverality_N integer, number of metacells in which a gene can be highly expressed (>1.4)
#'    to be considered for plotting, by default this is the total number of metacells
#' @param transv_excluded_mc character, metacells to be excluded in transversality calculation
#'    (default: NULL)
#'
#'
scp_plot_cmod_markers <- function(
	mc_object, mat_object, black_list=c(),
	output_file = NULL, plot_format = "png", height = 8000, width = 3000, res = NA,
	print_heatmap = FALSE, show_heatmap_legend = TRUE,
	sub_list_mc = NULL, plot_sc = FALSE, plot_sc_width_cex = 2, plot_sc_height_cex = 1,
	
	gene_list = NULL, order_genes = TRUE, highlight_genes = NULL,
	gene_annot_file = NULL, tfs_annot_file=NULL, annot_header = FALSE,
	show_gene_names = FALSE, gene_font_size = 4, gene_chr_limit = 100,
	
	clust_ord = NULL, clust_col = NULL, clust_bars = NULL, show_clust_borders=TRUE,
	clust_anno_size = unit(15,"mm"),
	show_mc_names = TRUE, mc_font_size = 4, mc_label_cex = 2,
	column_anno_padding=unit(10,"mm"), row_anno_padding=unit(5,"mm"), dimname_padding=unit(5,"mm"),
	
	per_clust_genes = 20, gene_min_fold = 2, smoothen = 5, max_expression_fc = 5, max_expression_fc_sc = 5,
	transverality_N = ncol(mc_object@mc_fp), transv_excluded_mc = NULL,
	verbose=FALSE
) {
	# load gene annotations
	if (!is.null(gene_annot_file) & "character" %in% class(gene_annot_file)) {
		annot = read.table(gene_annot_file, header=annot_header, sep="\t", fill=TRUE, quote="", row.names=1)
	} else if (!is.null(gene_annot_file) & "data.frame" %in% class(gene_annot_file)) {
		annot = gene_annot_file
	}
	
	# directory to save output files to
	if (is.null(output_file)) {
		outdir <- getwd()
	} else {
		outdir <- dirname(output_file)
	}
	
	# expression matrix
	if (is.null(sub_list_mc)){
		niche_geomean_n= mc_object@mc_fp
	}else{
		niche_geomean_n= mc_object@mc_fp[,sub_list_mc]
		clust_ord=sub_list_mc
	}
	
	# select genes for plotting
	if (is.null(gene_list)){
		# exclude genes with fc < gene_min_fold
		genes=unique(as.vector(unlist(apply(niche_geomean_n, 2, function(x) names(head(sort(-x [x > gene_min_fold]),n=per_clust_genes))))))
		# exclude blacklisted genes
		genes=setdiff(genes, black_list)
		# exclude genes with transversality > transverality_N
		transversal_genes=names(which(
			apply(
				niche_geomean_n[,setdiff(as.character(colnames(niche_geomean_n)),transv_excluded_mc)],
				1,
				function(x) sort(x,decreasing=TRUE)[transverality_N] > 1.4
			)
		))
		genes=setdiff(genes, transversal_genes)
	} else {
		genes=unique(as.vector(unlist(apply(niche_geomean_n, 2, function(x) names(x[x > gene_min_fold])))))
		# plot only genes in gene list, if it is specified
		genes=intersect(gene_list,genes)
	}
	
	genes=intersect(genes,rownames(niche_geomean_n))
	message("Will use ",length(genes)," genes")
	
	mat_niche <- niche_geomean_n[genes,]
	
	# if cluster order is not specified, do hierarchical clustering
	if (is.null(clust_ord)) {
		message("Recomputing cell ord")
		hc1 = hclust(dist(cor(mat_niche,method="pearson")), "ward.D2")
		clust_ord = as.character(hc1$order)
		# write.table(
		#   clust_ord,
		#   file=file.path(outdir, "tmp_cell_clusts_ordered_by_scp_markers_plot.txt"),
		#   quote=FALSE,col.names=FALSE,row.names=FALSE
		# )
		# png(
		#   file.path(outdir,"tmp_cell_clusts_ordered_by_scp_markes_tree.png"),
		#   height=500,width=1000
		# )
		# plot(hc1,xlab="",xaxt='n',hang=-1,ylab="",main="",cex=1)
		# dev.off()
		scr_tmp_niche_order <- as.character(hc1$order)
	}
	
	# barplot annotation
	if (!is.null(clust_bars)) {
		if (is.null(names(clust_bars)))
			names(clust_bars) <- as.character(clust_ord)
		anno_bar <- clust_bars[as.character(clust_ord)]
	}
	
	# if gene order is TRUE, order genes
	if (order_genes){
		message("Ordering genes")
		gene_ord = order(apply(mat_niche[,as.character(clust_ord)],1,function(x) which.max(rollmean(x,1))))
	} else {
		gene_ord= 1:length(genes)
	}
	gene_ord <- rev(gene_ord)
	# write.table(
	#   genes[gene_ord],
	#   file=file.path(outdir,"tmp_markers_ordered.txt"),
	#   quote=FALSE,col.names = FALSE,row.names=FALSE
	# )
	
	# gene labels
	if (!is.null(gene_annot_file)) {
		
		gene_labels_0 <- genes[gene_ord]
		
		message("Genes: ", head(gene_labels_0), "...")
		
		gene_labels_1 <- as.character(annot[genes[gene_ord],2])
		message("Gene labels: ", head(gene_labels_0), "...")
		bad_labels <- gene_labels_1 %in% c("","-"," ") | is.na(gene_labels_1)
		message(sum(bad_labels), " bad gene labels")
		gene_labels_1[bad_labels] <- gene_labels_0[bad_labels]
		long_labels <- nchar(gene_labels_1) > gene_chr_limit
		message(sum(long_labels), " long gene labels")
		gene_labels_1[long_labels] <- paste0(substr(gene_labels_1[long_labels],1,gene_chr_limit - 3),"...")
		names(gene_labels_1) <- genes[gene_ord]
		
		gene_labels_3 <- as.character(annot[genes[gene_ord],1])
		bad_labels <- gene_labels_3 %in% c("","-"," ") | is.na(gene_labels_3)
		gene_labels_3[bad_labels] <- genes[gene_ord][bad_labels]
		long_labels <- nchar(gene_labels_3) > gene_chr_limit
		gene_labels_3[long_labels] <- paste0(substr(gene_labels_3[long_labels],1,gene_chr_limit - 3),"...")
		gene_labels_2 <- ifelse(gene_labels_0 == gene_labels_3, gene_labels_0, paste(gene_labels_0,gene_labels_3,sep=" "))
		names(gene_labels_2) <- genes[gene_ord]
		
	} else {
		
		gene_labels_0 <- genes[gene_ord]
		gene_labels_1 <- genes[gene_ord]
		gene_labels_2 <- genes[gene_ord]
		
	}
	
	gene_font_col <- rep("black",length(gene_labels_0))
	if (!is.null(tfs_annot_file)) {
		tfs_ann <- fread(tfs_annot_file,header=FALSE)
		gene_font_col[gene_labels_0 %in% tfs_ann[[1]]] <- "red"
	}
	
	if (!is.null(highlight_genes)) {
		if (class(highlight_genes) == "logical") {
			if (highlight_genes == TRUE) {
				highlight_genes <- gene_labels_0[gene_labels_1 != gene_labels_0]
			}
		}
		highlight_ids <- match(highlight_genes,gene_labels_0)
		highlight_genes <- highlight_genes[gene_labels_1[highlight_ids] != gene_labels_0[highlight_ids]]
		
		gids <- which(gene_labels_0 %in% highlight_genes)
		if (length(gids > 0)) {
			gene_labels_1[!gene_labels_0 %in% highlight_genes] <- ""
			gene_labels_2[!gene_labels_0 %in% highlight_genes] <- ""
		} else {
			gids <- which(gene_labels_0)
		}
		
	}
	
	
	
	########################### PLOT METACELL PROFILE ###########################
	
	message("Plotting metacell expression")
	
	# heatmap
	ht_opt(
		COLUMN_ANNO_PADDING=column_anno_padding,
		ROW_ANNO_PADDING=row_anno_padding,
		DIMNAME_PADDING=dimname_padding
	)
	#print(ht_opt())
	mat1 <- pmin(log2(niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + 1),max_expression_fc)
	if (show_gene_names) {
		if (length(highlight_genes) > 1) {
			message("Gene annots highlights")
			row_ha_right <- HeatmapAnnotation(
				which = "row", simple_anno_size = unit(10,"mm"),
				gene = anno_mark(
					which="row", side="right", at=gids, labels=gene_labels_1[gids],
					labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]),
					extend=unit(10, "mm")
				)
			)
			row_ha_left <- HeatmapAnnotation(
				which = "row", simple_anno_size = unit(10,"mm"),
				gene = anno_mark(
					which="row", side="left", at=gids, labels=gene_labels_2[gids],
					labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]),
					extend=unit(10, "mm")
				)
			)
		} else {
			message("Gene annots")
			if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
			row_ha_right <- HeatmapAnnotation(
				which = "row", simple_anno_size = unit(10,"mm"),
				gene = anno_text(which = "row", gene_labels_1, location = 0.5, just = "left", gp = gpar(
					fontsize = gene_font_size, col = gene_font_col
				)), gap = unit(5, "mm")
			)
			if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
			row_ha_left <- HeatmapAnnotation(
				which = "row", simple_anno_size = unit(10,"mm"),
				gene = anno_text(which = "row", gene_labels_2, location = 0.5, just = "right", gp = gpar(
					fontsize = gene_font_size, col = gene_font_col
				)), gap = unit(5, "mm")
			)
		}
	} else {
		row_ha_left <- HeatmapAnnotation(
			which = "row", empty = anno_empty(which = "row", border = FALSE)
		)
		row_ha_right <- row_ha_left
	}
	
	message("Expreession colors...")
	col_fun <- colorRampPalette(c("white","white","orange","red","purple","black"))
	shades <- col_fun(1000)
	
	# mc labels
	message("Metacell labels...")
	collabs <- colnames(mat1)
	if (!show_mc_names) collabs <- rep("",length(collabs))
	column_lab_ha <- HeatmapAnnotation(
		which = "column",
		LAB = anno_text(which = "column", collabs, gp = gpar(fontsize = mc_font_size, rot=90))
	)
	top_column_ha <- c(column_lab_ha)
	bottom_column_ha <- c(column_lab_ha)
	empty_ha <- HeatmapAnnotation(
		empty = anno_empty(which = "column", border = FALSE, height = unit(5,"mm"))
	)
	
	# cluster colours
	if (!is.null(clust_col)){
		message("Columns...")
		if (is.null(names(clust_col))) {
			names(clust_col) <- clust_ord
		} else {
			if (!all(names(clust_col) %in% clust_ord))
				stop("Colour and cluster names do not match!")
			clust_col <- clust_col[clust_ord]
		}
		column_col_ha <- HeatmapAnnotation(
			which = "column",
			MC = colnames(mat1), col = list(MC = clust_col),
			border = c(TRUE), simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
		)
		top_column_ha <- c(column_col_ha,empty_ha,top_column_ha)
		bottom_column_ha <- c(bottom_column_ha,empty_ha,column_col_ha)
		
	}
	
	# barplots
	if (!is.null(clust_bars)) {
		
		column_bar_ha <- HeatmapAnnotation(
			which = "column",
			BAR = anno_barplot(
				anno_bar, height = 3 * clust_anno_size, bar_width = 0.9,
				gp = gpar(fill = clust_col, fontsize = mc_font_size),
				axis_param = list(gp = gpar(fontsize = mc_font_size))
			),
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
		)
		top_column_ha <- c(column_bar_ha,empty_ha,top_column_ha)
		
	}
	
	# expression heatmap
	h1 <- Heatmap(
		mat1, name = "expression", col = shades, use_raster = TRUE,
		cluster_rows = FALSE, cluster_columns = FALSE,
		show_column_names = FALSE, show_row_names = FALSE,
		#row_names_side = "left", row_labels = gene_labels_1, row_names_gp = gpar(fontsize = gene_font_size),
		right_annotation = row_ha_right, left_annotation = row_ha_left,
		column_names_gp = gpar(fontsize = mc_font_size),
		top_annotation = top_column_ha, bottom_annotation = bottom_column_ha,
		border = TRUE, show_heatmap_legend = show_heatmap_legend
	)
	
	hlist <- h1
	
	if (!is.null(clust_col) & show_clust_borders) {
		mat2 <- rbind(clust_col[match(clust_ord, names(clust_col))])
	} else {
		mat2 <- rbind(clust_ord)
	}
	
	# save figure
	if (!is.null(output_file)) {
		output_file <- stringr::str_remove(output_file,"\\.png$")
		output_file <- stringr::str_remove(output_file,"\\.pdf$")
		output_file_rds <- sprintf("%s.RDS",output_file)
		output_file_png <- sprintf("%s.png",output_file)
		output_file_pdf <- sprintf("%s.pdf",output_file)
		saveRDS(hlist, output_file_rds)
		
		# # ComplexHeatmap version
		if (plot_format == "png") {
			png(output_file_png, h=height, w=width, res=res)
		} else if (plot_format == "pdf") {
			pdf(output_file_pdf, h=height, w=width, useDingbats=TRUE)
		}
		draw(hlist, padding=unit(c(10, 100, 10, 50), "mm"), adjust_annotation_extension = FALSE)
		if (show_clust_borders) {
			change_clust <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body("expression", {
				for (i in change_clust) {
					grid.lines(x = i / ncol(mat2), y = c(0,1), gp = gpar(lty = 1, lwd = 0.25))
				}
			})
		}
		dev.off()
		
	}
	
	######################### PLOT SINGLE-CELL PROFILE ##########################
	if (plot_sc == TRUE) {
		message("Plotting single cell expression")
		cell_order=c()
		for (niche in clust_ord){
			cells=names(mc_object@mc[which(mc_object@mc == niche)])
			cell_order=c(cell_order,cells)
		}
		cluster_cell_count=as.matrix(table(mc_object@mc))
		n_cells_cluster=cluster_cell_count[clust_ord,1]
		cells_clusts <- unlist(mapply(rep, clust_ord, n_cells_cluster, USE.NAMES=FALSE))
		
		umis=as.matrix(mat_object@mat[,names(mc_object@mc)])
		mat = umis[genes, cell_order]
		totu = colSums(umis[, cell_order])
		mat = t(t(mat) / totu) * 800
		
		lus_1 = log2(1 + 7 * mat[genes[gene_ord], cell_order])
		lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
		lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x, smoothen, fill=0)))
		
		# heatmap
		mat1sc <- pmin(lus_smoo,max_expression_fc_sc)
		colnames(mat1sc) <- colnames(lus_smoo)
		shades <- colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
		col_fun <- colorRampPalette(c("white","white","orange","red","purple","black"))
		shades <- col_fun(1000)
		
		# mc labels
		top_column_ha <- HeatmapAnnotation(
			mclabstop = anno_empty(which = "column", border = FALSE, height = unit(15,"mm"))
		)
		bottom_column_ha <- HeatmapAnnotation(
			mclabsbottom = anno_empty(which = "column", border = FALSE, height = unit(15,"mm"))
		)
		# top_column_ha <- c(sc_mc_labels_ha)
		# bottom_column_ha <- c(sc_mc_labels_ha)
		
		if (!is.null(clust_col)){
			
			column_col_ha <- HeatmapAnnotation(
				which = "column",
				MC = cells_clusts, col = list(MC = clust_col),
				border = c(TRUE), simple_anno_size = 2 * clust_anno_size,
				show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
			)
			top_column_ha <- c(column_col_ha,empty_ha,top_column_ha)
			bottom_column_ha <- c(bottom_column_ha,empty_ha,column_col_ha)
			
		}
		
		# expression heatmap
		h1sc <- Heatmap(
			mat1sc, name = "sc_expression", col = shades, use_raster = TRUE,
			cluster_rows = FALSE, cluster_columns = FALSE,
			show_column_names = FALSE, show_row_names = FALSE,
			#row_names_side = "left", row_labels = gene_labels_1, row_names_gp = gpar(fontsize = gene_font_size),      show_column_names = FALSE,
			right_annotation = row_ha_right, left_annotation = row_ha_left,
			top_annotation = top_column_ha, bottom_annotation = bottom_column_ha,
			show_heatmap_legend = FALSE, border = TRUE
		)
		hlistsc <- h1sc
		
		# save figure
		if (!is.null(output_file)) {
			output_file <- stringr::str_remove(output_file,"\\.png$")
			output_file_sc_rds <- sprintf("%s_sc.RDS",output_file)
			output_file_sc_png <- sprintf("%s_sc.png",output_file)
			output_file_sc_pdf <- sprintf("%s_sc.pdf",output_file)
			saveRDS(hlistsc, output_file_sc_rds)
			if (plot_format == "png") {
				png(output_file_sc_png, h=height * plot_sc_height_cex, w=width * plot_sc_width_cex, res=res)
			} else if (plot_format == "pdf") {
				pdf(output_file_sc_pdf, h=height * plot_sc_height_cex, w=width * plot_sc_width_cex, useDingbats=TRUE)
			}
			# heatmap
			draw(hlistsc, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
			
			# add grid
			if (!is.null(clust_col)) {
				mat2sc <- rbind(clust_col[match(cells_clusts, names(clust_col))])
			} else {
				mat2sc <- rbind(cells_clusts)
				if (is.null(colnames(mat2sc)))
					colnames(mat2sc) <- cells_clusts
			}
			change_clust_sc <- which(sapply(2:ncol(mat2sc), function(i) mat2sc[,i] != mat2sc[, i - 1]))
			change_mc_sc <- c(
				which(sapply(2:ncol(mat2sc), function(i) colnames(mat2sc)[i] != colnames(mat2sc)[ i - 1])),
				ncol(mat2sc)
			)
			if (show_clust_borders) {
				decorate_heatmap_body("sc_expression", {
					for (i in change_clust_sc) {
						grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 1))
					}
				})
			}
			.add_mc_labels <- function(pos,labs) {
				for (j in 1:length(pos)) {
					i <- pos[j]
					iprev <- ifelse(j == 1,0,pos[j - 1])
					nt <- ncol(mat2sc)
					grid.text(label = labs[j], x = i / nt - (i / nt - iprev / nt) / 2, y = 0.5, gp = gpar(fontsize = mc_label_cex * mc_font_size), rot=90)
				}
			}
			if (show_mc_names) {
				decorate_annotation("mclabstop", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
				decorate_annotation("mclabsbottom", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
			}
			dev.off()
			ht_opt(RESET = TRUE)
		}
		
	}
	########################## RETURN METACELL PROFILE ##########################
	if (print_heatmap == TRUE) {
		draw(hlist)
	}
	
	message("Done.")
}