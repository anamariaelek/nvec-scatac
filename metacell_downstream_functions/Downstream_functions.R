# libraries
require("metacell")
require("tgconfig")
require("tgstat")
require("tglkmeans")
require("zoo")
require("data.table")
require("dplyr")
require("stringr")
require("ggplot2")
require("scales")
require("RColorBrewer")
require("ComplexHeatmap")
require("rasterpdf")
require("vioplot")
require("MASS")

# deactivate complexheatmap warnings
ht_opt$message = FALSE

# from  github/tanaylab/metacell/R/utils.r
# rowFunction is a function that works on rows, like rowMeans
# # much faster than tapply
.row_stats_by_factor = function (data, fact, rowFunction = rowMeans) {
	u = as.character(sort(unique(fact)))
	fact[is.na(fact)] = FALSE
	n=length(u)
	centers = matrix(NA,dim(data)[1], n, dimnames = list(rownames(data), u))
	# much faster than tapply
	for (i in u) {
		if (sum(fact == i, na.rm=TRUE) > 1) {
			centers[,i] = rowFunction(data[,fact == i,drop=FALSE])
		} else {
			centers[,i] = data[,fact == i]
		}
	}
	return(centers)
}

# Helper function to save or return plots
#' 
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param EXPR expression that produces the plot (if multiple lines, enclose it in `{ ... }`)
#' 
#' @return plot (if `output_file` is `NULL`), otherwise `NULL`
#' 
plotting_function <- function(output_file=NULL, width, height, res=NA, EXP) {
	if (!is.null(output_file)) {
		extension <- stringr::str_extract(output_file,"(png|pdf)$")
		if ( is.na(extension) ) {
			extension = "pdf"
		}
		# open graphics device
		if (extension == "png") {
			png(output_file, height = height, width = width, res=res)
		} else if (extension == "pdf" ) {
			pdf(output_file, height = height, width = width, useDingbats=TRUE)
		}
	}
	EXP
	if (!is.null(output_file)) dev.off()
}


#' Gene UMI count in metacell
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param T_totumi integer, filter out genes with less than `T_totumi` counts
#' @param grouping_vector a vector with grouping categories for individual cells
#'    (same format and length as `mc@mc`). NULL by default, which means that `mc@mc`
#'    is used as a grouping factor.
#'
#' @return matrix of gene UMI counts in metacells
#'
sca_mc_gene_counts = function(mc_object, mat_object, T_totumi=0, grouping_vector = NULL) {
	
	# get UMIs
	scr_umis = mat_object@mat
	keep_genes_bool = Matrix::rowSums(scr_umis) > T_totumi
	keep_genes = names(which(keep_genes_bool))
	
	# counts per metacell or grouping vector
	if (is.null(grouping_vector)) {
		grouping_vector = mc_object@mc
		grouping_vector_order = colnames(mc_object@mc_fp)
	} else if ("factor" %in% class(grouping_vector)) {
		grouping_vector_order = levels(grouping_vector)
	} else {
		grouping_vector_order = unique(grouping_vector)
	}
	niche_counts = .row_stats_by_factor( scr_umis[ keep_genes, names(mc_object@mc) ], grouping_vector, Matrix::rowSums )
	niche_counts = niche_counts[ , grouping_vector_order ]
	return(niche_counts)
	
}


#' Gene UMI count in metacell (object-independent)
#'
#' @param mat_object single expression matrix (e.g. mat@mat)
#' @param T_totumi integer, filter out genes with less than `T_totumi` counts
#' @param grouping_vector a vector with grouping categories for individual cells (same format and length as `mc@mc`).
#' @param cells_vector vector of cell names i, to include in counts (defaults to `colnames(mat)`; same length as `grouping_vector`).
#'
#' @return matrix of gene UMI counts in metacells
#'
sca_mc_gene_counts_noobj = function(mat, grouping_vector, T_totumi = 0, cells_vector = colnames(mat)) {
	
	# get UMIs
	keep_genes_bool = rowSums(mat) > T_totumi
	keep_genes = names(which(keep_genes_bool))
	
	# counts per metacell or grouping vector
	if ("factor" %in% class(grouping_vector)) {
		grouping_vector_order = levels(grouping_vector)
	} else {
		grouping_vector_order = unique(grouping_vector)
	}
	niche_counts = .row_stats_by_factor( mat[ keep_genes, cells_vector ], grouping_vector, rowSums )
	niche_counts = niche_counts[ , grouping_vector_order ]
	return(niche_counts)
	
}


#' Gene UMI fraction in metacell
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mc_counts matrix of gene UMI counts in metacells (output of `sca_mc_gene_counts()`)
#'
#' @return matrix of gene UMI fractions in metacells
#'
sca_mc_gene_umifrac = function(mc_object, mc_counts, multiplying_factor = 1000){
	niche_counts = mc_counts
	niche_totals = Matrix::colSums(niche_counts)
	niche_umifrac = t(apply(niche_counts,1,function(x) x * multiplying_factor / niche_totals))
	niche_umifrac = niche_umifrac[,colnames(mc_object@mc_fp)]
	return(niche_umifrac)
}


#' Gene UMI fraction in metacell (object-independent)
#'
#' @param mc_counts matrix of gene UMI counts in metacells (output of `sca_mc_gene_counts()`)
#' @param mc_columns vector indicating order of metacells (defaults to column names in `mc_counts`, another common value could be `colnames(mc@mc_fp)`)
#'
#' @return matrix of gene UMI fractions in metacells
#'
sca_mc_gene_umifrac_noobj = function(mc_counts, mc_columns = colnames(mc_counts), multiplying_factor = 1000){
	niche_counts = mc_counts
	niche_totals = colSums(niche_counts)
	niche_umifrac = t(apply(niche_counts,1,function(x) x * multiplying_factor / niche_totals))
	niche_umifrac = niche_umifrac[,mc_columns]
	return(niche_umifrac)
}


#' Compute cell type footprint based on a standard metacell->cell type
#' definition table, using same strategy as mc_fp.
#'
#' @param input_table data.frame with three columns: metacell, cell_type and color,
#'    or a path to tsv file with this table
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#'
#' @return analogous mc_object for cell types
#'
sca_cell_type_fp <- function( input_table, mc_object, mat_object, nbins=10L) {
	
	# load input table
	if (any(class(input_table) %in% "character")) {
		
		cell_type_table = read.table(input_table, header = TRUE, sep="\t", comment.char="")
		colnames(cell_type_table) = c("metacell", "cell_type", "color")
		cell_types_ordered = unique(cell_type_table$cell_type)
		rownames(cell_type_table) = cell_type_table$metacell
		
	} else if (any(class(input_table) %in% "data.frame")) {
		
		cell_type_table = input_table
		class(cell_type_table) <- "data.frame"
		colnames(cell_type_table) = c("metacell", "cell_type", "color")
		cell_types_ordered = unique(cell_type_table$cell_type)
		rownames(cell_type_table) = cell_type_table$metacell
		
	}
	
	sc_ct_label = as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
	names(sc_ct_label) = names(mc_object@mc)
	
	cells_cols = cell_type_table[as.character(mc_object@mc),"color"]
	cells_cols = as.character(cells_cols)
	names(cells_cols) = names(mc_object@mc)
	
	# filter low expression genes not included in the mc_fp
	umis = mat_object@mat
	umis = umis [ rownames(mc_object@mc_fp) , ]
	umis = umis [ , names(mc_object@mc) ]
	
	#ct_counts=t(apply(umis,1,function(x) tapply(x, sc_ct_label, sum)))
	#ct_size=colSums(ct_counts)
	#ct_umifrac=t(apply(ct_counts,1,function(x) x*1000/ct_size))
	#ct_umifrac_n=(0.1 + ct_umifrac)/apply(0.1 + ct_umifrac, 1, median)
	
	ct_geomean = tryCatch(
	  t(apply(umis, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1))),
	  error = function(e) {
	    warning(e)
	    message("Calculating geom mean for genes in ",nbins," bins")
	    umis_list <- vector("list",nbins)
	    sl <- split(1:nrow(umis),cut(1:nrow(umis),nbins))
	    for (i in seq_along(sl)) {
	      message(i," / ", nbins)
	      umisub <- umis[sl[[i]],]
	      ct_geomean_sub <- t(apply(umisub, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1)))
	      umis_list[[i]] <- ct_geomean_sub 
	    }
	    do.call(rbind, umis_list)
	  }
	)
	ct_geomean = ct_geomean[ , cell_types_ordered ]
	ct_meansize = tapply(Matrix::colSums(umis), sc_ct_label, mean)
	ideal_cell_size = pmin(1000,median(ct_meansize))
	g_fp = t(ideal_cell_size * t(ct_geomean) / as.vector(ct_meansize))
	fp_reg = 0.05
	g_fp_n = (fp_reg + g_fp) / apply(fp_reg + g_fp, 1, median)
	
	# for compatibility with other functions, return as a MC-like object
	ct_table=mc_object
	ct_table@mc_fp=g_fp_n
	ct_table@mc=sc_ct_label
	ct_table@colors=cells_cols
	ct_table@cell_names=names(sc_ct_label)
	return(ct_table)
	
}

#' Transfer cell annotations to metacells
#'
#' @param cell_df data.frame with two columns containing cell and cell type IDs;
#'    cell IDs should be the same as in `mc_df`
#' @param mc_df data.frame with two columns: cell and metacell IDs;
#'    cell IDs should be the same as in `cell_df`
#' @param NA_cell_type logical, should the missing cell type assignment (NA) be kept
#'    among final assigned cell types (default: FALSE)
#' @param niche_order character, metacell ids in the order in which they should be
#'    plotted (default: NULL)
#' @param color_df data.frame with at least two columns: cell type and color
#' @param named_colors vector of colors colors to use for reference cell types. Only 
#'    used in case no `color_df` data.frame is provided. Otherwise ignored.
#' @param do_plot whether to create a barplot or not
#' @param barplot_border barplot border color (default white, NA for none)
#' @param barplot_space barplot spacing (default 0.2)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res,xlim numeric, the width, height, resolution and maximum number of bars (metacells) in the plot. 
#'    Dimensions are in pixels if png, in inches if pdf.
#' @param interactive_plot logical, whether to output plotly interactive graph
#'    (default: FALSE). Note that this only works when working in an interactive 
#'    session i.e. when `interactive() == TRUE` and when `do_ggplot == TRUE`
#' @param legend_position charcter, one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"
#' @param legend_ncol integer, the number of columns in which to set the legend items (default: 3)
#' 
#' @return data.frame with metacells' cell_type annotations. Also saves the plot to disk.
#'
sca_transfer_annotations <- function(
	cell_df,  mc_df, NA_cell_type = FALSE,
	niche_order = NULL, color_df = NULL,
	named_colors = c("burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan", "darkgray"),
	barplot_border = "white", barplot_space = 0.2,
	do_plot = TRUE,
	output_file = NULL, width = 20, height = 4, res = NA, xlim = 200,
	interactive_plot = FALSE, 
	legend_position = "top", legend_ncol = 3
) {
	
	require(RColorBrewer)
	require(ggplot2)
	require(data.table)
	
	# input data frames as data.table
	dt <- as.data.table(copy(cell_df)[,1:2])
	mc_dt <- as.data.table(copy(mc_df)[,1:2])
	setnames(dt, c("cell", "cell_type"))
	setnames(mc_dt, c("cell", "mc"))
	dt[,cell:=as.character(cell)]
	mc_dt[,cell:=as.character(cell)]
	
	# add mc annotation to cells
	dt[mc_dt, on="cell", mc:=i.mc]
	# remove cells that are not assigned to any mc
	missing_cells <- is.na(dt$mc)
	dt <- dt[!missing_cells]
	message(sprintf(
		"Removed %s cells (out of %s) that are not in any metacell!",
		sum(missing_cells),length(missing_cells)
	))
	# for each mc, count how many cells in it are assigned to each ct
	dt[,cells_from_mc_in_ct:=.N,by=.(mc,cell_type)]
	
	# find most occuring ct for each mc
	if (NA_cell_type) {
		dt[,max_cells_ct:=max(cells_from_mc_in_ct),mc]
		dt[,max := { cells_from_mc_in_ct == max_cells_ct } ]
		max_ct <- unique(dt[max == TRUE,.(cell_type,mc)])
	} else {
		dtfilt <- dt[!is.na(cell_type)][,max_cells_ct:=max(cells_from_mc_in_ct),mc]
		dtfilt[,max := { cells_from_mc_in_ct == max_cells_ct } ]
		max_ct <- unique(dtfilt[max == TRUE,.(cell_type,mc)])
	}
	dt[max_ct, on="mc", assigned_cell_type:=i.cell_type]
	setorder(dt, assigned_cell_type, mc, na.last = TRUE)
	
	# order metacells
	if (is.null(niche_order)) {
		# if there's not predetermined order, convert metacells into an ordered factor
		mc_levels = sort(as.numeric(unique(dt$mc)))
		dt[,mc:=factor(mc, levels = mc_levels)]
	} else {
		# if there's a predetermined order (niche_order), use this
		dt[,mc:=factor(mc, levels = niche_order)]
	}
	
	# create vector of named colors mapped to the original cell types
	if (!is.null(color_df)) {
		col_dt <- as.data.table(copy(color_df)[,1:2])
		setnames(col_dt, c("cell_type","color"))
		cols <- col_dt$color
		names(cols) <- col_dt$cell_type
	} else {
		cpalette <- colorRampPalette(named_colors)
		cts <- sort(unique(dt$cell_type))
		cols <- cpalette(length(cts))
		names(cols) <- cts
		col_dt <- data.table("cell_type"=names(cols),"color"=cols)
	}
	
	# create factor for original cell types
	dt$cell_type = factor(dt$cell_type, levels=sort(unique(dt$cell_type)))
	dt$cell_type = factor(dt$cell_type, levels=sort(unique(dt$cell_type)))
	
	# plot
	do_ggplot = interactive_plot
	if (do_plot | do_ggplot) {
		plotting_function(output_file, width, height, res, EXP = {
			# plot: either ggplot or base R
			if ( do_ggplot ) {
				
				# plotly doesn't like factors, but this means the ordering of bars will change...
				if (interactive_plot) {
					dt[,cell_type:=as.character(cell_type)]
					dt[,mc:=as.character(mc)]
				}
				
				gp <- ggplot(dt, aes(x=mc, fill=cell_type), color="black") +
					geom_bar(position = "fill") +
					scale_x_discrete(breaks=dt$mc, labels=dt$mc) +
					scale_y_continuous(expand = c(0,0)) +
					scale_fill_manual(values = cols) +
					theme_minimal() + theme(panel.grid = element_blank()) +
					labs(x="metacells",y="percent of cells in mc") +
					theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
				
				if (interactive() & interactive_plot) {
					require(plotly)
					ggplotly(gp)
				}
				
				# close ggplot
				print(gp)
				
			} else {
				
				dt_crosstab = table(dt$cell_type, dt$mc)
				dt_crosstab_frac = apply(dt_crosstab, 2 , function(x) x / sum(x) )
				# plot fraction
				barplot(
					dt_crosstab_frac,
					col=cols[rownames(dt_crosstab)],
					border = barplot_border,
					space = barplot_space,
					xlim = NULL,
					las = 2,
					cex.names = 0.8,
					xlab = "metacells",
					ylab = "Frequency of annotations"
				)
				#legend(legend_position, legend = rownames(dt_crosstab), fill = cols[rownames(dt_crosstab)], cex=0.6, bty="n", border="white", ncol=3)
				# plot totals
				barplot(
					dt_crosstab,
					col=cols[rownames(dt_crosstab)],
					border = barplot_border,
					space = barplot_space,
					xlim = NULL,
					las = 2,
					cex.names = 0.8,
					xlab = "metacells",
					ylab = "Frequency of annotations"
				)
				#legend(legend_position, legend = rownames(dt_crosstab), fill = cols[rownames(dt_crosstab)], cex=0.6, bty="n", border="white", ncol=3)
				# add legend
				plot(NA,NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
				legend(legend_position, legend = rownames(dt_crosstab), fill = cols[rownames(dt_crosstab)], cex=0.6, bty="n", border = barplot_border, ncol=legend_ncol)
				
			}
			
		})
	}
	
	# return reannotated dataframe
	col_dt2 <- copy(col_dt)
	setnames(col_dt2,"cell_type","assigned_cell_type")
	dt_col <- dt[col_dt, on="cell_type", color:=i.color]
	out_dt <- dt_col[col_dt2, on="assigned_cell_type", assigned_color:=i.color]
	out_df <- copy(out_dt)
	class(out_df) <- "data.frame"
	return(out_df)
	
}


#' Record frequency of cell type annotations in metacells
#'
#' @param mc_vector vector of metacell assignments for each single cell, e.g. the `mc` column in 
#'    the output of `sca_transfer_annotations`.
#' @param assigned_annot_vector vector of cell type annotations for each cell, e.g. `cell_type` 
#'    column in the output of `sca_transfer_annotations`.
#' @param freq_threshold float; assign an annotation to the metacell if it is found at >= 
#'    this frequency (default: 0.667)
#' @param freq_min float; ignore annotations below this frequency (default: 0.333)
#' @param collapse logical; collapse annotations found above the `freq_min` value into a single
#'    vector (default is TRUE). If FALSE, most common annotation above `freq_min` is returned.
#' @param empty_value character; value to report if no annotation fulfills any criteria 
#'    (default: "None")
#' @return data.frame with metacells' cell_type annotations.
#'
sca_metacell_annotation_freq = function(
	mc_vector, assigned_annot_vector, 
	freq_threshold = 2 / 3, 
	freq_min = 1 / 3, 
	collapse = TRUE,
	empty_value = "None") {
	
	# cross-tabulate per-cell mc vector and per-cell assigned annotation vector
	xtab = table(mc_vector, as.character(assigned_annot_vector))
	xtaf = xtab / rowSums(xtab)
	xtaf[is.nan(xtaf) | is.na(xtaf)] <- 0
	
	# find most common assignments per metacell at a dertain frequency threshold
	# if more than one assignment exists at that frequency threshold, concatenate
	# assignments
	assignments = apply(xtaf, 1, function(x) {
		ix = which(x > freq_threshold)
		ix_max = which.max(x)
		ix_minor = which(x > freq_min)
		
		# find most common annotation
		a = colnames(xtab)[ix_max]
		f = sprintf("%.3f", x[ix_max])
		
		# if there is more than one annotation above threshold and we have to collapse:
		if (length(ix) == 0 & length(ix_minor) >= 1 & collapse) {
			a = paste(colnames(xtab)[ix_minor], collapse = " | ")
			f = paste(sprintf("%.3f",x[ix_minor]), collapse = " | " )
		}
		
		if (length(ix) == 0 & length(ix_minor) == 0) {
			a = empty_value
			f = sprintf("%.3f", x[ix_max])
		}
		
		# return list
		return(list(a, f))
	})
	
	# list to matrix
	assignments_mat = matrix(unlist(assignments), nrow=2)
	
	# return annotations
	anns = data.frame(
		mc = rownames(xtab),
		assigned_annot = assignments_mat[1,],
		annot_freqs = assignments_mat[2,]
	)
	
	return(anns)
	
}


#' Split metacells by annotation
#' 
#' Split metacells that have cells assigned to more than one cell type by individual cell annotations.
#' Only works for splitting every selected metacell in two (TBA: allow multiple splits)
#' 
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param new_mc_id character, name of the new metacell object to add to database;
#'    if NULL (default), the object is not added to metacell database
#' @param cell_df data.frame with two columns containing cell and cell type IDs
#' @param anns_df data.frame with at least two columns containing metacell IDs and cell type annotations, 
#'    e.g. the output of `sca_metacell_annotation_freq`; 
#'    cell type annotations should be the same as in `cell_df`, and combination of those separated by " | " 
#' @return metacell object (`gMCCov` class) with metacells that were assigned multiple annotations split
#'
sca_metacell_annotation_split <- function(
	mc_object, mat_object, 
	cell_df, anns_df,
	new_mc_id = NULL
) {
	# cell-metacell asignments from mc object
	mc_df <- data.table(
		cell = names(mc@mc),
		metacell = mc@mc
	)
	# mc vector to be updated
	new_mc <- mc@mc
	# assigned annotations
	anns_df <- anns_df[,1:2]
	setDT(anns_df)
	setnames(anns_df,c("mc","assigned_annot"))
	split_mcs <- anns_df[grep("\\|", assigned_annot)]$mc
	mcss <- sort(as.integer(as.character(split_mcs)))
	for (i in mcss) {
		message(sprintf("Splitting metacell %s into %s",i,2))
		# increase mc ids after the one to be split
		inc <- match(i,mcss) - 1
		mci <- i + inc
		new_mc[new_mc > mci] <- new_mc[new_mc > mci] + 1
		# get cell ids to split
		cts <- str_split(anns_df[mc == i]$assigned_annot, " \\| ")[[1]]
		cell_ids <- sapply(cts, function(ct)
			merge.data.table(mc_df[metacell == i], cell_df[cell_type == ct])$cell, 
			USE.NAMES = TRUE, simplify = FALSE
		)
		# assign new mc ids to cells
		new_mc[cell_ids[[1]]] <- mci
		new_mc[cell_ids[[2]]] <- mci + 1
	}
	new_mc <- tgMCCov(new_mc, mc@outliers, mat)
	if (!is.null(new_mc_id)) scdb_add_mc(new_mc_id, new_mc)
	return(new_mc)
}

#' Plot barplot of metacell sizes
#'
#' @param mc_object
#' @param mat_object
#' @param mc_annotation data.frame, metacell annotation; or character, path to metacell annotation file
#' @param mc_color vector of metacell colors, not used if `mc_annotation` is supplied
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#'
scp_metacell_size <- function(
    mc_object, 
    mat_object, 
    mc_annotation = NULL,
    mc_color = NULL,
    output_file = NULL, 
    width = 24, height = 12, res = NA
) {
  
  mc_cs <- data.table(table(mc_object@mc))
  setnames(mc_cs, c("metacell","N"))
  mc_cs[,metacell:=factor(metacell, levels=seq_along(unique(mc_cs$metacell)))]
  
  # add annotations
  if (!is.null(mc_annotation)) {
    if ("character" %in% class(mc_annotation)) {
      mc_ct <- fread(mc_annotation)
    } else if ("data.frame" %in% class(mc_annotation)) {
      mc_ct <- as.data.table(mc_annotation)
    }
    mc_ct[,metacell:=factor(metacell,levels=levels(mc_cs$metacell))]
    mc_ct[,cell_type:=factor(cell_type,levels=unique(mc_ct$cell_type))]
    mc_cs <- merge.data.table(mc_cs, mc_ct, by="metacell", sort=FALSE)
    ct_cols <- unique(mc_ct[,.(cell_type,color)])
    ct_cols <- structure(ct_cols$color, names=as.character(ct_cols$cell_type))
    cttab <- table(names(
      structure(mc_ct$metacell, names=as.character(mc_ct$cell_type))[mc_object@mc]
    ))
    cttab <- cttab[levels(mc_cs$cell_type)]
    gp <- ggplot(mc_cs, aes(metacell,N,label=N,fill=cell_type)) + 
      geom_bar(stat = "identity") +
      scale_fill_manual(
        values = ct_cols,
        limits = names(cttab),
        labels = sprintf("%s (%s)",names(cttab), format(cttab, big.mark=",", trim=TRUE))
      ) +
      theme(legend.position = "bottom", legend.title = element_blank())
  } else {
    if (is.null(mc_color))
      mc_color <- "grey"
    gp <- ggplot(mc_cs, aes(metacell,N,label=N)) + 
      geom_bar(stat = "identity", fill = mc_color)
  }
  
  # add counts
  gp <- gp + 
    geom_text(angle = 90, vjust = 0.5, hjust = 0, size = 2.6) +
    scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    labs(x = "metacells", y = "cells")
  
  
  # save plot
  plotting_function(output_file, width, height, res, EXP = {print(gp)})
  
}

#' Plot barplot of metacell batch distribution
#'
#' Plots barplot of batch distribution of cells in each metacell
#'
#' @param mc_object
#' @param mat_object
#' @param clust_order ordered metacells (all values should be in `mc_object@mc`)
#' @param batch_field character, name of the column in `mat_object@cell_metadata` containing batch information
#' @param color_field character, name of the column in `mat_object@cell_metadata` containing color information
#' @param update_color_file NULL
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#'
scp_batch_distribution <- function(
  mc_object, mat_object, clust_order=NULL,
	batch_field="dataset", color_field="color", update_color_file=NULL,
	output_file=NULL, width=1200, height=600, res=NA
){
	if (is.null(clust_order)) {
		clust_order=sort(unique(mc_object@mc))
	}

	x=sapply(clust_order, function(x) 
	  table(
	    factor(
	      mat_object@cell_metadata[names(which(mc_object@mc == x)),batch_field],
	      levels=unique(mat_object@cell_metadata[,batch_field])
	    )
	  )
	)
	colnames(x)=as.character(clust_order)
	# get batches
	categories=as.character(unique(mat_object@cell_metadata[,batch_field]))
	# check if there's a defined batch color
	if (is.null(mat_object@cell_metadata$color)){
		n_colors=length(categories)
		color=colorRampPalette(brewer.pal(8, "Dark2"))(n_colors)
		names(color)=categories
	} else {
		color=sort(sapply(categories,function(x) 
		  as.character(
		    unique(
		      mat_object@cell_metadata[which(mat_object@cell_metadata[,batch_field] == x),color_field]
		    )
		  )
		))
		color=color[rownames(x)]
	}
	if (!is.null(update_color_file)) {
		update_color_table=read.table(update_color_file,h=TRUE,row.names=1,sep="\t",comment.char="")
		update_color=as.character(unique(update_color_table$color))
		names(update_color)=unique(update_color_table[,batch_field])
		color=update_color[names(color)]
	}
	# batch frequencies
	x_colSums=colSums(x)
	x_freq=t(apply(x,1,function(y) y / x_colSums))
	bckg_freq=table(mat_object@cell_metadata[,batch_field]) / ncol(mat_object@mat)
	FC=apply(x_freq,2,function(y) y / bckg_freq)
	logFC=log2(FC + 1)
	
	# save
	if (is.null(.scfigs_base) | !file.exists(.scfigs_base)) {
		stop("figs directory at ", .scfigs_base, " is missing")
	}
	output_file_base <- sprintf("%s/%s.%s", .scfigs_base, mc_id, "batch_distribution")
	
	plotting_function(output_file, width, height, res, EXP = {
		barplot(x,las=2,col=color,cex.names=1,cex.axis=2)
		legend("topleft",legend=names(color),fill=color,box.lty=0,horiz=FALSE,ncol=round(length(color) / 3,0),cex=1.5)
	})
	
	return(x)
	
}

#' Plot barplot of metacell batch distribution
#'
#' @param mc_object
#' @param mat_object
#' @param batch_field character, name of the column in `mat_object@cell_metadata` containing batch information, or a custom vector (must match cell order in `mat`!)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#'
scp_batch_effects_per_mc = function(
    mc_object, 
    mat_object, 
    batch_field = "batch_set_id", 
    mc_color = NULL, 
    batch_color = NULL,
    output_file = NULL, 
    width = 24, height = 12, res = NA
) {
	
	# crosstabulate
	if (length(batch_field) == 1) {
		xtab = t(table(mc_object@mc, mat_object@cell_metadata [ names(mc_object@mc) , batch_field]))
	} else if (length(batch_field) == length(mc_object@mc)) {
		xtab = t(table(mc_object@mc, batch_field))
	} else {
		stop("`batch_field` should either be a column in `mat@cell_metadata` or a vector with cell-level batch information, matching `mc@mc`")
	}
	# keep metacell order from mc object
	xtab = xtab [ , colnames(mc_object@mc_fp) ]
	
	# colors for each batch
	if (is.null(batch_color)) {
	  color_palette = colorRampPalette(c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet"))
	  batch_color = color_palette(n=nrow(xtab))
	}
	
	# start plotting device
	plotting_function(output_file, width, height, res, EXP = {
		
		par(mfrow=c(2,1))
		
		# counts plot
		barplot(
			xtab, 
			las = 2,
			xlim = c( 1, ncol(mc_object@mc_fp) + 60 ),
			ylim = c(-5, max(colSums(xtab))),
			ylab = "# cells", xlab = "Metacells",
			border = "white",
			space = 0,
			col = batch_color)
		legend_t = sprintf("%s n = %i", rownames(xtab), rowSums(xtab))
		legend("topright", legend = legend_t, fill = batch_color, bty = "n", cex = 0.8, title = batch_field)
		if (!is.null(mc_color)) { 
			points(y = rep(-2, length(mc_color)), x = 1:length(mc_color) - 0.5, bg = mc_color, col = mc_color, pch = 22)
		}
		
		# fraction plot
		xtaf = t(t(xtab) / colSums(xtab))
		barplot(
			xtaf * 100, 
			las = 2,
			xlim = c( 1, ncol(mc_object@mc_fp) + 60 ),
			ylim = c(-5,100),
			ylab = "% cells", xlab = "Metacells",
			border = "white",
			space = 0,
			col = batch_color)
		legend_t = sprintf("%s n = %i", rownames(xtab), rowSums(xtab))
		legend("topright", legend = legend_t, fill = batch_color, bty = "n", cex = 0.8, title = batch_field)
		if (!is.null(mc_color)) { 
			points(y = rep(-2, length(mc_color)), x = 1:length(mc_color) - 0.5, bg = mc_color, col = mc_color, pch = 22)
		}

		# observed/expected plot per sample
		par(mfrow=c(4,1))
		expected = rowSums(xtab) / sum(xtab)
		for (batch in rownames(xtab)) {
			ob = xtaf[batch,]
			obex = ob / expected
			barplot(
				log2(obex), 
				xlim = c( 1, ncol(mc_object@mc_fp) + 60 ),
				ylim = c(-2,2),
				space = 0, las = 2, border = "white",
				main = batch, col = mc_color, 
				ylab = "log2(Observed/expected)")
			abline(h = 0, lty = 2)
		}
		
	})
	
}

#' Comparison of gene expression of two clusterig solutions
#'
#' @param matrix1,matrix2 gene fold change matrices, with genes in rows and
#'    clusters in columns; shold have both row names and column names
#' @param markers character, optional list of marker genes to use for
#'    correlation calculation
#' @param master_gene character
#' @param cor_method character, comparison metric to use, should be one of the following
#'    `c("overlap","jaccard","pearson","kendall","spearman")`
#' @param fc_thrs numeric, lfc threshold to use for binarizing gene expression in clusters
#'    (default: 2); used for calculation when `cor_method` is one of `c("overlap",jaccard")`,
#'    and for reporting overlapping genes
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param original_ordering_1st
#' @param annotation_file_1
#' @param reorder_by_ann1,reorder_by_ann2 logical, whether to plot metacells in order
#'   in which they appear in annotation files (default: FALSE)
#' @param cor_color character, colors to use for plotting
#' @param plot_type character, either `"heatmap"` (default) or `"dotplot"`
#' @param annotation_size numeric, height of the annotation color bar
#' @param label_font_size numeric, size of annotation labels
#' @param cor_max numeric, color scaling max value (default: 1)
#' @param cex_dot numeric, dot size scaling factor (default: 1),
#'    only used if `plot_type == "dotplot"`
#' @param grid logical
#' @param annotation_grid_1,annotation_grid_2 logical
#'
#' @return a list with following elements:
#'   1) `heatmap` complex heatmap object
#'   2) `cor_matrix` correlation matrix used for plotting
#'   3) `overlap_matrix` matrix with number of overlapping genes
#'   4) `overlapping_genes` list of overlapping genes, nested list where at
#'   the first level are the columns from the first matrix, and at the second
#'   level are the columns from the second matrix.
#'
sca_clustering_comparison <- function(
	matrix1, matrix2,
	markers=NULL, master_gene=NULL,
	cor_method="jaccard", fc_thrs=1.2,
	output_file=NULL, width=3000, height=3000, res=NA,
	name1=NULL, name2=NULL,
	original_ordering_1st=FALSE, original_ordering_2nd=FALSE,
	annotation_file_1=NULL, annotation_file_2=NULL,
	reorder_sps2=FALSE, reorder_by_ann1=FALSE, reorder_by_ann2=FALSE,
	cor_color = NULL, plot_type="heatmap", annotation_size = 10, label_font_size = 12, cor_max=1, cex_dot=1,
	grid = TRUE,  annotation_grid_1 = TRUE, annotation_grid_2 = TRUE
) {
	# select marker genes
	if (!is.null(markers)){
		markers=intersect(markers,rownames(matrix1))
		markers=intersect(markers,rownames(matrix2))
		mat1=matrix1[markers,]
		mat2=matrix2[markers,]
	} else {
		mat1=matrix1
		mat2=matrix2
	}
	# binarize expression
	intg <- intersect(rownames(mat1),rownames(mat2))
	mat1[mat1 < fc_thrs] <- 0
	mat2[mat2 < fc_thrs] <- 0
	mat1[!(mat1 < fc_thrs)] <- 1
	mat2[!(mat2 < fc_thrs)] <- 1
	mat1 <- data.matrix(mat1)[intg,]
	mat2 <- data.matrix(mat2)[intg,]
	# genes in common
	out_list <- vector("list",length = ncol(mat1))
	for (i in 1:ncol(mat1)) {
		intglist <- lapply(1:ncol(mat2), function(j) {
			ints <- which(mat1[,i] > 0 & mat2[,j] > 0)
			if (length(ints) == 0) {
				NA
			} else {
				rownames(mat1)[ints]
			}
		})
		names(intglist) <- colnames(mat2)
		out_list[[i]] <- intglist
	}
	names(out_list) <- colnames(mat1)
	
	# compare clusters
	if (cor_method %in% c("jaccard","overlap")) {
		cl_comp=get(cor_method)(mat1,mat2)
	} else if (cor_method %in% c("pearson","spearman","kendall")) {
		cl_comp=cor(mat1,mat2,method=cor_method)
	}
	rownames(cl_comp)=as.character(colnames(mat1))
	colnames(cl_comp)=as.character(colnames(mat2))
	
	# ordering (if order is not retained)
	if (is.null(master_gene) & !original_ordering_1st & !reorder_by_ann1){
		hc=hclust(dist(cor(t(cl_comp))),method="ward.D2")
		tmp_cor=cl_comp[hc$order,]
		rownames(tmp_cor)=seq(1,nrow(tmp_cor))
		order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
	} else if (original_ordering_1st & !original_ordering_2nd & !reorder_by_ann1){
		hc=list()
		hc$order=colnames(mat1)
		tmp_cor=cl_comp[,]
		rownames(tmp_cor)=seq(1,nrow(tmp_cor))
		order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
	} else if (!is.null(master_gene) & !reorder_by_ann1 & !reorder_by_ann2){
		order_cols=names(sort(matrix2[master_gene,],decreasing=TRUE))
		hc=hclust(dist(cor(t(cl_comp[,order_cols]))),method="ward.D2")
	} else {
		hc=list()
		hc$order=colnames(mat1)
		order_cols=colnames(mat2)
	}
	# plotting matrix
	cor_mat <- cl_comp[hc$order,order_cols]
	cor_name <- sprintf("%s \n",cor_method)
	col_ord <- colnames(cor_mat)
	row_ord <- rownames(cor_mat)
	
	# metacell annotations
	if (!is.null(annotation_file_1)){
		clust_anno_size  <- unit(annotation_size,"mm")
		annr <- read.table(annotation_file_1,header=TRUE,comment.char="",sep="\t")
		if (reorder_by_ann1 == TRUE) {
			row_ord <- as.character(annr$metacell)
			cor_mat <- cor_mat[row_ord,]
		}
		rid <- unlist(lapply(row_ord,match,table=annr$metacell))
		row_clusts <- annr[rid,]$cell_type
		row_clust_col <- annr[rid,]$color
		names(row_clust_col) <- row_clusts
		right_row_col_ha <- HeatmapAnnotation(
			which = "row",
			MC = row_clusts,
			lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
			col = list(MC = row_clust_col),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size / 3,"mm")
		)
		left_row_col_ha <- HeatmapAnnotation(
			which = "row",
			lab = anno_text(which = "row", row_ord, just = "right", location = unit(1, "npc"), gp = gpar(fontsize = label_font_size)),
			MC = row_clusts,
			col = list(MC = row_clust_col),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size / 3,"mm")
		)
	} else {
		right_row_col_ha <- HeatmapAnnotation(
			which = "row",
			lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
		)
		left_row_col_ha <- HeatmapAnnotation(
			which = "row",
			lab = anno_text(which = "row", row_ord, just = "right", location = unit(1, "npc"), gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
		)
	}
	
	if (!is.null(annotation_file_2)){
		clust_anno_size  <- unit(annotation_size,"mm")
		annc <- read.table(annotation_file_2, header=TRUE,comment.char="",sep="\t")
		if (reorder_by_ann2 == TRUE) {
			col_ord <- as.character(annc$metacell)
			cor_mat <- cor_mat[,col_ord]
		}
		cid <- unlist(lapply(col_ord,match,table=annc$metacell))
		col_clusts <- annc[cid,]$cell_type
		col_clust_col <- annc[cid,]$color
		names(col_clust_col) <- col_clusts
		top_column_col_ha <- HeatmapAnnotation(
			which = "column",
			lab = anno_text(which = "column", col_ord, just = "left", location = unit(0, "npc"), gp = gpar(fontsize = label_font_size)),
			MC = col_clusts,
			col = list(MC = col_clust_col),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size / 3,"mm")
		)
		bottom_column_col_ha <- HeatmapAnnotation(
			which = "column",
			MC = col_clusts,
			lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
			col = list(MC = col_clust_col),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size / 3,"mm")
		)
	} else {
		top_column_col_ha <- HeatmapAnnotation(
			which = "column",
			lab = anno_text(which = "column", col_ord, just = "left", location = unit(0, "npc"), gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
		)
		bottom_column_col_ha <- HeatmapAnnotation(
			which = "column",
			lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
		)
	}
	
	# intersect
	cl_int <- overlap(mat1,mat2)
	ovrl_mat <- cl_int[row_ord,col_ord]
	sf <- max(ovrl_mat,na.rm=TRUE)
	cor_mat_sc <- ovrl_mat / sf
	
	# heatmap
	if (is.null(name1)) {
		row_title=""
	} else {
		row_title=name1
	}
	
	if (is.null(name2)) {
		column_title=""
	} else {
		column_title=name2
	}
	
	if (plot_type == "dotplot") {
		if (is.null(cor_color))
			cor_color <- circlize::colorRamp2(c(0, cor_max), c("white", "red"))
		hm <- Heatmap(
			cor_mat, col = cor_color, name = cor_name, border = TRUE, rect_gp = gpar(type = "none"),
			cell_fun = function(j, i, x, y, width, height, fill) {
				if (grid == TRUE)
					grid.rect(
						x = x, y = y, width = width, height = height,
						gp = gpar(col = "gray50", fill = NA, lty = 1, lwd = 0.1)
					)
				grid.circle(
					x = x, y = y, r = abs(cor_mat_sc[i, j]) / 2 * max(unit.c(width, height)) * cex_dot,
					gp = gpar(fill = cor_color(cor_mat[i, j]), col = "white")
				)
			},
			cluster_rows = FALSE, cluster_columns = FALSE,
			show_row_names = FALSE, show_column_names = FALSE,
			row_title = row_title, column_title = column_title,
			right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,
			top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
			heatmap_legend_param = list(
				#col_fun = cor_color, at = cols_range,
				title = cor_name, border = TRUE,
				legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
				title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
				labels_gp = gpar(fontsize = label_font_size)
			)
		)
		
	} else if (plot_type == "heatmap") {
		if (is.null(cor_color))
			cor_color <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
		hm <- Heatmap(
			pmin(cor_mat,cor_max), col = cor_color, name = cor_name, border = TRUE,
			rect_gp = gpar(col = "gray50", lwd = 0.2),
			cluster_rows = FALSE, cluster_columns = FALSE,
			show_row_names = FALSE, show_column_names = FALSE,
			row_title = row_title, column_title = column_title,
			right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,
			top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
			heatmap_legend_param = list(
				title = cor_name, border = TRUE,
				legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
				title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
				labels_gp = gpar(fontsize = label_font_size)
			)
		)
	}
	
	plotting_function(output_file, width, height, res, EXP = {
		
		ht_opt(
			COLUMN_ANNO_PADDING=unit(15,"mm"),
			ROW_ANNO_PADDING=unit(15,"mm"),
			DIMNAME_PADDING=unit(5,"mm")
			#HEATMAP_LEGEND_PADDING=unit(10,"mm")
		)
		draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
		
		if (!is.null(annotation_file_1) & annotation_grid_1 == TRUE){
			mat2 <- rbind(row_clust_col)
			change_clust_row <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body(cor_name, {
				for (i in change_clust_row) {
					grid.lines(x = c(0,1), y = 1 - i / ncol(mat2),  gp = gpar(col = "gray50", lwd = 0.05))
				}
			})
		}
		if (!is.null(annotation_file_2) & annotation_grid_2 == TRUE){
			mat2 <- rbind(col_clust_col)
			change_clust_col <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body(cor_name, {
				for (i in change_clust_col) {
					grid.lines(x = i / ncol(mat2), y = c(0,1), gp = gpar(col = "gray50", lwd = 0.05))
				}
			})
		}
		
	})
	
	ht_opt(RESET=TRUE)
	
	return(list(
		heatmap=hm, comp=cor_mat, overlap_mat=ovrl_mat, overlapping_genes=out_list
	))
	
}

#' Correlation between gene expression of two clusterig solutions
#'
#' @param matrix1,matrix2 matrices, with genes in rows and clusters in columns;
#'    shold have both row names and column names
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param markers character, optional list of marker genes to use for
#'    correlation calculation
#' @param master_gene
#'
sca_clustering_correlation <- function(
	matrix1, matrix2, markers=NULL, master_gene=NULL,
	pmin=1, pmax=0, cor_method="pearson",
	original_ordering_1st=FALSE, original_ordering_2nd=FALSE,
	annotation_file_1=NULL, annotation_file_2=NULL,
	output_file=NULL, width=3000, height=3000, res=NA,
	annotation_size = 10,  label_font_size = 12
	# cex_row=1, cex_column=1, x_lab_rot=FALSE
){
	
	if (!is.null(markers)){
		markers=intersect(markers,rownames(matrix1))
		markers=intersect(markers,rownames(matrix2))
		mat1=matrix1[markers,]
		mat2=matrix2[markers,]
	}
	if (is.null(markers)){
		mat1=matrix1
		mat2=matrix2
	}
	
	#colnames(scr_markers_km$centers)
	cl_cor=cor(mat1,mat2,method=cor_method)
	rownames(cl_cor)=as.character(colnames(mat1))
	colnames(cl_cor)=as.character(colnames(mat2))
	
	
	
	if (is.null(master_gene) & !original_ordering_1st){
		hc=hclust(dist(cor(t(cl_cor))),method="ward.D2")
		tmp_cor=cl_cor[hc$order,]
		rownames(tmp_cor)=seq(1,nrow(tmp_cor))
		order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
	} else if (original_ordering_1st & !original_ordering_2nd){
		hc=list()
		hc$order=colnames(mat1)
		tmp_cor=cl_cor[,]
		rownames(tmp_cor)=seq(1,nrow(tmp_cor))
		order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
	}else if (original_ordering_1st & original_ordering_2nd){
		hc=list()
		hc$order=colnames(mat1)
		order_cols=colnames(matrix2)
	} else if (!is.null(master_gene)){
		order_cols=names(sort(matrix2[master_gene,],decreasing=TRUE))
		hc=hclust(dist(cor(t(cl_cor[,order_cols]))),method="ward.D2")
	}
	
	cor_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
	
	ht_opt(
		COLUMN_ANNO_PADDING=unit(15,"mm"),
		ROW_ANNO_PADDING=unit(15,"mm"),
		DIMNAME_PADDING=unit(5,"mm"),
		HEATMAP_LEGEND_PADDING=unit(15,"mm")
	)
	
	cor_mat <- (pmin(pmax(cl_cor[hc$order,order_cols],pmax),pmin))
	cor_name <- sprintf("%s \n",cor_method)
	
	col_ord <- colnames(cor_mat)
	row_ord <- rownames(cor_mat)
	
	if (!is.null(annotation_file_1)){
		clust_anno_size  <- unit(annotation_size,"mm")
		annr <- read.table(annotation_file_1,header=TRUE,quote="",sep="\t",comment.char="")
		rid <- unlist(lapply(row_ord,match,table=annr$metacell))
		row_clusts <- annr[rid,]$cell_type
		row_clust_col <- annr[rid,]$color
		names(row_clust_col) <- row_clusts
		row_col_ha <- HeatmapAnnotation(
			which = "row",
			MC = row_clusts, col = list(MC = row_clust_col),
			lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size / 3,"mm")
		)
	} else {
		row_col_ha <- HeatmapAnnotation(
			which = "row",
			lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
		)
	}
	
	
	if (!is.null(annotation_file_2)){
		clust_anno_size  <- unit(annotation_size,"mm")
		annc <- read.table(annotation_file_2, header=TRUE,quote="",sep="\t",comment.char="")
		cid <- unlist(lapply(col_ord,match,table=annc$metacell))
		col_clusts <- annc[cid,]$cell_type
		col_clust_col <- annc[cid,]$color
		names(col_clust_col) <- col_clusts
		column_col_ha <- HeatmapAnnotation(
			which = "column",
			MC = col_clusts, col = list(MC = col_clust_col),
			lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size / 3,"mm")
		)
	} else {
		column_col_ha <- HeatmapAnnotation(
			which = "column",
			lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
		)
	}
	
	hm <- Heatmap(
		cor_mat, col = cor_color, name = cor_name, border = TRUE,
		rect_gp = gpar(col = "gray50", lwd = 0.2),
		cluster_rows = FALSE, cluster_columns = FALSE,
		show_row_names = FALSE, show_column_names = FALSE,
		right_annotation = row_col_ha, left_annotation = row_col_ha,
		top_annotation = column_col_ha, bottom_annotation = column_col_ha,
		heatmap_legend_param = list(
			title = cor_name, border = TRUE,
			legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
			title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
			labels_gp = gpar(fontsize = label_font_size)
		)
	)
	
	plotting_function(output_file, width, height, res, EXP = {
		
		draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
		if (!is.null(annotation_file_1)){
			mat2 <- rbind(row_clust_col)
			change_clust_row <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body(cor_name, {
				for (i in change_clust_row) {
					grid.lines(x = c(0,1), y = 1 - i / ncol(mat2), gp = gpar(col = "gray50", lwd = 0.2))
				}
			})
		}
		if (!is.null(annotation_file_2)){
			mat2 <- rbind(col_clust_col)
			change_clust_col <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body(cor_name, {
				for (i in change_clust_col) {
					grid.lines(x = i / ncol(mat2), y = c(0,1), gp = gpar(col = "gray50", lwd = 0.2))
				}
			})
		}
		
	})
	ht_opt(RESET=TRUE)
	
	return(cl_cor[hc$order,order_cols])
	
}



#' Prepare heatmap of gene expression: select marker genes & create annotations
#'
#' Prepares heatmap of gene expression fold change for metacells and single cells
#' (no plotting done, returns list with selected markers and prepared annotations),
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param black_list character, blacklisted genes
#' @param sub_list_mc
#' @param gene_list character, list of genes to plot, if NULL (default) ...
#' @param order_genes logical, whether to cluster genes (default: TRUE)
#' @param gene_annot_file charcter, file path to the file containig gene annotations,
#'    it should have three tab separated columns containing gene ID, pfam architecture
#'    and any additional annotation in the last column
#' @param annot_header logical, gene annotation file has column names?
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param clust_ord character, metacells in the order in which they should be
#'    plotted; if cluster order is not specified (default: NULL), it is
#'    determined by hierarchical clustering
#' @param per_clust_genes integer, how many genes per cluster to aim to show in the heatmap
#'    (default: 20)
#' @param gene_min_fold numeric, minimum fold change for a gene to be considered for plotting
#'    (default: 2)
#' @param transverality_N integer, number of metacells in which a gene can be highly expressed (>transverality_thr, by default >1.4)
#'    to be considered for plotting, by default this is the total number of metacells
#' @param transverality_thr integer, expression threshold for the transverality_N filter (>1.4)
#' @param transv_excluded_mc character, metacells to be excluded in transversality calculation
#'    (default: NULL)
#' @param output_file optionally, a path to RDS file to which the function output will be saved
#'
scp_plot_cmod_markers_select <- function(
	mc_object,
	black_list = c(),
	sub_list_mc = NULL,
	gene_list = NULL,
	order_genes = TRUE,
	gene_annot_file = NULL,
	annot_header = FALSE,
	gene_font_size = 4,
	clust_ord = NULL,
	per_clust_genes = 20,
	gene_min_fold = 2,
	transverality_N = ncol(mc_object@mc_fp),
	transverality_thr = 1.4,
	transv_excluded_mc = NULL,
	output_file=NULL
) {
	
	# load gene annotations
	if (!is.null(gene_annot_file) & "character" %in% class(gene_annot_file)) {
		annot = read.table(gene_annot_file, header=annot_header, sep="\t", fill=TRUE, quote="", row.names=1)
	} else if (!is.null(gene_annot_file) & "data.frame" %in% class(gene_annot_file)) {
		annot = gene_annot_file
	}
	
	
	# expression matrix
	if (is.null(sub_list_mc)) {
		niche_geomean_n= mc_object@mc_fp
	} else {
		niche_geomean_n= mc_object@mc_fp[,sub_list_mc]
		clust_ord=sub_list_mc
	}
	
	# exclude genes with fc < gene_min_fold
	genes=unique(as.vector(unlist(apply(niche_geomean_n, 2, function(x) names(head(sort(-x[x > gene_min_fold]),n=per_clust_genes))))))
	
	# select genes for plotting
	# if (!is.null(black_list)) 
	#   black_list <- vector("character")
	if (is.null(gene_list)){
		# exclude blacklisted genes
		genes=setdiff(genes, black_list)
		message(sprintf("Excluded %s blacklisted genes", length(black_list)))
		# exclude genes with transversality > transverality_N
		transversal_genes=names(which(
			apply(
				niche_geomean_n[,setdiff(as.character(colnames(niche_geomean_n)),transv_excluded_mc)],
				1,
				function(x) sort(x,decreasing=TRUE)[transverality_N] > transverality_thr
			)
		))
		genes=setdiff(genes, transversal_genes)
	} else {
		# plot only genes in gene list, if it is specified
		genes <- gene_list[gene_list %in% genes]
	}
	genes=genes[genes %in% rownames(niche_geomean_n)]
	message("Will use ",length(genes)," genes")
	
	mat_niche <- niche_geomean_n[genes,]
	
	# if cluster order is not specified, do hierarchical clustering
	if (is.null(clust_ord)) {
		message("Recomputing cell ord")
		hc1 = hclust(dist(cor(mat_niche,method="pearson")), "ward.D2")
		clust_ord = as.character(hc1$order)
		scr_tmp_niche_order <- as.character(hc1$order)
	}
	
	# if gene order is TRUE, order genes
	if (order_genes){
		message("Ordering genes")
		gene_ord = order(apply(mat_niche[,as.character(clust_ord)],1,function(x) which.max(rollmean(x,1))))
	} else {
		gene_ord= 1:nrow(mat_niche)
	}
	gene_ord <- rev(gene_ord)
	
	# gene labels
	if (!is.null(gene_annot_file)) {
		
		gene_labels_0 <- genes[gene_ord]
		
		message("Genes: ", head(gene_labels_0), "...")
		
		gene_labels_1 <- as.character(annot[genes[gene_ord],2])
		message("Gene labels: ", head(gene_labels_0), "...")
		bad_labels <- gene_labels_1 %in% c("","-"," ") | is.na(gene_labels_1)
		message(sum(bad_labels), " bad gene labels")
		gene_labels_1[bad_labels] <- gene_labels_0[bad_labels]
		# OMIT GENE ANNOTATION SHORTENING: truncation + padding done in the actual plotting functions
		# long_labels <- nchar(gene_labels_1)>gene_chr_limit
		# message(sum(long_labels), " long gene labels")
		# gene_labels_1[long_labels] <- paste0(substr(gene_labels_1[long_labels],1,gene_chr_limit-3),"...")
		names(gene_labels_1) <- genes[gene_ord]
		
		gene_labels_3 <- as.character(annot[genes[gene_ord],1])
		bad_labels <- gene_labels_3 %in% c("","-"," ") | is.na(gene_labels_3)
		gene_labels_3[bad_labels] <- genes[gene_ord][bad_labels]
		# long_labels <- nchar(gene_labels_3)>gene_chr_limit
		# gene_labels_3[long_labels] <- paste0(substr(gene_labels_3[long_labels],1,gene_chr_limit-3),"...")
		gene_labels_2 <- ifelse(gene_labels_0 == gene_labels_3, gene_labels_0, paste(gene_labels_0,gene_labels_3, sep=" "))
		names(gene_labels_2) <- genes[gene_ord]
		
	} else {
		
		gene_labels_0 <- genes[gene_ord]
		gene_labels_1 <- genes[gene_ord]
		gene_labels_2 <- genes[gene_ord]
		
	}
	
	# return objects necessary for plotting
	marker_data_list = list(
		genes = genes,
		gene_ord = gene_ord,
		clust_ord = clust_ord,
		niche_geomean_n = niche_geomean_n,
		gene_labels_1 = gene_labels_1,
		gene_labels_2 = gene_labels_2
	)
	
	if (!is.null(output_file)) saveRDS(marker_data_list, output_file)
	
	return(marker_data_list)
	
}


#' Convenience function to create the vectors of gene name color and indexes to modify (used by the heatmap plotting functions)
#'
scp_highlight_genes_function = function(gene_labels, highlight_genes, highlight_color = "blue") {
	
	# create vector of colors (same length as input labels)
	gene_font_col = rep("black",length(gene_labels))
	highlight_ixs = which(gene_labels %in% highlight_genes)
	gene_font_col [ highlight_ixs ] = highlight_color
	
	return(list(gene_font_col = gene_font_col, highlight_ixs = highlight_ixs ))
	
}



#' Plot heatmap of gene expression for metacells
#'
#' Plots heatmap of gene expression fold change for metacells. Requires selecting markers and ordering first with `scp_plot_cmod_markers_select`
#'
#' @param marker_data_list character, list object produced by the heatmap pre-processing function `scp_plot_cmod_markers_select`
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param show_gene_names logical, show gene names as heatmap row names
#'    (default: FALSE)
#' @param highlight_genes vector of gene names to highlight in the heatmap (in a different color).
#' @param highlight_color highlight color for genes in `highlight_genes`; default is "blue"
#' @param omit_unhighlighted if TRUE, do not plot any other gene name other than the ones in `highlight_genes` (default: FALSE)
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param gene_chr_limit numeric, limit gene annotations to given number of charaters
#' @param clust_row data.frame with gene colour annotations in columns,
#'    rownames should be `marker_data_list$gids`; 
#'    if `NULL` (default), row colour annotation bar is not printed
#' @param clust_row_color named list of color mappings for annotations in `clust_row`; 
#'    names should be `colnames(clust_row)`
#' @param clust_col either a character vector or data.frame with cell colour annotations;
#'    if `clust_col` is a vector, it should specify colors to assign to metacells in the order given by `gene_list$clust_ord`,
#'    and if it is named, the names should be metacell names;
#'    if `clust_col` is a data.frame, rownames should be metacell names, annotations (not colors) should be in columns, 
#'    and color mappings should be given by `clust_col_color` argument;
#'    if `NULL` (default), cluster colour annotation bar is not printed
#' @param clust_col_color named list of color mappings for annotations in `clust_col`; 
#'    names should be `colnames(clust_col)`
#' @param clust_bars numeric, optional values to be ploted as bars annotation on top of
#'    heatmap columns, length should be same as `length(mc@mc_fp)` (default:  NULL)
#' @param clust_bars_color` character, color for barplots, can be either single color
#'    or a (named) vector (names should be `gene_list$clust_ord`)
#' @param clust_anno_size unit, height of the column annotation bar (default: `unit(1,"mm")`)
#' @param show_mc_names logical, show metacell names as heatmap column names (default: TRUE) 
#' @param use_raster logical, whether to rasterise heatmap bodies (default: TRUE) 
#' @param mc_font_size numeric, size of the metacell names (default: 4)
#' @param heatmap_colors vector of colors to use for heatmap coloring function (low to high; default: `c("white","gray99","orange","orangered2","#520c52")`)
#' @param max_expression_fc,min_expression_fc numeric, max and minimum expression values to scale metacell heatmap coloring to (default: 5 and 0). It's applied AFTER any transformation of the data.
#' @param add_expression_constant numeric, a psuedocount to add to the expresion matrix (default 1; useful if `use_log2` = TRUE). It's applied BEFORE transformation of the data.
#' @param transformation_fun, the name of a function used to transform the data (default is `log2`; to omit, set to NULL)
#'
scp_plot_cmod_markers_mc <- function(
	marker_data_list,
	mat_label = "FC",
	output_file = NULL,
	height = 10,
	width = 5,
	res = NA,
	show_heatmap_legend = TRUE,
	show_gene_names = FALSE,
	highlight_genes = NULL,
	highlight_color = "blue",
	omit_unhighlighted = FALSE,
	gene_font_size = 5,
	clust_row = NULL, clust_row_color=NULL,
	clust_col = NULL, clust_col_color=NULL,
	clust_bars = NULL, clust_bars_color="grey",
	clust_anno_size = unit(4,"mm"), clust_anno_gap = unit(1,"mm"),
	show_mc_names = TRUE, mc_font_size = 5,
	heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
	max_expression_fc = 5,
	min_expression_fc = 0,
	add_expression_constant = 1,
	transformation_fun = log2,
	gene_chr_limit = 70,
	use_raster = TRUE,
	verbose=FALSE,
	print_border=TRUE,
	show_clust_borders=TRUE
) {
	
	require("ComplexHeatmap")
	require("stringr")
	
	# PLOT METACELL PROFILE
	message("Plotting metacell expression")
	
	# get variables necessary for hm definition (from previous function call)
	genes = marker_data_list$genes
	gene_ord = marker_data_list$gene_ord
	clust_ord = marker_data_list$clust_ord
	niche_geomean_n = marker_data_list$niche_geomean_n
	gene_labels_1 = marker_data_list$gene_labels_1
	gene_labels_2 = marker_data_list$gene_labels_2
	
	# truncate left-side annotations
	gene_labels_1 = stringr::str_trunc(gene_labels_1, gene_chr_limit)
	gene_labels_1 = stringr::str_pad(gene_labels_1, gene_chr_limit, side="right")
	# truncate right annotations
	gene_labels_2 = stringr::str_trunc(gene_labels_2, gene_chr_limit)
	gene_labels_2 = stringr::str_pad(gene_labels_2, gene_chr_limit, side="left")
	
	# define matrix per metacell, based on geometric means
	if (!is.null(transformation_fun)) {
		mat1 = transformation_fun( niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + add_expression_constant )
		mat1 = pmin( mat1, max_expression_fc )
		mat1 = pmax( mat1, min_expression_fc )
	} else {
		mat1 = niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + add_expression_constant
		mat1 = pmin( mat1, max_expression_fc )
		mat1 = pmax( mat1, min_expression_fc )
	}
	
	# create gene annotations
	# if we have a vector of
	if (!is.null(highlight_genes)) {
		highlight_list = scp_highlight_genes_function(gene_labels = marker_data_list$gene_labels_1, highlight_genes = highlight_genes, highlight_color = highlight_color)
		gids = highlight_list$highlight_ixs
		gene_font_col = highlight_list$gene_font_col
	} else {
		gene_font_col = rep("black", length(marker_data_list$gene_labels_1))
	}
	
	if (show_gene_names) {
		if (!is.null(highlight_genes) & omit_unhighlighted) {
			# if we have a vector of genes to highlight and we don't want to plot the other genes...
			message("Gene annots highlights")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"), 
				gene = anno_mark(which="row", side="right", at=gids, labels=gene_labels_1[gids], labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"),
				gene = anno_mark(which="row", side="left", at=gids, labels=gene_labels_2[gids],  labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
		} else {
			# default scenario: plot all genes (and if there are genes to highlight, paint them in a different color)
			message("Gene annots")
			if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_1, location = 0, just = "left", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
			if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_2, location = 1, just = "right", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
		}
	} else {
		row_ha_left = ComplexHeatmap::HeatmapAnnotation(
			which = "row", empty = anno_empty(which = "row", border = FALSE)
		)
		row_ha_right = row_ha_left
	}
	
	# gene colour labels
	if (!is.null(clust_row)) {
		message("Row colour annots...")
		if (!is.null(clust_row_color)) {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1),,drop=FALSE], col = clust_row_color,
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		} else {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1),,drop=FALSE], 
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		}
		row_ha_left = c(row_ha_left, row_col_ha, gap = clust_anno_gap)
		row_ha_right = c(row_col_ha, row_ha_right, gap = clust_anno_gap)
	}
	
	
	message("Expression colors...")
	col_fun = circlize::colorRamp2(
		breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)), 
		colors = heatmap_colors)
	# shades = col_fun(seq(from = min_expression_fc, to = max_expression_fc, length.out = 40))
	
	# mc labels
	message("Metacell labels...")
	collabs <- colnames(mat1)
	if (show_mc_names) { 
		#collabs <- rep("",length(collabs))
		column_lab_ha = ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			LAB = anno_text(which = "column", collabs, gp = gpar(fontsize = mc_font_size, rot=90)),
			height=unit(2,"mm")
		)
	} else {
		column_lab_ha = ComplexHeatmap::HeatmapAnnotation(
			which = "column", empty = anno_empty(which = "column", border = FALSE),
			height = clust_anno_gap 
		)
	}
	top_column_ha = c(column_lab_ha)
	bottom_column_ha = c(column_lab_ha)
	
	# column color annotation
	if (!is.null(clust_col)) {
		message("Column colour annots...")
		
		if (class(clust_col) == "character") {
			
			if (is.null(names(clust_col))) {
				names(clust_col) <- clust_ord
			} else {
				if (!all(names(clust_col) %in% clust_ord))
					stop("Colour and cluster names do not match!")
				clust_col <- clust_col[clust_ord]
			}
		  annot_labs = as.character(unique(names(clust_col)))
			column_col_ha = ComplexHeatmap::HeatmapAnnotation(
				which = "column",
				"cluster" = colnames(mat1),
				col = list("cluster" = clust_col),
				border = TRUE,
				simple_anno_size = clust_anno_size,
				height = unit(1,"mm"),
				show_annotation_name = TRUE, show_legend = FALSE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size),
				annotation_legend_param = list(cluster = list(
				  at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120))
				))
			)
			
		} else if ("data.frame" %in% class(clust_col)) {
		  #annotation_legend_order <- sapply(colnames(clust_col_df), function(cn) {
		  #  annot_labs = as.character(unique(clust_col_df[,cn]))
		  #  list(at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120)))
		  #}, USE.NAMES = TRUE, simplify = FALSE)
			if (!is.null(clust_col_color)) {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col, col = clust_col_color,
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
					annotation_name_gp = gpar(fontsize = mc_font_size)
					#annotation_legend_param = annotation_legend_order
				)
			} else {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col, 
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
					annotation_name_gp = gpar(fontsize = mc_font_size)
					#annotation_legend_param = annotation_legend_order
				)
			}
			
		}
		
		top_column_ha <- c(column_col_ha,top_column_ha, gap = clust_anno_gap)
		bottom_column_ha = c(bottom_column_ha,column_col_ha, gap = clust_anno_gap)
	}
	
	# barplots
	if (!is.null(clust_bars)) {
		
		if (is.null(names(clust_bars)))
			names(clust_bars) <- as.character(clust_ord)
		anno_bar <- clust_bars[as.character(clust_ord)]
		baxl <- range(clust_bars)
		
		column_bar_ha <- ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			BAR = anno_barplot(
				anno_bar, height = 3 * clust_anno_size, bar_width = 0.9,
				gp = gpar(fill = clust_bars_color, col = clust_bars_color, fontsize = mc_font_size),
				axis_param = list(gp = gpar(fontsize = mc_font_size), at = baxl, labels = baxl)
			),
			show_annotation_name = FALSE, show_legend = FALSE, gap = clust_anno_gap
		)
		top_column_ha = c(column_bar_ha,top_column_ha, gap = clust_anno_gap)
		
	}
	
	# expression heatmap
	h1 = ComplexHeatmap::Heatmap(
		mat1, name = mat_label, col = col_fun, use_raster = use_raster,
		cluster_rows = FALSE, cluster_columns = FALSE,
		width = width,
		height = height,
		column_title = sprintf( "%i metacells", ncol(mat1) ),
		row_title = sprintf( "%i marker genes", nrow(mat1) ),
		show_column_names = FALSE,
		show_row_names = FALSE,
		right_annotation = row_ha_right,
		left_annotation = row_ha_left,
		top_annotation = top_column_ha,
		bottom_annotation = bottom_column_ha,
		column_names_gp = gpar(fontsize = mc_font_size),
		show_heatmap_legend = show_heatmap_legend,
		border = print_border
	)
	
	# save figure
	plotting_function(output_file, width, height, res, EXP = {
		
		# draw heatmap
		draw(h1)
		
		# add cluster borders
		if (!is.null(clust_col) & show_clust_borders & class(clust_col) == "character") {
			mat2 <- rbind(clust_col[1][match(clust_ord, rownames(clust_col))])
		} else if (!is.null(clust_col) & show_clust_borders & class(clust_col) %in% "data.frame") {
			mat2 <- rbind(clust_col[,1][match(clust_ord, rownames(clust_col))])
		} else {
			mat2 <- rbind(clust_ord)
		}
		if (show_clust_borders) {
			change_clust <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body(mat_label, {
				for (i in change_clust) {
					grid.lines(x = i / ncol(mat2), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
				}
			})
		}
	})
	
	message("Metacell heatmap done")
	
}




#' Plot heatmap of gene expression for single cells
#'
#' Plots heatmap of gene expression fold change for metacells. Requires selecting markers and ordering first with `scp_plot_cmod_markers_select`
#'
#' @param marker_data_list character, list object produced by the heatmap pre-processing function `scp_plot_cmod_markers_select`
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param show_gene_names logical, show gene names as heatmap row names (default: FALSE)
#' @param highlight_genes vector of gene names to highlight in the heatmap (in a different color).
#' @param highlight_color highlight color for genes in `highlight_genes`; default is "blue"
#' @param omit_unhighlighted if TRUE, do not plot any other gene name other than the ones in `highlight_genes` (default: FALSE)
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param gene_chr_limit numeric, limit gene annotations to given number of charaters
#' @param clust_row data.frame with gene colour annotations in columns,
#'    rownames should be `marker_data_list$gids`; 
#'    if `NULL` (default), row colour annotation bar is not printed
#' @param clust_row_color named list of color mappings for annotations in `clust_row`; 
#'    names should be `colnames(clust_row)`
#' @param clust_col either a character vector or data.frame with cell colour annotations;
#'    if `clust_col` is a vector, it should be of the same length and in the same order as `gene_list$clust_ord`, 
#'    and if it is named, the names should be metacell names;
#'    if `clust_col` is a data.frame, it should have single cells as rownames, annotations (not colors) should be in columns, 
#'    and color mappings should be given by `clust_col_color` argument;
#'    if `NULL` (default), cluster colour annotation bar is not printed
#' @param clust_col_color named list of color mappings for annotations in `clust_col` when it is a data.frame; 
#'    names should be `colnames(clust_col)`
#' @param clust_bars numeric, optional values to be ploted as bars annotation on top of
#'    heatmap columns, length should be same as `length(mc@mc)` (default:  NULL)
#' @param clust_bars_color` character, color for barplots, can be either single color
#'    or a (named) vector (names should be `gene_list$clust_ord`)
#' @param clust_anno_size unit, height of the column annotation bar (default: `unit(1,"mm")`)
#' @param show_mc_names logical, show metacell names as heatmap column names
#'    (default: TRUE) ~NOT IMPLEMENTED~
#' @param use_raster logical, whether to rasterise heatmap bodies (default: TRUE) 
#' @param mc_font_size numeric, size of the metacell names (default: 4)
#' @param heatmap_colors vector of colors to use for heatmap coloring function (low to high; default: `c("white","white","orange","red","purple","black")`)
#' @param min_expression_fc,max_expression_fc numeric, max and min expression values to scale single cell heatmap coloring to (default: 0 and 5)
#'
scp_plot_cmod_markers_sc <- function(
	marker_data_list,
	mc_object,
	mat_object,
	output_file = NULL,
	height = 10,
	width = 5,
	res = NA,
	show_heatmap_legend = TRUE,
	show_gene_names = FALSE,
	highlight_genes = NULL,
	highlight_color = "blue",
	omit_unhighlighted = FALSE,
	gene_font_size = 5,
	clust_row = NULL, clust_row_color = NULL,
	clust_col = NULL, clust_col_color = NULL,
	clust_bars = NULL, clust_bars_color="grey",
	clust_anno_size = unit(4,"mm"), clust_anno_gap = unit(1,"mm"),
	show_mc_names = TRUE, mc_font_size = 5,
	heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
	gene_chr_limit = 70,
	verbose=FALSE,
	smoothen = 5,
	min_expression_fc = 0,
	max_expression_fc = 5,
	use_raster = TRUE,
	print_border=TRUE,
	show_clust_borders = TRUE
	
) {
	
	# get variables necessary for hm definition (from previous function call)
	genes = marker_data_list$genes
	gene_ord = marker_data_list$gene_ord
	clust_ord = marker_data_list$clust_ord
	niche_geomean_n = marker_data_list$niche_geomean_n
	gene_labels_1 = marker_data_list$gene_labels_1
	gene_labels_2 = marker_data_list$gene_labels_2
	
	# truncate left-side annotations
	gene_labels_1 = stringr::str_trunc(gene_labels_1, gene_chr_limit)
	gene_labels_1 = stringr::str_pad(gene_labels_1, gene_chr_limit, side="right")
	# truncate right annotations
	gene_labels_2 = stringr::str_trunc(gene_labels_2, gene_chr_limit)
	gene_labels_2 = stringr::str_pad(gene_labels_2, gene_chr_limit, side="left")
	
	# directory to save output files to
	if (is.null(output_file)) {
		outdir <- getwd()
	} else {
		outdir <- dirname(output_file)
	}
	
	
	# PLOT SINGLE-CELL PROFILE
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
	
	# heatmap per metacell
	mat1sc <- pmin(lus_smoo, max_expression_fc)
	mat1sc <- pmax(lus_smoo, min_expression_fc)
	colnames(mat1sc) <- colnames(lus_smoo)
	
	# define matrix per metacell
	mat1 <- niche_geomean_n[genes[gene_ord],as.character(clust_ord)]
	
	# colors in heatmap
	col_fun = circlize::colorRamp2(
		breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)), 
		colors = heatmap_colors)
	# shades = col_fun(seq(from = min_expression_fc, to = max_expression_fc, length.out = 40))
	
	# mc labels
	top_column_ha <- HeatmapAnnotation(
		mclabstop = anno_empty(which = "column", border = FALSE, height = unit(2,"mm"))
	)
	bottom_column_ha <- HeatmapAnnotation(
		mclabsbottom = anno_empty(which = "column", border = FALSE, height = unit(2,"mm"))
	)
	
	# Column colour annotation
	if (!is.null(clust_col)) {
		message("Column colour annots...")
		
		if (class(clust_col) == "character") {
			
			if (is.null(names(clust_col))) {
				names(clust_col) = clust_ord
			}
		  annot_labs = as.character(unique(names(clust_col)))
			column_col_ha <- HeatmapAnnotation(
				which = "column",
				"cluster" = cells_clusts, col = list("cluster" = clust_col),
				border = c(TRUE),
				simple_anno_size = clust_anno_size,
				show_annotation_name = TRUE, show_legend = FALSE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size),
				annotation_legend_param = list(cluster = list(
				  at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120))
				))
			)
			top_column_ha <- c(column_col_ha,top_column_ha, gap = clust_anno_gap)
			bottom_column_ha <- c(bottom_column_ha,column_col_ha, gap = clust_anno_gap)
			
		} else if ("data.frame" %in% class(clust_col)) {
		  annotation_legend_order <- sapply(colnames(clust_col_df), function(cn) {
		    annot_labs = as.character(unique(clust_col_df[,cn]))
		    list(at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120)))
		  }, USE.NAMES = TRUE, simplify = FALSE)
			if (!is.null(clust_col_color)) {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col[colnames(mat1sc),,drop=FALSE], col = clust_col_color,
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
					annotation_name_gp = gpar(fontsize = mc_font_size),
					annotation_legend_param = annotation_legend_order
				)
			} else {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col[colnames(mat1sc),,drop=FALSE], 
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
					annotation_name_gp = gpar(fontsize = mc_font_size),
					annotation_legend_param = annotation_legend_order
				)
			}
			top_column_ha <- c(column_col_ha,top_column_ha, gap = clust_anno_gap)
			bottom_column_ha = c(bottom_column_ha,column_col_ha, gap = clust_anno_gap)
		}
		
	}	
	# barplot annotation
	if (!is.null(clust_bars)) {
		
		cell_cols <- unlist(mapply(rep, clust_bars_color, n_cells_cluster, USE.NAMES=FALSE))
		cell_col <- structure(cell_cols, names=cell_order)
		
		if (is.null(names(clust_bars)))
			names(clust_bars) <- as.character(cell_order)
		
		anno_bar <- clust_bars[as.character(cell_order)]
		baxl <- range(clust_bars)
		
		column_bar_ha <- ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			BAR = anno_barplot(
				anno_bar, height = 3 * clust_anno_size, bar_width = 1, 
				gp = gpar(fill = cell_col, col = cell_col, fontsize = mc_font_size),
				axis_param = list(gp = gpar(fontsize = mc_font_size), at = baxl, labels = baxl)
			),
			show_annotation_name = FALSE, show_legend = FALSE, gap = clust_anno_gap
		)
		top_column_ha = c(column_bar_ha,top_column_ha)
		
	}
	
	# create gene annotations
	# if we have a vector of genes to highlight...
	if (!is.null(highlight_genes)) {
		highlight_list = scp_highlight_genes_function(gene_labels = names(marker_data_list$gene_labels_1), highlight_genes = highlight_genes, highlight_color = highlight_color)
		gids = highlight_list$highlight_ixs
		gene_font_col = highlight_list$gene_font_col
	} else {
		gene_font_col = rep("black", length(marker_data_list$gene_labels_1))
	}
	
	if (show_gene_names) {
		if (!is.null(highlight_genes) & omit_unhighlighted) {
			# if we have a vector of genes to highlight and we don't want to plot the other genes...
			message("Gene annots highlights")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"), 
				gene = anno_mark(which="row", side="right", at=gids, labels=gene_labels_1[gids], labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"),
				gene = anno_mark(which="row", side="left", at=gids, labels=gene_labels_2[gids],  labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
		} else {
			# default scenario: plot all genes (and if there are genes to highlight, paint them in a different color)
			message("Gene annots")
			if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_1, location = 0, just = "left", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
			if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_2, location = 1, just = "right", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
		}
	} else {
		row_ha_left = ComplexHeatmap::HeatmapAnnotation(
			which = "row", empty = anno_empty(which = "row", border = FALSE)
		)
		row_ha_right = row_ha_left
	}
	
	# gene colour labels
	if (!is.null(clust_row)) {
		message("Row colour annots...")
		if (!is.null(clust_row_color)) {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1sc),,drop=FALSE], col = clust_row_color,
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		} else {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1sc),,drop=FALSE], 
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		}
		row_ha_left = c(row_ha_left, row_col_ha, gap = clust_anno_gap)
		row_ha_right = c(row_col_ha, row_ha_right, gap = clust_anno_gap)
	}
	
	# expression heatmap
	h1sc <- Heatmap(
		mat1sc, name = "sc_expression", col = col_fun, use_raster = use_raster,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		show_column_names = FALSE,
		show_row_names = FALSE,
		width = width,
		height = height,
		column_title = sprintf( "%i cells", ncol(mat1sc) ),
		row_title = sprintf( "%i marker genes", nrow(mat1sc) ),
		right_annotation = row_ha_right,
		left_annotation = row_ha_left,
		top_annotation = top_column_ha,
		bottom_annotation = bottom_column_ha,
		show_heatmap_legend = show_heatmap_legend,
		border = print_border
	)
	hlistsc <- h1sc
	
	# save figure
	
	plotting_function(output_file, width, height, res, EXP={
		# heatmap
		# draw(hlistsc, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
		# drop all this artificial padding...
		draw(hlistsc)
		
		# add labels of metacells
		if (!is.null(clust_col)) {
			if (class(clust_col) == "character") {
				mat2sc <- rbind(clust_col[match(cells_clusts, names(clust_col))])
			} else if ("data.frame" %in% class(clust_col)) {
				mat2sc <- rbind(clust_col[colnames(mat1sc),colnames(clust_col)[1]])
			}
		} else {
			mat2sc <- rbind(cells_clusts)
		}
		if (is.null(colnames(mat2sc)))
			colnames(mat2sc) <- cells_clusts
		change_clust_sc <- which(sapply(2:ncol(mat2sc), function(i) mat2sc[,i] != mat2sc[,i - 1]))
		change_mc_sc <- c(
			which(sapply(2:ncol(mat2sc), function(i) colnames(mat2sc)[i] != colnames(mat2sc)[i - 1])),
			ncol(mat2sc)
		)
		if (show_clust_borders) {
			decorate_heatmap_body("sc_expression", {
				for (i in change_clust_sc) {
					grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
				}
				for (i in change_mc_sc) {
					grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
				}
			})
		}
		
		.add_mc_labels <- function(pos,labs) {
			for (j in 1:length(pos)) {
				i <- pos[j]
				iprev <- ifelse(j == 1,0,pos[j - 1])
				nt <- ncol(mat2sc)
				grid.text(label = labs[j], x = i / nt - (i / nt - iprev / nt) / 2, y = 0.5, gp = gpar(fontsize = mc_font_size), rot=90)
			}
		}
		if (show_mc_names) {
			decorate_annotation("mclabstop", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
			decorate_annotation("mclabsbottom", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
		}
	})
	ht_opt(RESET = TRUE)
	
	message("Single-cell heatmap done")
	
}



#' Add cell annotations to metacell object
#'
#' @param input_table path to cell type annotation table (tab-separated,)
#' @param mc_object loaded metacell object (`gMCCov` class)
#'
#' @return data.frame with metacells' cell_type annotations. Also saves the plot to disk.
#'
sca_load_cell_type_table = function(input_table,mc_object){
	
	cell_type_table = read.table(input_table,h=TRUE,sep="\t",comment.char="",stringsAsFactors=FALSE)
	rownames(cell_type_table)=cell_type_table$metacell
	
	sc_ct_label=as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
	names(sc_ct_label)=names(mc_object@mc)
	
	clust_cols=as.character(cell_type_table[,"color"])
	names(clust_cols)=rownames(cell_type_table)
	
	cells_cols=cell_type_table[as.character(mc_object@mc),"color"]
	cells_cols=as.character(cells_cols)
	names(cells_cols)=names(mc_object@mc)
	
	mc_object@colors=clust_cols
	return(list(ct_table=cell_type_table,mc_color=clust_cols,sc_ct_label=sc_ct_label,sc_color=cells_cols))
	
}


sca_reorder_mc2d = function(mc, mc2d) {
	
	# keep mc order from mc
	mc2d_r = mc2d
	
	# reorder mc coordinates
	mc2d_r@mc_x = mc2d_r@mc_x [ colnames(mc@mc_fp) ]
	mc2d_r@mc_y = mc2d_r@mc_y [ colnames(mc@mc_fp) ]
	
	# # reorder mc graph
	# mc_dict = colnames(mc@mc_fp)
	# names(mc_dict) = 1:length(colnames(mc@mc_fp))
	# mc2d_r@graph$mc1 = as.character(mc_dict [ mc2d_r@graph$mc1 ])
	# mc2d_r@graph$mc2 = as.character(mc_dict [ mc2d_r@graph$mc2 ])
	
	return(mc2d_r)
	
}



#' Plot single cell/metacell 2D projection
#'
#' @param mc2d 2D projection from `mc2d` object
#' @param mc loaded metacell object (`gMCCov` class)
#' @param mc_colors optional named vector of metacell colors; if `NULL` (default), `mc@colors` are used
#' @param cell_colors optional named vector of individual cell colors; if `NULL` (default), metacell color is used
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' 
#' @return data.frame with metacells' cell_type annotations. Also saves the plot to disk.
#'
scp_plot_sc_2d = function(
	mc2d,
	mc,
	mc_colors=NULL,
	cell_colors=NULL,
	plot_edges=FALSE,
	plot_mcs=FALSE,
	plot_mc_name=FALSE,
	width=12,
	height=12,
	res=NA,
	output_file=NULL,
	cex_mc=3,
	cex_sc=0.75,
	alpha_mc=0.8,
	alpha_sc=0.4,
	do_axes=FALSE) {
	
	# get colors of metacells
	if (is.null(mc_colors)) mc_colors = mc@colors
	if (is.null(names(mc_colors))) names(mc_colors) <- colnames(mc@mc_fp)
	# get colors of individual cells
	if (is.null(cell_colors)) cell_colors = mc_colors[ as.character(mc@mc[names(mc2d@sc_x)]) ]
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		# determine plot max/min
		if (plot_mcs) {
			
			xlim=c(min(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE), max(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE), max(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE))
			
		} else {
			
			xlim=c(min(mc2d@sc_x, na.rm = TRUE), max(mc2d@sc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, na.rm = TRUE), max(mc2d@sc_y, na.rm = TRUE))
			
		}
		
		# plot individual cells
		plot(
			mc2d@sc_x,
			mc2d@sc_y,
			pch = 19,lwd = 0,
			cex = cex_sc,
			col = alpha(cell_colors,alpha_sc),
			xlim = xlim,
			ylim = ylim,
			xlab = NA, ylab = NA, axes=do_axes
		)
		
		# plot edges between metacells?
		if (plot_edges) {
			fr = as.character(mc2d@graph$mc1)
			to = as.character(mc2d@graph$mc2)
			segments(
				x0=mc2d@mc_x[fr],
				y0=mc2d@mc_y[fr],
				x1=mc2d@mc_x[to],
				y1=mc2d@mc_y[to],
				col="gray70")
		}
		
		# plot metacells?
		if (plot_mcs) {
			points(
				mc2d@mc_x,
				mc2d@mc_y,
				pch=19,lwd=0.5,
				cex=cex_mc,
				col=alpha(mc_colors,alpha_mc))
		}
		
		# plot metacell ids?
		if (plot_mc_name) {
			text(
				mc2d@mc_x,
				mc2d@mc_y,
				#cex=cex_mc,
				labels=names(mc2d@mc_x))
		}
		title(
			sub = sprintf(
				"2D projection\nn = %i cells | n = %i metacells", 
				length(mc2d@sc_x),
				length(mc2d@mc_x)
			)
		)
		
	})
}



#' Plot heatmaps and barplots for predetermined gene sets
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param mc_counts UMI counts per metacell (`dataframe` class, genes are rows)
#' @param markers_file path to file containing gene annotations (first column is gene name, others are optional annotations). You can supply a vector of gene names as well (which will then be unnannotated) or a dataframe fitting the markers_file description.
#' @param output_file_heatmap,output_file_barplot paths to file to which the plots will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param do_log_fp whether to use log2 scale for heatmaps (default: FALSE)
#' @param use_raster logical, whether to rasterise heatmap bodies (default: TRUE) 
#' @param max_expression_fc,min_expression_fc values to cap colorscale
#' @param pmin fc values higher than this will be capped to this value (a bit improper, you'd be better off capping the colorscale with max_expression_fc instead).
#' @param set_pmin whether to find a good `pmin` value or not (default is FALSE; a bit improper, you'd be better off capping the colorscale with max_expression_fc instead).
#' 
scp_barplot_heatmap_markers = function(
	mc_object,
	mat_object,
	mc_counts,
	markers_file,
	output_file_heatmap = NULL,
	output_file_barplot = NULL,
	width = 10, height = NULL, res = NA,
	use_raster = FALSE,
	heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
	T_totumi = 15,
	do_log_fp = FALSE,
	clust_order= NULL,
	mc_color= NULL,
	addumidata = TRUE,
	print_barplots= FALSE,
	# pmin = 4,
	min_gene_fc = 2,
	# set_pmin = FALSE,
	column_names = NULL,
	table_file = NULL,
	sort_markers = TRUE,
	gene_font_size = 5,
	mc_font_size = 5,
	min_expression_fc = 0,
	max_expression_fc = 4,
	gene_chr_limit = 70,
	clust_anno_size = unit(2,"mm")
	
) {
	
	# if no cluster reordering, get from mc object
	if (is.null(clust_order)) {
		clust_order=colnames(mc_object@mc_fp)
	}
	
	clusts=mc_object@mc
	footprint_table=mc_object@mc_fp[,clust_order]
	
	if (do_log_fp) {
		footprint_table = log2(footprint_table)
	}
	
	# define markers
	# first, load full table of markers
	# either a dataframe where genes are indicated as row names, or a path to such a table
	if (class(markers_file) == "character" & length(markers_file) == 1) {
		markers_table = read.table(markers_file, header=FALSE, row.names=1, comment.char="", sep="\t")
	} else if (class(markers_file) == "character" & length(markers_file) > 1) {
		markers_table = data.frame(row.names = markers_file, annot = markers_file)
	} else if (class(markers_file) == "data.frame") {
		markers_table = markers_file
	} else {
		message("`markers_file` is not a path to a file or a data.frame, but I'll proceed anyway and see what comes out of it.")
	}
	# format
	if (ncol(markers_table) > 1) {
		collapsed_markers_annotations = apply(markers_table, 1, function(x) paste(x, collapse = " | "))
		markers_table = data.frame(row.names = rownames(markers_table), annotation = collapsed_markers_annotations)
	}
	# keep markers present in the footprint table
	markers = intersect(rownames(markers_table),rownames(footprint_table))
	
	# keep markers with min expression and FC
	f_marker_cov = rowSums(mc_counts[markers,]) > T_totumi
	f_marker_fc = apply(mc_object@mc_fp[markers,],1,max) > min_gene_fc
	markers = intersect(names(which(f_marker_cov)),names(which(f_marker_fc)))
	markers = intersect(markers,rownames(footprint_table))
	
	# define heatmap height based on number of markers
	# height = pmin(length(markers)*0.7,10),10)
	if (is.null(height)) {
		height = round(length(markers) / 20) + 4
	}
	
	# calculate umifrac per metacell
	mc_umifrac = sca_mc_gene_umifrac_noobj(mc_counts = mc_counts, multiplying_factor = 10000)
	
	# sort markers based on number of UMIs
	marker_fp = footprint_table[markers,as.character(clust_order)]
	if (sort_markers & length(markers) > 1) {
		markers_sorted=rev(markers[as.numeric(order(apply(marker_fp, 1,function(x) which.max(rollmean(x,1)))))])
	} else {
		markers_sorted=rev(markers)
	}
	
	if (length(markers_sorted) > 1) {
		
		markers_tot_umi = rowSums(mc_counts[markers_sorted,])
		marker_max_umi_frac = round(apply(mc_counts[markers_sorted,] / markers_tot_umi,1,max),2)
		message(sprintf("%i markers pass filters", length(markers_sorted)))
		
		# if (set_pmin){
		# 	pmin=pmin(quantile(footprint_table[markers_sorted,clust_order],0.995),4)
		# }
		# message(sprintf("I'm using PMIN %i", pmin))
		# marker_fp_to_plot = pmin(footprint_table[markers_sorted,clust_order],pmin)
		# marker_fp_to_plot_ordered = marker_fp_to_plot[markers_sorted,]
		marker_fp_to_plot = footprint_table[markers_sorted,clust_order]
		marker_fp_to_plot_ordered = marker_fp_to_plot[markers_sorted,]
		
		# prepare heatmap
		# first, create color map
		col_fun = circlize::colorRamp2(
			breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)), 
			colors = heatmap_colors)
		
		# second, get row annotations
		# left-side annotations (gene ids)
		if (addumidata) {
			labels_left = paste(markers_tot_umi, marker_max_umi_frac, markers_sorted, sep=" | ")
		} else {
			labels_left = markers_sorted
		}
		labels_left = stringr::str_trunc(labels_left, gene_chr_limit)
		labels_left = stringr::str_pad(labels_left, gene_chr_limit, side="left")
		row_ha_left = ComplexHeatmap::HeatmapAnnotation(
			which = "row", simple_anno_size = unit(1,"mm"), gap = unit(5,"mm"),
			log10_tot_umi = anno_barplot(log10(markers_tot_umi), which = "row", axis = TRUE, gp = gpar(col = NA, fill = "slategray3"), baseline = 0, border = FALSE),
			max_umifrac = anno_barplot(marker_max_umi_frac,      which = "row", axis = TRUE, gp = gpar(col = NA, fill = "slategray4"), ylim = c(0,1), border = FALSE),
			gene = anno_text(which = "row", labels_left, gp = gpar(fontsize = gene_font_size, col = "black"))
		)
		# right-side annotations (gene names)
		labels_right = markers_table[ markers_sorted, 1 ]
		labels_right = stringr::str_trunc(labels_right, gene_chr_limit)
		labels_right = stringr::str_pad(labels_right, gene_chr_limit, side="right")
		row_ha_right = ComplexHeatmap::HeatmapAnnotation(
			which = "row", simple_anno_size = unit(1,"mm"),
			gene = anno_text(which = "row", labels_right, gp = gpar(fontsize = gene_font_size, col = "black"))
		)
		
		# third, plot column annotations (cluster colors names and colors)
		# start with labels
		labels_columns = clust_order
		top_column_ha = ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			metacell = anno_text(which = "column", labels_columns, just = "left", location = 0, gp = gpar(fontsize = mc_font_size, rot=90))
		)
		bot_column_ha = top_column_ha
		# add color strip to column
		if (!is.null(mc_color)) {
			labeled_colors = mc_color
			names(labeled_colors) = labels_columns
			top_column_ha_color = ComplexHeatmap::HeatmapAnnotation(
				which = "column",
				MC = labels_columns,
				col = list(MC = labeled_colors),
				border = TRUE,
				simple_anno_size = clust_anno_size,
				height = unit(5,"mm"),
				show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
			)
			bot_column_ha = c(bot_column_ha, top_column_ha_color, gap = unit(1, "mm"))
			top_column_ha = c(top_column_ha_color, top_column_ha, gap = unit(1, "mm"))
		}
		
		
		# finally, draw heatmap
		message("Plotting annotated heatmap")
		matrix_title = ifelse(do_log_fp, "log2(fp)", "fp")
		h1 <- ComplexHeatmap::Heatmap(
			matrix = marker_fp_to_plot_ordered,
			name = matrix_title,
			col = col_fun,
			cluster_rows = FALSE,
			cluster_columns = FALSE,
			width = width,
			height = height,
			show_column_names = FALSE,
			show_row_names = FALSE,
			right_annotation = row_ha_right,
			left_annotation = row_ha_left,
			top_annotation = top_column_ha,
			bottom_annotation = bot_column_ha,
			column_names_gp = gpar(fontsize = mc_font_size),
			use_raster = use_raster,
			show_heatmap_legend = TRUE,
			border = TRUE
		)
		plotting_function(output_file_heatmap, width, height, res, EXP = {
			draw(h1)
		})
		
		# output table with annotations and per-metacell counts
		if (!is.null(table_file)) {
			
			message("Writing footprint table")
			# create dataframe
			m_to_write = cbind(
				markers_table[markers_sorted,1],
				round(footprint_table[markers_sorted,clust_order],3)
			)
			colnames(m_to_write) = c(colnames(markers_table)[1], clust_order)
			
			# write
			write.table(
				m_to_write,
				file = table_file,
				col.names=TRUE,
				row.names=TRUE,
				sep="\t",
				quote=FALSE)
			
		}
		
	} else {
		
		message("Omit footprint tables and heatmaps: one marker or less")
		
	}
	
	if (print_barplots) {
		plotting_function(output_file_barplot, width, 3, res=NA, EXP = {
			# print barplots of individual marker genes
			message("Plotting per-gene UMI frac barplots")
			for (nm in markers) {
				barplot(
					mc_umifrac[nm,],
					main=sprintf("%s\n%s",markers_table[nm,1], nm),
					border="white", col=mc_color,
					xlim = c( 1, ncol(mc_umifrac) ),
					las=2, space=0, ylab = ("UMI per 10k"))
			}
		})
	}
	
	
	if (length(markers_sorted) > 1) {
		return(marker_fp_to_plot_ordered)
	} else {
		return(NULL)
	}
	
}



#' Plot metacell size and UMI counts
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param mc_counts UMI counts per metacell (`dataframe` class, genes are rows)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param return_table logical; whether to return a data.frame with per-metacell size stats
#'
scp_plot_mc_size_counts = function(mc_object, mat_object, mc_counts = NULL, mc_color="blue", output_file = NULL, width = 24, height = 6, res = NA, return_table = FALSE) {
	
	# UMI counts per metacell
	if (is.null(mc_counts)) {
		mc_counts = sca_mc_gene_counts(mc_object, mat_object, 0)
	}
	# num cells per metacell
	mc_n_cells = table(mc_object@mc)
	# median cell size
	if ("dgCMatrix" %in% class(mat_object@mat)) {
		sc_counts = Matrix::colSums(mat_object@mat[,names(mc_object@mc)])
	} else {
		sc_counts = colSums(as.matrix(mat_object@mat[,names(mc_object@mc)]))
	}
	mc_median_cell_size = tapply(sc_counts, mc_object@mc, median)
	mc_mean_cell_size = tapply(sc_counts, mc_object@mc, mean)
	
	# keep metacell order from mc object
	mc_counts = mc_counts [ , colnames(mc_object@mc_fp) ]
	mc_n_cells = mc_n_cells [ colnames(mc_object@mc_fp) ]
	mc_median_cell_size = mc_median_cell_size [ colnames(mc_object@mc_fp) ]
	mc_object_mc_factor = factor(mc_object@mc, levels = colnames(mc_object@mc_fp))
	
	plotting_function(output_file, width, height, res, EXP = {
		
		# metacell size (num cells)
		barplot(
			mc_n_cells,
			col = mc_color, 
			border="white", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space= 0,
			ylab = "# cells" ,
			main = "# cells per metacell")
		
		boxplot(
			sc_counts + 1 ~ mc_object_mc_factor,
			col = mc_color, 
			border= "black", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space = 0, 
			ylab = "# UMI" ,
			log = "y",
			main = "# UMI per metacell")
		abline(h=0, lty=2, col="gray")
		
		# median umi per metacell
		barplot(
			mc_median_cell_size,
			col = mc_color, 
			border="white", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space=0, 
			ylab = "# UMI" ,
			main = "Median # UMI per metacell")
		
		# total umi per metacell
		barplot(
			colSums(mc_counts),
			col = mc_color, 
			border="white", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space=0, 
			ylab = "# UMI" ,
			main = "Total # UMI per metacell")
		
		# distrubution of foldchange values per metacell
		vioplot::vioplot(
			log2(apply(mc@mc_fp, 2, function(x) x[x > 0] )),
			col = mc_color, 
			border= "black", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space = 0, 
			ylab = "log2 fold change" ,
			main = "fold changes per metacell")
		abline(h=log2(c(1,1.5,2)), lty=2, col="gray")
		
	})
	
	# return table?
	if (return_table) {
		
		sta = data.frame(
			metacell = colnames(mc@mc_fp),
			num_cells = as.vector(mc_n_cells),
			tot_umis = colSums(mc_counts),
			cell_size_mean   = tapply(sc_counts, mc_object@mc, mean),
			cell_size_Q1     = tapply(sc_counts, mc_object@mc, function(i) quantile(i, 0.25)),
			cell_size_Q2     = mc_median_cell_size,
			cell_size_Q3     = tapply(sc_counts, mc_object@mc, function(i) quantile(i, 0.75)),
			log2_fc_mean     = log2(apply(mc@mc_fp, 2, function(i) mean(i[i > 1]))),
			log2_fc_Q1       = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.25))),
			log2_fc_Q2       = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.50))),
			log2_fc_Q3       = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.75))),
			log2_fc_Qp95     = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.95)))
			
		)
		
	}
	
}



#' Plot single cell/metacell 2D projection with single gene expression
#'
#' @param mc2d 2D projection from `mc2d` object
#' @param mc loaded metacell object (`gMCCov` class)
#' @param mat loaded single cell matrix object (`tgScMat` class)
#' @param gene_id string, gene id (rownames in `mat@mat` object)
#' @param sc_vector vector of single cell-level values to map. If NULL (default), the # of UMIs is taken from the `mat` object.
#' @param sc_scale indicate wether cell-level values are  unidirectional (`unidir`) or bidirectional (`bidir`); defaults to `unidir`
#' @param sc_transform a function to transform sc vector (default is NULL, for no transformation)
#' @param sc_min min value of sc vector (lower values are raised to this threshold); should be in the transformed scale if transformation is used
#' @param sc_max max value of sc vector (higher values are lowered to this threshold); should be in the transformed scale if transformation is used
#' @param sc_max_quant if sc_max is NULL, use this quantile to find the upper threshold (default is 0.75)
#' @param sc_label legend title label for the sc vector. Defaults to "cell UMI"
#' @param sc_zero_color a color to be used for zero values in the sc vector (useful to distinguish zero values from the rest); default to light gray.
#' @param do_umifrac_sc boolean; if set to TRUE, UMI fraction per 10k is calculated, and `sc_label` is set to "UMI/10k"
#' @param mc_vector vector of metacell-level values to map. If NULL (default), footprint values from `mc@mc_fp` are used and log2-transformed.
#' @param mc_scale indicate wether mc-level values are  unidirectional ("unidir") or bidirectional ("bidir"); defaults to "unidir"
#' @param mc_transform a function to transform mc vector (default is `log` transformation)
#' @param mc_min min value of mc vector (lower values are raised to this threshold); should be in the transformed scale if transformation is used
#' @param mc_max max value of mc vector (higher values are lowered to this threshold); should be in the transformed scale if transformation is used
#' @param mc_max_quant if mc_max is NULL, use this quantile to find the upper threshold (default is 0.75)
#' @param mc_label legend title label for the mc vector. Defaults to "mc log(fp)"
#' @param mc_zero_color a color to be used for zero values in the mc vector (useful to distinguish zero values from the rest); default to NULL
#' @param unidir_color_scale vector of colors that will be used for unidirectional color scales
#' @param bidir_color_scale vector of colors that will be used for bidirectional color scales
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param plot_mcs whether to plot metacells (default FALSE)
#' @param plot_mc_name whether to plot metacell names (default FALSE)
#' @param plot_edges whether to plot edges between metacells (default FALSE)
#' @param plot_cells_as_2d_density instead of plotting expression in single cells, create a 2d density plot from the coordinates of single cells (can work well with UMI counts/fraction data; default FALSE).
#' @param cex_mc size of metacell points (default 3)
#' @param cex_sc size of single cell points (default 0.75)
#' @param alpha_mc alpha value for metacell points (default 0.9)
#' @param alpha_sc alpha value for single cell points (default 0.7)
#' @param alpha_sc_2d alpha value for single cell 2d density (default 1)
#' @param do_legend_sc,do_legend_mc whether to add legends with expression scale for metacells and single cells (default TRUE)
#' @param do_axes whether to plot axes for the 2D projection (default FALSE)
#'
scp_plot_sc_2d_gene_exp = function(
	mc2d,
	mc,
	mat,
	gene_id,
	sc_vector = NULL,
	sc_scale = "unidir",
	sc_transform = NULL,
	sc_min = 0,
	sc_max = NULL,
	sc_max_quant = 0.75,
	sc_label = "cell UMI",
	sc_zero_color = "lightblue1",
	do_umifrac_sc = FALSE,
	mc_vector = NULL,
	mc_scale = "unidir",
	mc_transform = log2,
	mc_min = 0,
	mc_max = 2,
	mc_max_quant = 0.75,
	mc_label = "mc log2(fp)",
	mc_zero_color = NULL,
	unidir_color_scale = c("gray90","orange","orangered2","#520c52"),
	bidir_color_scale = c("midnightblue","dodgerblue3","deepskyblue","#b8e0ed","gray90","#eccac0","#ff8d36","#f34312","#8e0631"),
	plot_mcs=FALSE,
	plot_mc_name=FALSE,
	plot_edges=FALSE,
	plot_cells_as_2d_density=FALSE,
	title=NULL,
	width=12,
	height=12,
	res=NA,
	output_file=NULL,
	cex_mc=3,
	cex_sc=0.75,
	alpha_mc=0.9,
	alpha_sc=0.7,
	alpha_sc_2d=1,
	do_legend_sc=TRUE,
	do_legend_mc=TRUE,
	do_axes=FALSE) {

	# get expression vectors (if not already provided)
	if (is.null(sc_vector)) {
		sc_vector = mat@mat[gene_id,]
	}
	if (is.null(mc_vector)) {
		mc_vector = mc@mc_fp[gene_id,]
	}


	# apply transformations
	if (do_umifrac_sc) {
		sc_vector = sc_vector / apply(mat@mat, 2, function(x) sum(x, nar.rm = TRUE)) * 10000
		sc_label = "UMI/10k"
	}
	if (!is.null(sc_transform)) {
		sc_vector = sc_transform(sc_vector)
	}
	if (!is.null(mc_transform)) {
		mc_vector = mc_transform(mc_vector)
	}


	# if sc_max or mc_max values are NULL, get them from a high quantile in the mc_vector
	if (is.null(sc_max)) {
		sc_max = round(quantile(sc_vector[sc_vector > 0], sc_max_quant, na.rm=TRUE))
	}
	if (is.null(mc_max)) {
		mc_max = round(quantile(mc_vector[mc_vector > 0], mc_max_quant, na.rm=TRUE))
	}

	# apply min/max to sc vector
	sc_vector [ sc_vector < sc_min ] = sc_min
	sc_vector [ sc_vector > sc_max ] = sc_max
	# apply min/max to mc vector
	mc_vector [ mc_vector < mc_min ] = mc_min
	mc_vector [ mc_vector > mc_max ] = mc_max

	# create color palettes for single cells
	if (sc_scale == "unidir") {
		sc_color_fun = scales::col_numeric(palette=unidir_color_scale, domain=c(sc_min, sc_max))
	} else if (sc_scale == "bidir") {
		sc_color_fun = scales::col_numeric(palette=bidir_color_scale, domain=c(sc_min, sc_max))
	} else {
		message("`sc_scale` should be either `unidir` or `bidir`")
	}
	sc_color = sc_color_fun(sc_vector)
	sc_color_labels = seq(from=sc_min, to=sc_max, length.out=9)
	sc_color_legend = sc_color_fun(sc_color_labels)
	sc_color_labels = sprintf("%.2f",sc_color_labels)
	sc_color_labels [ length(sc_color_labels) ] = paste(">=", sc_color_labels [ length(sc_color_labels) ], sep="")
	sc_color_labels [ 1 ] = paste("<=", sc_color_labels [ 1 ], sep="")
	# zero color
	if (!is.null(sc_zero_color)) {
		sc_color [ sc_vector == 0 ] = sc_zero_color
		sc_color_legend [ sc_color_labels == 0 ] = sc_zero_color
	}

	# color palettes for metacells
	if (mc_scale == "unidir") {
		mc_color_fun = scales::col_numeric(palette=unidir_color_scale, domain=c(mc_min, mc_max))
	} else if (mc_scale == "bidir") {
		mc_color_fun = scales::col_numeric(palette=bidir_color_scale, domain=c(mc_min, mc_max))
	} else {
		message("`sc_scale` should be either `unidir` or `bidir`")
	}
	mc_color = mc_color_fun(mc_vector)
	mc_color_labels = seq(from=mc_min, to=mc_max, length.out=9)
	mc_color_legend = mc_color_fun(mc_color_labels)
	mc_color_labels = sprintf("%.2f",mc_color_labels)
	mc_color_labels [ length(mc_color_labels) ] = paste(">=", mc_color_labels [ length(mc_color_labels) ], sep="")
	mc_color_labels [ 1 ] = paste("<=", mc_color_labels [ 1 ], sep="")
	# zero color
	if (!is.null(mc_zero_color)) {
		mc_color [ mc_vector == 0 ] = mc_zero_color
		mc_color_legend [ mc_color_labels == 0 ] = mc_zero_color
	}

	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		# determine plot max/min
		if (plot_mcs) {
			xlim=c(min(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE), max(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE), max(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE))
		} else {
			xlim=c(min(mc2d@sc_x, na.rm = TRUE), max(mc2d@sc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, na.rm = TRUE), max(mc2d@sc_y, na.rm = TRUE))
		}
		xlim[1] = xlim[1] - 0.2 * abs(xlim[1] - xlim[2])
		xlim[2] = xlim[2] + 0.2 * abs(xlim[1] - xlim[2])
		ylim[1] = ylim[1] - 0.2 * abs(ylim[1] - ylim[2])
		ylim[2] = ylim[2] + 0.2 * abs(ylim[1] - ylim[2])
		
		# get single cell coordinates
		mc2d_sc_x = mc2d@sc_x [ names(mc2d@sc_x) %in% names(sc_vector) ]
		mc2d_sc_y = mc2d@sc_y [ names(mc2d@sc_y) %in% names(sc_vector) ]
		sc_vector = sc_vector [ names(sc_vector) %in% names(mc2d_sc_x) ]
		
		# plot individual cells
		# by default, plot individual cells
		if (! plot_cells_as_2d_density ) {
			plot(
				mc2d_sc_x,
				mc2d_sc_y,
				pch = 19,lwd = 0,
				cex = cex_sc,
				col = alpha(sc_color, alpha_sc),
				xlim = xlim,
				ylim = ylim,
				xlab = NA, ylab = NA, axes=do_axes
			)
		
		# if else, plot 2d density (works with UMI counts or UMI frac, probably with anythin that's unidirectional too)
		} else {
			
			# create fake mc2d positions where each cell appears in the 
			# same positionas many times as its UMI counts
			den_mc2d_x = unlist(plyr::alply(1:length(sc_vector), 1, function(i) {
				if (sc_vector[i] > sc_min) { rep(mc2d_sc_x[i], sc_vector[i]) }
			}))
			den_mc2d_y = unlist(plyr::alply(1:length(sc_vector), 1, function(i) {
				if (sc_vector[i] > sc_min) { rep(mc2d_sc_y[i], sc_vector[i]) }
			}))
			
			# drop NA?
			den_mc2d_x = den_mc2d_x [ !is.na(den_mc2d_x) ]
			den_mc2d_y = den_mc2d_y [ !is.na(den_mc2d_y) ]
			
			# # add a bit of jitter
			# den_mc2d_x = jitter(den_mc2d_x, factor = 1)
			# den_mc2d_y = jitter(den_mc2d_y, factor = 1)
			
			# create color scale
			den_mc2d_colvec = unidir_color_scale
			den_mc2d_colvec[1] = NA
			den_mc2d_colfun = colorRampPalette(den_mc2d_colvec)
			den_mc2d_col = den_mc2d_colfun(40)
			
			# create interploated 2d density plot
			den_mc2d_k = MASS::kde2d(den_mc2d_x, den_mc2d_y, n = 200, lims = c(xlim, ylim))
			
			# plot
			image(den_mc2d_k, col= alpha(den_mc2d_col, alpha_sc_2d), xlim = xlim, ylim = ylim, axes = do_axes, add = FALSE, useRaster = TRUE)
			points(
				mc2d@sc_x,
				mc2d@sc_y,
				pch = 1, lwd = 0.5,
				cex = cex_sc,
				col = alpha("black", alpha_sc),
				xlim = xlim,
				ylim = ylim,
				xlab = NA, ylab = NA, axes=do_axes
			)
			
		}
		
		# plot legend
		if (do_legend_sc) {
			legend("topleft",fill=sc_color_legend, legend=sc_color_labels, cex=0.4, title=sc_label)
		}
		
		# plot edges between metacells?
		if (plot_edges) {
			fr = mc2d@graph$mc1
			to = mc2d@graph$mc2
			segments(
				x0=mc2d@mc_x[fr],
				y0=mc2d@mc_y[fr],
				x1=mc2d@mc_x[to],
				y1=mc2d@mc_y[to],
				col="gray70")
		}
		
		# plot metacells?
		if (plot_mcs) {
		points(
			mc2d@mc_x,
			mc2d@mc_y,
			pch=19,lwd=0.5,
			cex=cex_mc,
			col=alpha(mc_color,alpha_mc))
			if (do_legend_mc) {
				legend("topright",fill=mc_color_legend, legend=mc_color_labels, cex=0.4, title=mc_label)
			}
		}
		
		# plot metacell ids?
		if (plot_mc_name) {
			text(
				mc2d@mc_x,
				mc2d@mc_y,
				#cex=cex_mc,
				labels=names(mc2d@mc_x))
		}
		
		# plot title
		if (is.null(title)) {
			title = gene_id
		}
		title(
			main = title,
			sub = sprintf(
				"2D projection\nn = %i cells | n = %i metacells", 
				length(mc2d@sc_x),
				length(mc2d@mc_x)
			)
		)
		
	}
	)

}


#' Binarise gene expression matrix
#'
#' @param counts a matrix with gene counts (rows) per cell/metacell/else (columns). 
#'    By default, the appropriate threshold is obtained from the first valley in the 
#'    counts distribution (excluding zero values).
#' @param quantile numeric, a quantile of expression used to obtain a threshold for 
#'    binarisation (default is NULL, i.e. threshold is obtained from distribution).
#' @param threshold numeric, a hard threshold for binarisation (default is NULL, 
#'    i.e. threshold is obtained from distribution).
#' @param apply_on_nonzero logical, whether to apply quantiles/thresholds on non-zero
#'    values, or the whole dataset (default is TRUE)
#' 
#' @return matrix with binarised counts
#' 
sca_binarise_expression = function(
	counts,
	quantile = NULL,
	threshold = NULL,
	apply_on_nonzero = TRUE
	
) {
	
	# get distribution of nonzero values
	if (apply_on_nonzero) {
		counts_vec = counts [ counts > 0 ]
	} else {
		counts_vec = counts
	}
	
	if (is.null(threshold) & is.null(quantile)) {
		# if no quantile is specified, get threshold from first
		# valley in the counts distribution
		counts_vec_hist = hist(log10(counts_vec), plot=FALSE)
		threshold_ix = find_peaks(x = 1 / (counts_vec_hist$density + 1), m = 2)[1]
		if (is.null(threshold_ix)) {
			threshold_ix = 1
		}
		threshold = 10 ^ (counts_vec_hist$breaks [ threshold_ix ])
		message(sprintf("Binarising counts at >= %.2f threshold (first valley)", threshold))
	} else if (is.null(threshold) & is.numeric(quantile)) {
		# else, use quantile to define threshold
		threshold = quantile(counts_vec, quantile)
		message(sprintf("Binarising counts at >= %.2f threshold (quantile %.2f)", threshold, quantile))
	} else if (is.numeric(threshold)) {
		message(sprintf("Binarising counts at >= %.2f threshold (given)", threshold))
	}
	
	# apply threshold
	counts_out = (counts >= threshold) * 1
	rownames(counts_out) = rownames(counts)
	colnames(counts_out) = colnames(counts)
	
	# output
	return(counts_out)
	
}


#' Find peaks in a vector of data
#'
#' @param x vector of numeric values
#' @param m numeric, number of data points at each side of the peak that
#'     need to be lower than the peak (default is 3)
#' 
#' @return indexes of all peaks
find_peaks = function (x, m = 3){
	shape = diff(sign(diff(x, na.pad = FALSE)))
	pks = sapply(which(shape < 0), FUN = function(i){
		z = i - m + 1
		z = ifelse(z > 0, z, 1)
		w = i + m + 1
		w = ifelse(w < length(x), w, length(x))
		if (all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) {
			return(i + 1) 
		} else {
			return(numeric(0))
		}
	})
	pks = unlist(pks)
	pks
}


#' Merge multiple metacell solution in a single object
#' 
#' @param mat_ids character, matrix ids in database to merge (they should br in the same order as corresponding `mc_ids`)
#' @param mc_ids character, metacell ids in database to merge (they should br in the same order as corresponding `mat_ids`)
#' @param new_mat_id character, name of the new matrix object to add to database;
#'    if NULL (default), the object is not added to metacell database
#' @param new_mc_id character, name of the new metacell object to add to database;
#'    if NULL (default), the object is not added to metacell database
#' 
#' @return metacell object (`gMCCov` class)
#' Cell names will have corresponding original matrix ids prepended to them, 
#' and metacell ids will have "mc_" and matrix id prepended to the original ids.
#' 
sca_merge_metacell_solutions <- function(
	mat_ids, 
	mc_ids, 
	new_mat_id = NULL, 
	new_mc_id = NULL,
	use_tgMCCov = TRUE
) {
	
	# merge UMI matrices
  message(Sys.time(), " Loading mats")
	mat_list <- lapply(mat_ids, function(x){
	  message(x)
	  mat <- scdb_mat(x)
	  # add original mat ids to cell names in UMI matrix
	  mat@cell_metadata$original_mat_id <- x
	  mat@cell_metadata$cell <- rownames(mat@cell_metadata)
	  rownames(mat@cell_metadata) <- paste(x, rownames(mat@cell_metadata), sep="_")
	  colnames(mat@mat) <- rownames(mat@cell_metadata)[match(colnames(mat@mat), mat@cell_metadata$cell)]
	  mat@cells <- rownames(mat@cell_metadata)[match(mat@cells, mat@cell_metadata$cell)]
	  mat
	})
	message(Sys.time(), " Merging mats")
	
	mat <- sca_merge_mats(mat_list)

	# list cells in metacells and metacell colors
	message(Sys.time(), " Loading mcs")
	mc_list <- lapply(1:length(mc_ids), function(i) {
	  x <- mc_ids[i]
	  y <- mat_ids[i]
	  message(x)
		mc <- scdb_mc(x)
		mc_vector <- sort(mc@mc)
		if (is.null(names(mc@colors)) | length(names(mc@colors)) == 0)
			names(mc@colors) <- seq_along(mc@colors)
		list(
			mc_vector = structure(paste(y, mc_vector, sep="_mc"), names=paste(y, names(mc_vector), sep="_")),
			mc_cols = structure(mc@colors, names=paste(y, names(mc@colors), sep="_mc"))
		)
	})
	names(mc_list) <- NULL
	sc_mc_vector <- unlist(lapply(mc_list, function(x) x[[1]]), use.names = TRUE)
	mc_cols <-  unlist(lapply(mc_list, function(x) x[[2]]), use.names = TRUE)
	if (any(is.na(sc_mc_vector))) {
		sc_mc_vector <- sc_mc_vector[!is.na(sc_mc_vector)]
		mc_cols <- mc_cols[!is.na(sc_mc_vector)]
	}
	
	# update matrix to match mc
	message(Sys.time(), " Updating mat to match mc")
	mat@mat <- mat@mat[,names(sc_mc_vector)]
	mat@cells <- colnames(mat@mat)
	mat@ncells <- length(mat@cells)
	mat@cell_metadata <- mat@cell_metadata[mat@cells,]
	
	# add new matrix to database
	if (!is.null(new_mat_id)) {
		scdb_add_mat(new_mat_id, mat)
		message(Sys.time(), " Added new matrix object (", new_mat_id, ") to database")
	}
	
	# create mc object
	message(Sys.time(), " Creating new metacell object")
	sc_mc_vector_int <- as.integer(as.factor(sc_mc_vector))
	names(sc_mc_vector_int) <- names(sc_mc_vector)
	if (use_tgMCCov == TRUE) {
	  mc <- new("tgMCCov", mc = sc_mc_vector_int, outliers = vector("character"), scmat = mat)
	  mc@mc <- sc_mc_vector_int
	} else {
	  # calculate mc footprint manually
	  ct_geomean=t(apply(mat, 1,  function(x) tapply(x, sc_mc_vector_int, function(y) exp(mean(log(1+y)))-1)))
	  ct_meansize=tapply(colSums(mat), sc_mc_vector_int, mean)
	  ideal_cell_size=pmin(1000,median(ct_meansize))
	  g_fp=t(ideal_cell_size * t(ct_geomean)/as.vector(ct_meansize))
	  fp_reg=0.05
	  g_fp_n=(fp_reg + g_fp)/apply(fp_reg + g_fp, 1, median)
	  # create mc objrect manually
	  mc=scdb_mc(mc_ids[1])
	  mc@mc_fp=g_fp_n
	  mc@mc=sc_mc_vector_int
	  mc@cell_names=names(sc_mc_vector_int)
	}
	
	# add new mc to databasr
	if (!is.null(new_mc_id)){
		scdb_add_mc(new_mc_id, mc)
		message(Sys.time(), " Added new metacell object (", new_mc_id, ") to database")
	}
	
	return(mc)
}

#' Single cell-level identification of batch-correlated genes
#'
#' @param mat_object `mat` object with expression data (`mat@mat`) and batch information 
#' @param cor_thr_q quantile of top correlation values to use as a threshold (default = 0.99, ie top 1% of genes are excluded)
#' @param cor_thr correlation value to use as a threshold (default is NULL, ie use the quantile `cor_thr_q` value instead)
#' @param bidirectional default FALSE; whether to remove anti correlated genes as well
#' @param method correlation method (from `cor()`, default is "pearson")
#' @param batch_field character, name of the column in `mat_object@cell_metadata` containing batch information (default is "batch_set_id"), or a custom vector (must match cell order in `mat`!)
#' 
#' @return metacell object (`gMCCov` class)
#'
sca_batch_correlated_genes_sc = function(mat_object, cor_thr_q = 0.99, cor_thr = NULL, bidirectional = FALSE, method = "pearson", batch_field = "batch_set_id") {
	
	# expression and batch data per cell
	mat_e = as.matrix(mat_object@mat)
	if (length(batch_field) == 1) {
		mat_b = as.numeric(factor(mat_object@cell_metadata[colnames(mat_e), batch_field ]))
	} else if (length(batch_field) == length(mc_object@mc)) {
		mat_b = as.numeric(factor(batch_field))
	}
	
	# correlation
	gen_bat_cor = apply(mat_e, 1, function(g) cor(g, mat_b, method = method) )
	
	# get anticorrelated genes as well?
	if (bidirectional) {
		gen_bat_cor = abs(gen_bat_cor)
	}
	
	# find genes
	if (!is.null(cor_thr)) {
		gen_bat_ids = names(which(gen_bat_cor >= cor_thr))
		message(sprintf("Batch-correlated genes | n = %i | %s cor >= %.3f | quantile = [given]", length(gen_bat_ids), method, cor_thr ))
	} else {
		cor_thr = quantile(gen_bat_cor, cor_thr_q, na.rm = TRUE)
		gen_bat_ids = names(which(gen_bat_cor >= cor_thr))
		message(sprintf("Batch-correlated genes | n = %i | %s cor >= %.3f | quantile = %.3f", length(gen_bat_ids), method, cor_thr, cor_thr_q ))
	}
	
	# hist(gen_bat_cor, breaks = 60, xlab = method, main = "Gene UMI ~ batch correlation", border = NA)
	
	# return	
	return(list(genes = gen_bat_ids, cor = gen_bat_cor))
	
}


#' Metacell-level identification of batch-correlated genes
#'
#' @param mc_object `mc` object with expression data (`mat@mc_fp`)
#' @param mat_object `mat` object with batch information 
#' @param cor_thr_q quantile of top correlation values to use as a threshold (default = 0.99, ie top 1% of genes are excluded)
#' @param cor_thr correlation value to use as a threshold (default is NULL, ie use the quantile `cor_thr_q` value instead)
#' @param method correlation method (from `cor()`, default is "pearson")
#' @param batch_field character, name of the column in `mat_object@cell_metadata` containing batch information (default is "batch_set_id"), or a custom vector (must match cell order in `mat`!)
#' 
#' @return list with two elements: `genes` correlated with batch distribution, and  `cor`, a correlation matrix
#'
sca_batch_correlated_genes_mc = function(
	mc_object, mat_object, 
	cor_thr_q = 0.999, 
	cor_thr_p = NULL,
	cor_thr = NULL, 
	method = "pearson",
	batch_field = "batch_set_id",
	do_plots = FALSE,
	output_file = NULL, width = 4, height = 4, res = NA
	) {
	
	# get batch frequencies per metacell
	if (length(batch_field) == 1) {
		mc_t = t(table(mc_object@mc, mat_object@cell_metadata [ names(mc_object@mc) , batch_field]))
	} else if (length(batch_field) == length(mc_object@mc)) {
		mc_t = t(table(mc_object@mc, batch_field))
	} else {
		stop("`batch_field` should either be a column in `mat@cell_metadata` or a vector with cell-level batch information, matching `mc@mc`")
	}
	# keep metacell order from mc object
	mc_t = mc_t [ , colnames(mc_object@mc_fp) ]
	mc_f = t(t(mc_t) / colSums(mc_t))

	# expression and batch data per cell
	mc_e = as.matrix(mc_object@mc_fp)
	
	# correlation
	gen_bat_cor = t(apply(mc_e, 1, function(g) cor(g, t(mc_f), method = method) ))
	colnames(gen_bat_cor) = rownames(mc_f)

	# max values per gene
	gen_bat_cor_max = apply(abs(gen_bat_cor), 1, max)
	
	# find genes
	if (!is.null(cor_thr)) {
		gen_bat_ids = names(which(gen_bat_cor_max >= cor_thr))
		message(sprintf("Batch-correlated genes | n = %i | %s cor >= %.3f | quantile = [given]", length(gen_bat_ids), method, cor_thr ))
	} else if (!is.null(cor_thr_p)) {
		gen_bat_cor_max_p = as.vector(pnorm(scale(gen_bat_cor_max)))
		names(gen_bat_cor_max_p) = names(gen_bat_cor_max)
		cor_thr = min(gen_bat_cor_max [ gen_bat_cor_max_p >= cor_thr_p ])
		gen_bat_ids = names(which(gen_bat_cor_max >= cor_thr))
		message(sprintf("Batch-correlated genes | n = %i | %s cor >= %.3f | probability = %.3f", length(gen_bat_ids), method, cor_thr, cor_thr_p ))
	} else {
		cor_thr = quantile(gen_bat_cor_max, cor_thr_q, na.rm = TRUE)
		gen_bat_ids = names(which(gen_bat_cor_max >= cor_thr))
		message(sprintf("Batch-correlated genes | n = %i | %s cor >= %.3f | quantile = %.3f", length(gen_bat_ids), method, cor_thr, cor_thr_q ))
	}
	
	
	# diagnostic plots
	if (do_plots) {
		plotting_function(output_file, width, height, res, EXP = {
			
			par(mfrow = c(2,2))
			# histogram
			hist(
				gen_bat_cor, breaks = 60, xlab = method,
				main = "distribution of fp ~ batch correlation values",
				border = NA, xlim = c(-1,1))
			
			# max correlation per gene
			plot(
				sort(gen_bat_cor_max), col = "blue", cex = 0.5,
				ylab = "Max correlation per gene",
				sub = sprintf("thr = %.3f | n = %i genes", cor_thr, length(gen_bat_ids)),
				main = "Max correlation per gene"
			)
			abline(h=cor_thr, lty = 2)
			
			# pairwise scatterplots
			for (i in 1:ncol(gen_bat_cor)) {
				for (j in 1:ncol(gen_bat_cor)) {
					if (i < j) {
						p_cor = cor(gen_bat_cor[,i], gen_bat_cor[,j], method = method)
						plot(
							gen_bat_cor[,i], gen_bat_cor[,j], 
							xlab = colnames(gen_bat_cor)[i],
							ylab = colnames(gen_bat_cor)[j], 
							cex = 0.5, col = alpha("blue",0.6),
							xlim = c(-1,1), ylim = c(-1,1),
							main = sprintf("Per-gene batch correlation in batch pairs\n%s (%i) ~ %s (%i)", colnames(gen_bat_cor)[i], i, colnames(gen_bat_cor)[j], j),
							sub = sprintf("%s = %.3f", method, p_cor),
							cex.main = 0.9)
					}
				}
			}
		})
	} 
	# end diagnostic plots
	
	# return	
	return(list(genes = gen_bat_ids, cor = gen_bat_cor))
	
}


#' Identify low-quality metacells to drop based on UMI distribution
#'
#' @param mc_object `mc` object with expression data (`mat@mc_fp`)
#' @param mat_object `mat` object with batch information 
#' @param p_thr_median,p_thr_total probability threshold below which cells with very few median or total UMIs (respectively) are dropped (after scaling the log-transformed metacell-level total or median UMI distribution).
#' @param do_global_filter if TRUE (default), apply filters to global distribution (ie all metacells at the same time). All criteria are additive.
#' @param mc_clusters a vector indicating clusters of metacells (same order as columns in mc@mc_fp). Default is NULL. If specified, thresholds are applied to cluster-wise scaled distributions, which can identify cells within a cluster that have very few UMIs. All criteria are additive.
#' 
#' @return vector with metacell names to drop.
#'
sca_drop_mcs_by_umis = function(
	mc_object, mat_object, mc_color = NULL,
	do_global_filter = TRUE, 
	mc_clusters = NULL,
	p_thr_median = 0.01, 
	p_thr_total = 0.01, 
	mc_counts = NULL, 
	min_median_umis = 0, 
	min_total_umis = 0,
	do_plots = FALSE,
	output_file = NULL, width = 6, height = 6, res = NA
	) {

	# init vector
	drop_mcs = NULL
	drop_mcs_cw = NULL
	
	# colors for mcs?
	if (is.null(mc_color)) {
		mc_color = rep("blue", ncol(mc_object@mc_fp))
	}

	# common stats
	if (is.null(mc_counts)) {
		mc_counts = sca_mc_gene_counts(mc_object, mat_object, 0)
	}
	sc_counts = Matrix::colSums(mat@mat[,names(mc_object@mc)])
	
	# global filters
	if (do_global_filter) {
		
		# discard metacells based on very few total UMIs
		mc_counts_zscore = scale(log10(colSums(mc_counts)))[,1]
		drop_mcs_i = union( names(which(pnorm(mc_counts_zscore) < p_thr_total)) , names(which(colSums(mc_counts) < min_total_umis)) )
		message(sprintf("drop metacells | mcs to drop due to low total UMIs  = %s (n=%i)", paste(sort(as.numeric(drop_mcs_i)), collapse = ","), length(drop_mcs_i)))
		drop_mcs = union(drop_mcs, drop_mcs_i)
		
		# discard metacells based on very few median UMIs per cell
		mc_median_cell_size = tapply(sc_counts, mc_object@mc, median)
		mc_median_cell_size_zscore = scale(log10(mc_median_cell_size))[,1]
		drop_mcs_i = union( names(which(pnorm(mc_median_cell_size_zscore) < p_thr_median)) , names(which( mc_median_cell_size < min_median_umis )))
		message(sprintf("drop metacells | mcs to drop due to low median UMIs = %s (n=%i)", paste(sort(as.numeric(drop_mcs_i)), collapse = ","), length(drop_mcs_i)))
		drop_mcs = union(drop_mcs, drop_mcs_i)
		
	}

	# cluster-wise filters
	if ( !is.null(mc_clusters) ) {

		# discard metacells based on very few UMIs compared to other similar metacells
		drop_mcs_io_t = NULL
		drop_mcs_io_m = NULL
		for (cti in unique(mc_clusters)) {

			cti_ixs = which(mc_clusters == cti)
			# total counts
			mc_counts_i = colSums(mc_counts)[ cti_ixs ]
			mc_counts_i_zscore = scale(log10(mc_counts_i))[,1]
			drop_mcs_i_t = union( names(which(pnorm(mc_counts_i_zscore) < p_thr_total)) , names(which(mc_counts_i < min_total_umis)) )
			# median counts
			mc_median_cell_size_i = tapply(sc_counts [ mc_object@mc %in% cti_ixs ], mc_object@mc [ mc_object@mc %in% cti_ixs ], median)
			mc_median_cell_size_i_zscore = scale(log10(mc_median_cell_size_i))[,1]
			drop_mcs_i_m = union( names(which(pnorm(mc_median_cell_size_i_zscore) < p_thr_median)),  names(which( mc_median_cell_size_i < min_median_umis ))  )
			# concatenate
			drop_mcs_i = union(drop_mcs_i_t, drop_mcs_i_m)
			drop_mcs_io_t = c(drop_mcs_io_t, drop_mcs_i_t)
			drop_mcs_io_m = c(drop_mcs_io_m, drop_mcs_i_m)
			
		}
		message(sprintf("drop metacells | mcs to drop due to low cluster-wise total UMIs  = %s (n=%i)", paste(sort(as.numeric(drop_mcs_io_t)), collapse = ","), length(drop_mcs_io_t)))
		message(sprintf("drop metacells | mcs to drop due to low cluster-wise median UMIs = %s (n=%i)", paste(sort(as.numeric(drop_mcs_io_m)), collapse = ","), length(drop_mcs_io_m)))
		drop_mcs_cw = union(drop_mcs_io_m, drop_mcs_io_t)
		drop_mcs = union(drop_mcs, drop_mcs_cw)

	}
	

	# log
	message(sprintf("drop metacells by UMI | mcs to drop = %s (n=%i)", paste(sort(as.numeric(drop_mcs)), collapse = ","), length(drop_mcs)))
	message(sprintf("drop metacells by UMI | total cell to drop = %i out of %i (%.2fp)", sum(mc_object@mc %in% drop_mcs),          length(mc_object@mc),  100 * sum(mc_object@mc %in% drop_mcs) / length(mc_object@mc) ))
	message(sprintf("drop metacells by UMI | total UMIs to drop = %i out of %i (%.2fp)", sum(colSums(mc_counts)[drop_mcs]), sum(mc_counts), 100 * sum(colSums(mc_counts)[drop_mcs]) / sum(mc_counts) ))
	gc()
	
	# do diagnostic plots?
	if (do_plots & length(drop_mcs) > 0) {
		plotting_function(output_file, width, height, res, EXP = {
		plot(
			colSums(mc_counts),
			mc_median_cell_size,
			col = mc_color, 
			cex = ifelse(names(mc_median_cell_size) %in% drop_mcs, 1, 0.7),
			pch = ifelse(names(mc_median_cell_size) %in% drop_mcs_cw, 17, ifelse(names(mc_median_cell_size) %in% drop_mcs, 19, 20)),
			log = "xy",
			ylim = c(min(mc_median_cell_size),max(mc_median_cell_size)),
			xlim = c(min(colSums(mc_counts)),max(colSums(mc_counts))),
			ylab = "median cell size in metacell", xlab = "metacell size",
			main = sprintf("Discarded cells"))
		legend("topleft", c("kept", "discarded (global dist)", "discarded (cluster-wise dists)"), pch = c(20,19,17), col = "black", bty = "n")
		text(
			colSums(mc_counts) [ as.numeric(drop_mcs) ], 
			mc_median_cell_size [ as.numeric(drop_mcs) ], 
			names(mc_median_cell_size) [ as.numeric(drop_mcs) ],
			col = "black")
		})
	}

	
	# return
	return(drop_mcs)

}


#' Identify metacells to drop based on fold change distribution
#'
#' @param mc_object `mc` object with expression data (`mat@mc_fp`)
#' @param min_fc,min_markers thresholds_ find at least n=`min_markers` genes with this fp<`min_fc`
#' @param key_markers use only gene markers from this list (e.g. TFs, to filter out metacells without any specific TF).
#' 
#' @return vector with metacell names to drop.
#'
sca_drop_mcs_by_fc = function(mc_object, min_fc = 5, min_markers = 10, key_markers = NULL) {

	mc_fp = mc_object@mc_fp
	if (!is.null(key_markers)) {
		mc_fp = mc_object@mc_fp [ key_markers, ]
		message(sprintf("drop metacells by FC | focus on n=%i markers", length(key_markers)))
	}

	drop_mcs = names(which(apply(mc_fp, 2, function(i) length(which(sort(i, decreasing = TRUE) > min_fc)) < min_markers )))
	message(sprintf("drop metacells by FC | mcs to drop due to <%i markers with fc<%.2f = %s (n=%i)", min_markers, min_fc, paste(sort(as.numeric(drop_mcs)), collapse = ","), length(drop_mcs)))
	message(sprintf("drop metacells by FC | total cell to drop = %i out of %i (%.2fp)", sum(mc_object@mc %in% drop_mcs),          length(mc_object@mc),  100 * sum(mc_object@mc %in% drop_mcs) / length(mc_object@mc) ))
	# message(sprintf("drop metacells by FC | total UMIs to drop = %i out of %i (%.2fp)", sum(colSums(mc_counts)[drop_mcs]), sum(mc_counts), 100 * sum(colSums(mc_counts)[drop_mcs]) / sum(mc_counts) ))

	return(drop_mcs)

}

#' Identify doublet cells based on the `scDblFinder::findDoubletClusters` method, with parallelisation
#'
#' @param mat, UMI matrix (e.g. `mat@mat`)
#' @param clusters, vector with clusters to which each gene belongs (cell types, metacells...), same order as columns in `mat`
#' @param threshold, p-value threshold to define DE genes between clusters, default is 0.01
#' @param num_threads, num threads for doublet detection loop, default is 4
#' @param BPPARAM, `BiocParallel` method for parallelisation of `scran::findMarkers` (slowest step)
#' 
#' @return dataframe with putative donor clusters for each query cluster. See https://rdrr.io/github/plger/scDblFinder/man/findDoubletClusters.html
#'
sca_find_doublet_clusters = function(
	mat, 
	clusters, 
	subset_genes = NULL, 
	threshold = 0.01, 
	num_threads = 4,
	BPPARAM = BiocParallel::MulticoreParam(4, progressbar = TRUE, log = FALSE)
	) {
	if (length(unique(clusters)) < 3L) {
		stop("need at least three clusters to detect doublet clusters")
	}

	# Computing normalized counts using the library size (looking for compositional differences!)
	message("doublet detection | create sce...")
	sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
	message("doublet detection | size factors...")
	SingleCellExperiment::sizeFactors(sce) = scuttle::librarySizeFactors(mat, subset_row = subset_genes)
	message("doublet detection | log-normalise...")
	sce <- scuttle::logNormCounts(sce)

	# find cluster-specific markers (can take a long time if many clusters)
	message(sprintf("doublet detection | find markers for n=%i clusters, with n=%i threads...", length(unique(clusters)), num_threads))
	degs <- scran::findMarkers(sce, clusters, subset.row = subset_genes, full.stats = TRUE, BPPARAM = BPPARAM)
	message("doublet detection | split...")
	med.lib.size <- vapply(split(SingleCellExperiment::sizeFactors(sce), clusters), FUN=median, FUN.VALUE=0)
	n.cluster <- table(clusters)/length(clusters)

	# Setting up the output.
	all.clusters <- names(degs)
	collected.top <- collected.all <- vector("list", length(all.clusters))
	names(collected.top) <- names(collected.all) <- all.clusters

	# Running through all pairs of clusters and testing against the third cluster.
	fun_loop_ref = function(rn) {

        # get ref
		ref = all.clusters[rn]
		ref.stats <- degs[[ref]]
		remnants <- setdiff(all.clusters, ref)
		
        # log progress
		if (rn == 1 | rn %% 10 == 0) {
			message(sprintf("doublet detection | score cluster %i/%i...", rn, length(all.clusters)))
		}

        # work on ref
		num <- length(remnants) * (length(remnants) - 1L)/2L
		all.N <- med.N <- all.gene <- all.parent1 <- all.parent2 <- integer(num)
		all.p <- numeric(num)
		idx <- 1L

		for (i1 in seq_along(remnants)) {
			stats1 <- ref.stats[[paste0("stats.", remnants[i1])]]
			for (i2 in seq_len(i1-1L)) {

				stats2 <- ref.stats[[paste0("stats.", remnants[i2])]]

				# Obtaining the IUT and setting opposing log-fold changes to 1.
				max.log.p <- pmax(stats1$log.p.value, stats2$log.p.value)
				max.log.p[sign(stats1$logFC) != sign(stats2$logFC)] <- 0

				# Correcting pval across genes
				log.adj.p <- scran:::.logBH(max.log.p)
				best.gene <- which.min(max.log.p)[1]

				all.N[idx] <- sum(log.adj.p <= log(threshold), na.rm=TRUE)
				all.gene[idx] <- best.gene
				all.p[idx] <- exp(log.adj.p[best.gene])
				all.parent1[idx] <- i1
				all.parent2[idx] <- i2
				idx <- idx + 1L
			}
		}

		# Formatting the output.
		parent1 <- remnants[all.parent1]
		parent2 <- remnants[all.parent2]

		stats <- data.frame(
			source1=parent1,
			source2=parent2,
			num.de=all.N,
			median.de=rep(0, length(all.N)), # placeholder, see below.
			best=rownames(ref.stats)[all.gene],
			p.value=all.p,
			lib.size1=unname(med.lib.size[parent1]/med.lib.size[ref]),
			lib.size2=unname(med.lib.size[parent2]/med.lib.size[ref])
		)

		o <- order(all.N, -all.p)
		top <- cbind(stats[o[1],], prop=n.cluster[[ref]])
		med.de <- median(all.N)
		top$median.de <- med.de
		rownames(top) <- ref
		collected.top[[ref]] <- top

		return(collected.top)
	}

	# loop over query cell types
	doParallel::registerDoParallel(cores = num_threads)
	outs = plyr::alply(.data = 1:length(all.clusters), .margins = 1,  .fun = fun_loop_ref, .parallel = TRUE, .inform = TRUE)
	names(outs) = all.clusters
	message(sprintf("doublet detection | score cluster %i/%i...", length(all.clusters), length(all.clusters)))

	# Returning the DataFrame of compiled results.
	message(sprintf("doublet detection | get top pair per source cluster..."))
	out <- do.call(rbind, outs)
	out = as.data.frame(t(sapply(1:length(outs), function(v)( unlist(outs[[v]])))))
	colnames(out) = c("source1","source2","num.de","median.de","best","p.value","lib.size1","lib.size2","prop")
	for (cc in c("num.de","median.de","p.value","lib.size1","lib.size2","prop")) {
		out[,cc] = as.numeric(out[,cc])
	}
	out$query = all.clusters
	rownames(out) = out$query
	out = out [ , c("query","source1","source2","num.de","median.de","best","p.value","lib.size1","lib.size2","prop") ]
	out = out[order(out$num.de, -out$prop, c(out$lib.size1+out$lib.size2)),]
	return(out)
	message(sprintf("doublet detection | done!"))
	
}


#' Identify doublet cells based on the `scDblFinder::findDoubletClusters` method, with parallelisation
#'
#' @param mc_counts, mc counts
#' @param mc_annotation, mc annotations (including color column)
#' @param recluster_mcs, whether to recluster metacells into louvain clusters (default FALSE, UNIMPLEMENTED)
#' @param louvain_resolution, resolution for louvain clustering of metacells (UNIMPLEMENTED)
#' @param cluster_algorithm, algorithm for clustering of metacells (UNIMPLEMENTED)
#' @param k_weight_sct, k weight for SCT-transformation in Seurat (default 50)
#' @param num_pcs, num PCs to use for NN calculations
#' @param num_pcs_for_umap, num PCs to retain for UMAP (if UMAP is used for NN)
#' @param build_tree_on_reduction, which reduction to use for tree building; either "pca" (default) or "umap"
#' @param p_val_thr, p-value threshold to determine differentiall expressed genes at each tree node
#' @param p_val_field, whether to use default or adjusted p-values ("p_val" or "p_val_adjusted")
#' @param min_num_de_genes, num DE genes in either direction (pos/neg) for an internal node to be considered valid (default is 10)
#' @param min_num_de_genes_positive, num overexpressed genes (pos DE) for an internal node to be considered valid (default is 0, i.e. ignore this parameter)
#' 
#' @return dataframe with putative donor clusters for each query cluster. See https://rdrr.io/github/plger/scDblFinder/man/findDoubletClusters.html
#'
sca_cluster_metacells_by_nn = function(
	mc_counts,
	mc_annotation,
	gene_annotations = NULL,
	background_mcs = "rest",
	recluster_mcs = FALSE,
	louvain_resolution = 3,
	cluster_algorithm = 1,
	k_weight_sct = 50,
	num_pcs = 50,
	num_pcs_for_umap = 30,
	build_tree_on_reduction = "pca",
	p_val_thr = 0.001,
	p_val_field = "p_val",
	min_num_de_genes = 10,
	min_num_de_genes_positive = 0,
	output_file = NULL
	) {
	
	# create seurat object
	message(sprintf("cluster mcs | create Seurat..."))
	if (is.null(mc_annotation)) {
		mc_annotation = data.frame(metacell = colnames(mc_counts), color = viridisLite::turbo(n = ncol(mc_counts), begin = 0.1, end = 0.9))
	}
	seu = Seurat::CreateSeuratObject(mc_counts, meta.data = mc_annotation)
	

	# dictionary to reconstruct gene names
	gg_v = rownames(mc_counts)
	names(gg_v) = gsub("_","-", rownames(mc_counts))
	
	# sct-transform to find markers
	message(sprintf("cluster mcs | transform and find markers..."))
	seu = Seurat::SCTransform(
		seu,
		do.scale = TRUE, 
		do.center = TRUE,
		k.weight = k_weight_sct,
		return.only.var.genes = FALSE,
		verbose = FALSE
	)
	# mat_seu = Seurat::PrepSCTFindMarkers(mat_seu, verbose = FALSE)
	
	# dimensionality reduction in mc space
	if (ncol(seu) < num_pcs) {
		message(sprintf("cluster mcs | decrease PC number to match number of cells..."))
		num_pcs = ncol(seu) - 1
		num_pcs_for_umap = ncol(seu) - 1
	}
	message(sprintf("cluster mcs | run PCA, n=%i PCs...", num_pcs))
	seu = Seurat::RunPCA(seu, npcs = num_pcs, verbose = FALSE)

	# if (dge_on == "counts") {
	# 	if (!any(is.null(mat_object) & is.null(mc_object))) {
	# 		stop("You need to define a `mat_object` and a `mc_object`")
	# 	} else {
	# 		mc_annotation_for_scs = mat_object@cell_metadata [ names(mc_object@mc), ]
	# 		mc_annotation_for_scs$metacell = mc@mc [rownames(mc_annotation_for_scs)]
	# 		sem = Seurat::CreateSeuratObject(as.matrix(mat_object@mat[,names(mc_object@mc)]), meta.data = mc_annotation_for_scs)
	# 		SeuratObject::Idents(sem) = sem@meta.data[,"metacell"]
	# 	}
	# }
	
	if (build_tree_on_reduction == "umap") {
		message(sprintf("cluster mcs | run UMAP, use n=%i PCs...", num_pcs_for_umap))
		seu = Seurat::RunUMAP(seu, dims = 1:num_pcs_for_umap, reduction = "pca", reduction.name = "umap", verbose = FALSE)
		num_pcs_for_tree = 2
	} else {
		num_pcs_for_tree = num_pcs_for_umap
	}

	# nearest neighbours graph for mcs
	message(sprintf("cluster mcs | find nearest neigbours..."))
	seu = Seurat::FindNeighbors(seu, dims = 1:num_pcs_for_tree, verbose = FALSE, reduction = build_tree_on_reduction, return.neighbor = FALSE)

	# find clusters with louvain
	if (recluster_mcs) {
		message(sprintf("cluster mcs | find clusters at resolution = %.1f...", louvain_resolution))
		seu = Seurat::FindClusters(seu, resolution = louvain_resolution, cluster.name = "cluster", verbose = FALSE, algorithm = cluster_algorithm)
		cluster_for_tree = "cluster"
	} else {
		cluster_for_tree = "metacell"
	}
	
	# tree of clusters
	message(sprintf("cluster mcs | build tree on %s(s)...", cluster_for_tree))
	SeuratObject::Idents(seu) = seu@meta.data[,cluster_for_tree]
	seu = Seurat::BuildClusterTree(seu, assay = "SCT", reduction = build_tree_on_reduction, dims = 1:num_pcs_for_tree, verbose = TRUE)
	phy = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	phy = ape::makeNodeLabel(phy)
	
	
	# find markers for internal nodes
	message(sprintf("cluster mcs | traverse %s tree: get ancestors of n=%i nodes...", cluster_for_tree, phy$Nnode - 1))
	descendants_from_nodes = adephylo::listTips(phy)
	list_labels = c(phy$tip.label, phy$node.label)
	list_labels_except_root = phy$node.label[!phy$node.label%in%phy$node.label[1]]
	
	# loop through tree: init table
	message(sprintf("cluster mcs | traverse %s tree: init bipartition markers table...", cluster_for_tree))
	markers_per_bipartition = data.frame()
	significant_markers_per_node = data.frame(node = phy$node.label, positive = 0, negative = 0, ancestor = "", row.names = phy$node.label)
	
	# if (dge_on == "counts") {
	# 	set = sem
	# 	layer = "counts"
	# } else {
	# 	set = seu
	# 	layer = "SCT"
	# }
	
	# loop through tree: loop
	for (focus_node in list_labels_except_root) {
		
		# get sister and ancestral nodes
		focus_node_ix = which(list_labels == focus_node)
		if(focus_node %in% phy$tip.label) {
			focus_tips = focus_node
		} else {
			focus_tips = as.character(sort(as.numeric(names(descendants_from_nodes[[focus_node]]))))
		}
		ancestor_node = list_labels[phangorn::Ancestors(phy, focus_node_ix)][1]
		
		# find markers in two ways: ingroup v rest or ingroup v sister
		if (background_mcs == "rest") {
			sister_tips = "rest"
			markers_focus_node = tryCatch(
				Seurat::FindMarkers(seu, ident.1 = focus_tips, layer = layer, min.cells.group = 1),
				error = function(e) {
					markers_focus_node = data.frame(p_val = 1, avg_log2FC = 0, pct.1 = 0, pct.2 = 0, p_val_adj = 1)
				}
			)
		} else if (background_mcs == "sister") {
			sister_tips = as.character(sort(as.numeric(names(descendants_from_nodes[[ancestor_node]]))))
			sister_tips = sister_tips [ ! sister_tips %in% focus_tips ]
			markers_focus_node = tryCatch(
				Seurat::FindMarkers(seu, ident.1 = focus_tips, ident.2 = sister_tips, layer = layer, min.cells.group = 1),
				error = function(e) {
					markers_focus_node = data.frame(p_val = 1, avg_log2FC = 0, pct.1 = 0, pct.2 = 0, p_val_adj = 1)
				}
			)
		}
		markers_focus_node$focus_node = focus_node
		markers_focus_node$ident.1.tips  = paste(focus_tips, collapse = ",")
		markers_focus_node$ident.2.tips = paste(sister_tips, collapse = ",")
		
		# keep significant only
		markers_focus_node = markers_focus_node [ markers_focus_node[,p_val_field] < p_val_thr, ]
		
		# count significant
		num_sigs = table(factor(sign(markers_focus_node$avg_log2FC), levels = c(1,-1)))
		significant_markers_per_node[focus_node,"positive"] = num_sigs[1]
		significant_markers_per_node[focus_node,"negative"] = num_sigs[2]
		significant_markers_per_node[focus_node,"ancestor"] = ancestor_node
		message(sprintf("cluster mcs | traverse %s tree: at %s node, n=%i markers (%i/%i pos/neg)", cluster_for_tree, focus_node, sum(num_sigs), num_sigs[1], num_sigs[2]))
		
		# store for later
		markers_focus_node$gene = rownames(markers_focus_node)
		markers_per_bipartition = rbind(markers_per_bipartition, markers_focus_node)
		
	}
	markers_per_bipartition$focus_node = factor(markers_per_bipartition$focus_node, levels = list_labels_except_root)
	
	# recompose gene names
	markers_per_bipartition$gene = as.character(gg_v [ markers_per_bipartition$gene ])
	rownames(markers_per_bipartition) = NULL
	
	# add gene annotations
	if (!is.null(gene_annotations)) {
		gene_annotations_v = gene_annotations[,2]
		names(gene_annotations_v) = gene_annotations[,1]
		markers_per_bipartition$gene_name = gene_annotations_v [ markers_per_bipartition$gene ]
	} else { 
		markers_per_bipartition$gene_name = ""
	}
	
	# plot tree of clusters
	width = ncol(seu) / 60 + 5
	height = ncol(seu) / 10 + 5
	plotting_function(output_file, width = width, height = height, res = NA, EXP = {
		ape::plot.phylo(phy, tip.color = mc_annotation$color, font = 2, edge.color = "darkgray", show.node.label = FALSE, cex = 0.7)
		node_labels = as.character(apply(significant_markers_per_node, 1, function(v) { sprintf("%s\n+%i/-%i", v[1], as.numeric(v[2]), as.numeric(v[3])) } ))
		node_labels_color = c("darkblue",scales::alpha("slateblue2",0.4)) [ factor(grepl("+0/-0",node_labels, fixed = TRUE), levels = c(FALSE, TRUE)) ]
		node_labels_colbg = c("lightblue",scales::alpha("white",0)) [ factor(grepl("+0/-0",node_labels, fixed = TRUE), levels = c(FALSE, TRUE)) ]
		ape::nodelabels(text = node_labels, bg = node_labels_colbg, cex = 0.7, col = node_labels_color, frame = "none")
		ape::add.scale.bar(y = -10, lcol = "darkgray")
		title(sub = sprintf("DGE genes in each node at %s < %.1E using %s as background", p_val_field, p_val_thr, background_mcs))
	})

	# filter good nodes
	valid_nodes = significant_markers_per_node [ 
		significant_markers_per_node$positive >= min_num_de_genes_positive &
		significant_markers_per_node$positive + significant_markers_per_node$negative >= min_num_de_genes,
	]
	
	# return	
	return(list(valid_nodes = valid_nodes, markers_per_bipartition = markers_per_bipartition, tree = phy))
	
}


#' Summary of UMI cell counts for metacell matrices.
#' 
#' This function summarizes UMI counts per cell and gene for metacell matrices.
#' It loads the matrices from the metacell database, applies optional filtering
#' based on UMI counts and gene counts, and returns a list of data.tables with
#' the summarized statistics.
#' 
#' @param scdb_dir metacell database directory
#' @param mat_ids list of metacell matrix ids; if NULL, all matrices are loaded
#' @param blacklist_genes list of genes to be excluded from the analysis
#' @param blacklist_genes_pattern list of gene patterns to be excluded from the analysis
#' @param gene_list list of gene sets to be included in the analysis; each element
#'  of the list should be a vector of gene names;
#' @param cell_UMI_thr threshold for filtering cells based on UMI counts; 
#'  If length 1, it is interpreted as minimum threshold above which to keep cells;
#'  if length 2, its elements are interpreted as lower and upper thresholds.
#' @param cell_genes_thr threshold for filtering cells based on gene counts
#'  If length 1, it is interpreted as minimum threshold above which to keep cells;
#'  if length 2, its elements are interpreted as lower and upper thresholds.
#' @param overwrite if TRUE, overwrite the original matrix with the optionally filtered
#'  matrix with updated metadata
#' 
#' @return list of data.tables with summarized stats:
#' \itemize{
#' \item{sum_cells} {cell UMI counts}
#' \item{sum_genes} {gene UMI counts}
#' \item{sum_median} {median UMI counts}
#' }
#' 
mc_summarize_counts <- function(
    scdb_dir, 
    mat_ids = NULL, 
    blacklist_genes = c(), 
    blacklist_genes_pattern = c(),
    gene_list = list(),
    cell_UMI_thr = 10,
    cell_genes_thr = 10,
    overwrite = FALSE
) {
  
    # initialize metacell database
    scdb_init(scdb_dir, force_reinit=TRUE)

    # load UMI matrices
    if (is.null(mat_ids)) {
        mat_ids <- str_remove(str_remove(
            list.files(scdb_dir, pattern = "mat.*", full.names = FALSE),
            "mat."
        ), ".Rda")
    }


    # loop over matrices
    sum_list <- sapply(mat_ids, function(mat_id) {

        # load UMI matrix
        message(mat_id)
        mat <- scdb_mat(mat_id)


        # blacklist genes
        genes <- rownames(mat@mat)
        if (length(blacklist_genes) > 0) {
            genes <- setdiff(genes, blacklist_genes)
        }
        if (length(blacklist_genes_pattern) > 0) {
            genes <- grep(paste(blacklist_genes_pattern, collapse = "|"), genes, invert = TRUE, value = TRUE)
        }
		if (any(length(blacklist_genes_pattern)>0, length(blacklist_genes) > 0)) 
			mat@mat <- mat@mat[genes, , drop = FALSE] 

        # count total UMIs per cell
        cell_sizes <- Matrix::colSums(mat@mat)

        # count UMIs for sets of genes
        cell_gene_sizes <- list()
        if (length(gene_list) > 0) {
          for (i in seq_along(gene_list)) {
            genes <- gene_list[[i]]
            genes <- intersect(genes, rownames(mat@mat))
            cell_gene_sizes[[i]] <- Matrix::colSums(mat@mat[genes, , drop = FALSE])
          }
		  if (!is.null(names(gene_list))) {
		    names(cell_gene_sizes) <- paste0(names(gene_list), "_UMI")
		  } else {
		    names(cell_gene_sizes) <- paste0("gene_list_", seq_along(cell_gene_sizes), "_UMI")
		  }
        }

        # merge UMI counts results
        if (!"dataset" %in% colnames(mat@cell_metadata)) mat@cell_metadata$dataset <- mat_id
        cell_smpls <- mat@cell_metadata[colnames(mat@mat), ]$dataset
        cell_batch <- mat@cell_metadata[colnames(mat@mat), ]$seq_batch_id
        cdt <- data.table(
            cell=colnames(mat@mat),
            cell_UMI=cell_sizes,
            seq_batch_id=cell_batch,
            dataset=cell_smpls
        )
        for (gene_set in names(cell_gene_sizes)) {
            cdt[[gene_set]] <- cell_gene_sizes[[gene_set]]
        }

        # count genes per cell
        cell_genes <- Matrix::colSums(mat@mat>0)
        gdt <- data.table(cell=colnames(mat@mat), cell_genes=cell_genes)
        cdt <- merge.data.table(cdt, gdt, by = "cell")
        setcolorder(cdt, c("cell", "seq_batch_id", "dataset", "cell_UMI", "cell_genes"))

        # count UMIs per gene
        gene_umis_dataset <- tapply(
            colnames(mat@mat), cell_smpls, function(x) 
                Matrix::rowSums(mat@mat[, x, drop = FALSE])
        )
        gene_umis_batch <- tapply(
          colnames(mat@mat), cell_batch, function(x)
            Matrix::rowSums(mat@mat[, x, drop = FALSE])
        )
        gdt_dataset <- rbindlist(
            lapply(names(gene_umis_dataset), function(x) data.table(
                gene=rownames(mat@mat), 
                gene_UMI=gene_umis_dataset[[x]],
                dataset=x
            ))
        )
        gdt_batch <- rbindlist(
            lapply(names(gene_umis_batch), function(x) data.table(
                gene=rownames(mat@mat), 
                gene_UMI=gene_umis_batch[[x]],
                seq_batch_id=x
            ))
        )
        gdt <- cbind(gdt_dataset, gdt_batch[, .(seq_batch_id)])

        # optional filtering
        if (!is.null(cell_UMI_thr)) {
            if (length(cell_UMI_thr) == 1) {
                cdt <- cdt[!cell_UMI<cell_UMI_thr]
            } else if (length(cell_UMI_thr) == 2) {
                cdt <- cdt[!cell_UMI<cell_UMI_thr[1] & !cell_UMI>cell_UMI_thr[2]]
            }
        }
        if (!is.null(cell_genes_thr)) {
            if (length(cell_genes_thr) == 1) {
                cdt <- cdt[!cell_genes_thr<cell_genes_thr]
            } else if (length(cell_genes_thr) == 2) {
                cdt <- cdt[!cell_genes_thr<cell_genes_thr[1] & !cell_genes_thr>cell_genes_thr[2]]
            }
        }

        # median values
        med <- cdt[,lapply(.SD,median),.(seq_batch_id,dataset),.SDcols=colnames(cdt)[-c(1:3)]]

		# save
		suff <- mat_id
		fwrite(cdt, file.path(scdb_dir, sprintf("cells_umi_counts_%s.tsv.gz", suff)), sep = "\t", compress = "gzip", quote = FALSE)
		fwrite(gdt, file.path(scdb_dir, sprintf("genes_umi_counts_%s.tsv.gz", suff)), sep = "\t", compress = "gzip", quote = FALSE)
		fwrite(med, file.path(scdb_dir, sprintf("median_umi_counts_%s.tsv", suff)), sep = "\t", quote = FALSE)

        # update matrix
        if (overwrite) {
            meta <- as.data.table(mat@cell_metadata, keep.rownames = "cell")
            for (cn in c("cell_UMI.x", "cell_genes.x", "cell_UMIs.y", "cell_genes.y")) {
                if (cn %in% colnames(meta)) {
                    meta[, (cn) := NULL]
                }
            }
            meta <- merge(meta, cdt, by = c("cell", "seq_batch_id", "dataset"), all.x = FALSE, sort = FALSE)
            class(meta) <- "data.frame"
            rownames(meta) <- meta$cell
            meta$cell <- NULL
            mat@cell_metadata <- meta
            mat@mat <- mat@mat[, rownames(meta)]
            mat@cells <- rownames(meta)
            mat@ncells <- ncol(mat@mat)
            scdb_add_mat(mat_id, mat)
        }
        
        list(cdt, gdt, med)
			
    }, simplify = FALSE, USE.NAMES = TRUE)
	names(sum_list) <- mat_ids
	
    # combine
    sum_cells <- rbindlist(lapply(sum_list, function(x) x[[1]]), idcol = "mat_id")
    sum_cells[, mat_id:=factor(mat_id, levels = mat_ids)]
    sum_genes <- rbindlist(lapply(sum_list, function(x) x[[2]]), idcol = "mat_id")
    sum_genes[, mat_id:=factor(mat_id, levels = mat_ids)]
    sum_median <- rbindlist(lapply(sum_list, function(x) x[[3]]), idcol = "mat_id")
    sum_median[, mat_id:=factor(mat_id, levels = mat_ids)]

    return(
        list(
            sum_cells = sum_cells,
            sum_genes = sum_genes,
            sum_median = sum_median
        )
    )

}