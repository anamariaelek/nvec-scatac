require(XML)
require(plyr)
require(dplyr)
require(DT)
require(data.table)
require(stringr)
require(ggplot2)
require(ggiraph)
require(patchwork)
require(rtracklayer)

# (1) loading datasets and general tools ---------------------------------------

# read rds
.readRDS <- function(file) {
  tryCatch(
    readRDS(file),
    error = function(e)
      stop("file:", file, "\n", e)
  )
}

# import metacell matrix
mc_mat_import <- function(file) {
  mcmat <- read.table(file=file, header=TRUE, row.names=1, sep ="\t")
  colnames(mcmat) <- gsub("[^[:digit:].]", "", colnames(mcmat))
  return(mcmat)
}

# import gene annotation
fread_gene_annotation <- function(
  file,
  cols=1:4,
  colnames=c("gene_id","gene name","PFAM domain","orthogroup")
) {
  GENE_ANNT <- fread(file=file, header=TRUE, sep ="\t", fill=TRUE, select=cols)
  setnames(GENE_ANNT, colnames)
  GENE_ANNT[,search_id:=do.call(paste,.SD),.SDcols=colnames]
  return(GENE_ANNT)
}

# import cell annotation
fread_cell_annotation <- function(
  file, 
  select=1:3, 
  colnames=c("mc","cell_type","color")
) {
  CELL_ANNT <- tryCatch(
    fread(file, sep = "\t", fill = TRUE, select = select),
    error = function(e)
      stop("file:", file, "\n", e)
  )
  setnames(CELL_ANNT, colnames)
  
  return(CELL_ANNT)
}

# reduce vector of metacells 
red_mc_vector <- function(x,range_sep=":") {
  all_mcs <- sort(as.integer(x))
  all_ir <- IRanges::IRanges(start = all_mcs, end = all_mcs)
  red_ir <- IRanges::reduce(all_ir)
  starts <- IRanges::start(red_ir)
  ends <- IRanges::end(red_ir)
  
  endsout <- ends
  for(i in 1:length(starts)) {
    if (starts[i] == ends[i])
      endsout[[i]] <- ""
  }
  outb <- paste(starts,endsout,sep=range_sep)
  outv <- str_remove(outb,paste0(range_sep,"$"))
  paste(outv,collapse=",")
}

# summarize annotations
summarize_cell_annotation <- function(annt) {
  
  tanns <- tapply(annt$mc, annt$cell_type, red_mc_vector)
  
  dt <- data.table(
    'cell type' = unique(annt$cell_type),
    metacells = tanns[unique(annt$cell_type)],
    cols = annt$color[match(unique(annt$cell_type),annt$cell_type)]
  ) 
  dt[,colshex := col2hex(cols)]
  dt %>%
    mutate( metacells = cell_spec(
      metacells, "html", color = "white", align = "c", 
      background = factor(`cell type`, dt$`cell type`, dt$colshex)
    )) %>%
    select(c("cell type","metacells")) %>%
    knitr::kable(escape = FALSE, align = "c")
}
summarize_dataset <- function(cols) {
  
  dt <- data.table(
    'dataset' = names(cols),
    'color' = cols
  ) 
  dt[,colshex := col2hex(color)]
  dt %>%
    mutate( dataset = cell_spec(
      dataset, "html", color = "white", align = "c", 
      background = factor(dataset, dt$dataset, dt$colshex)
    )) %>%
    select("dataset") %>%
    knitr::kable(escape = FALSE, align = "c")
}
summarize_cell_type <- function(cols) {
  
  dt <- data.table(
    'cell_type' = names(cols),
    'color' = cols
  ) 
  dt[,colshex := col2hex(color)]
  dt %>%
    mutate( cell_type = cell_spec(
      cell_type, "html", color = "white", align = "c", 
      background = factor(cell_type, dt$cell_type, dt$colshex)
    )) %>%
    select("cell_type") %>%
    knitr::kable(escape = FALSE, align = "c")
}

# (2) plotting tools ----------------------------------------------------------------------------------

#' Plot accessibility heatmap
#' 
plot_accessibility_heatmap <- function(
  mat,
  gene_list = NULL,
  heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
  show_cluster_names = TRUE,
  show_gene_names = TRUE,
  gene_font_size = 5,
  cluster_font_size = 5,
  clust_anno_height = unit(60,"mm"),
  clust_anno_width = unit(80,"mm"),
  clust_anno_gap = unit(1,"mm"),
  height = 10,
  width = 5,
  min_accessibility = NULL,
  max_accessibility = NULL,
  use_raster = TRUE,
  verbose = TRUE
) {
  
  # subset genes
  if (!is.null(gene_list))
    mat <- mat[rownames(mat) %in% gene_list,]
  # order columns
  hc_col <- hclust(dist(cor(mat,method="pearson")),method="ward.D2")
  mat <- mat[,hc_col$order]
  # order rows
  # hc_row <- hclust(dist(cor(t(mat),method="pearson")),method="ward.D2")
  gene_ord <- order(apply(mat,1,function(x) which.max(x)))
  mat <- mat[rev(gene_ord),]
  message("Plotting ",nrow(mat)," genes")
  # accessibility limits
  if (is.null(min_accessibility)) 
    min_accessibility <- min(mat[!is.infinite(mat)], na.rm = TRUE)
  mat <- pmax(mat,min_accessibility)
  if (is.null(max_accessibility)) 
    max_accessibility <- max(mat[!is.infinite(mat)], na.rm = TRUE)
  mat <- pmin(mat,max_accessibility)
  # colors
  col_fun = colorRampPalette(colors = heatmap_colors)
  shades = col_fun(40)
  # cluster labels
  collabs <- colnames(mat)
  if (!show_cluster_names) collabs <- rep("",length(collabs))
  column_ha_top = ComplexHeatmap::HeatmapAnnotation(
    which = "column",
    annotation_height = clust_anno_height,
    cluster = anno_text(
      which = "column", collabs, location = 0, just = "left",
      gp = gpar(fontsize = cluster_font_size, rot=90)
    )
  )
  column_ha_bottom = ComplexHeatmap::HeatmapAnnotation(
    which = "column",
    annotation_height = clust_anno_height,
    cluster = anno_text(
      which = "column", collabs, 
      gp = gpar(fontsize = cluster_font_size, rot=90)
    )
  )
  top_column_ha = c(column_ha_top, gap=clust_anno_gap)
  bottom_column_ha = c(column_ha_bottom, gap=clust_anno_gap)
  # gene labels
  gene_labels <- rownames(mat)
  if (!show_gene_names) gene_labels <- rep("",length(gene_labels))
  row_ha_right = ComplexHeatmap::HeatmapAnnotation(
    which = "row",
    annotation_width = clust_anno_width,
    gene = anno_text(
      which = "row", gene_labels, location = 0, just = "left",
      gp = gpar(fontsize = gene_font_size)
    )
  )
  row_ha_left = ComplexHeatmap::HeatmapAnnotation(
    which = "row",
    annotation_width = clust_anno_width,
    gene = anno_text(
      which = "row", gene_labels, location = 1, just = "right",
      gp = gpar(fontsize = gene_font_size)
    )
  )
  
  # heatmap
  h1 <- ComplexHeatmap::Heatmap(
    mat, name = "accessibility", col = shades, use_raster = use_raster,
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_title = "",
    row_title = sprintf( "%i marker genes", nrow(mat) ),
    right_annotation = row_ha_right,
    left_annotation = row_ha_left,
    top_annotation = top_column_ha,
    bottom_annotation = bottom_column_ha,
    height = 16,
    width = 8,
    border = TRUE
  )
  
  ComplexHeatmap::draw(h1)
  
}

#' Plot scATAC UMAP
#' 
#' @param df data.frame with two columns with UMAP coordinates, rownames should be cell names
#' @param color vector of colors of same length as `nrow(df)`
#' 
plot_umap <- function(
  df, 
  color = NULL,
  groups = NULL,
  labelMeans = TRUE,  
  defaultColor = "lightGrey",
  highlightPoints = NULL,
  size = 1,
  xlim = NULL, 
  ylim = NULL, 
  extend = 0.05, 
  xlabel = "UMAP1", 
  ylabel = "UMAP2",
  title = "", 
  downsample = TRUE,
  randomize = FALSE, 
  seed = 1,
  alpha = 1, 
  pointBorder = FALSE,
  showLabels = TRUE,
  showLegend = FALSE,
  labelAxisSize = 8,
  labelLegendSize = 8,
  labelLegendTrim = 28,
  legendSize = 4,
  legendPosition = "bottom",
  legendDirection = "vertical",
  legendRows = 3,
  ratioYX = 1, 
  bgWidth = 1,
  labelSize = 3,
  rastr = FALSE, 
  dpi = 300,
  ...
  
) {
  set.seed(seed)
  
  if(downsample==TRUE) {
    ds_id <- sample(1:nrow(df), size = pmin(1e4,0.2*nrow(df)))
    df <- df[ds_id,]
    color <- color[ds_id]
  }
  
  x=df[,1]
  y=df[,2]
  
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }
  
  if(randomize){
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  }else{
    idx <- seq_along(x)
  }
  
  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(color)){
      color <- color[include]
    }
  }
  
  if(is.null(xlim)){
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  
  if(is.null(ylim)){
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  

  if (!is.null(color)) {
    
    stopifnot(length(color) == nrow(df))
    
    # highlights
    if(!is.null(highlightPoints)){
      if(downsample==TRUE)
        highlightPoints <- intersect(highlightPoints,ds_id)
      if(length(highlightPoints) < length(color)){
        color[-highlightPoints] <- NA
        idx <- c(idx[-highlightPoints], idx[highlightPoints])
      }
    }
    
    df$color <- color
    df$group <- names(color)
    
    if (is.null(groups)) {
      group_values <- sort(unique(names(color)))
    } else {
      group_values <- groups
    }
    pallete <- color[match(group_values, names(color))]
    names(pallete) <- group_values
    df$group <- factor(df$group, levels = group_values)
        
  } else {
    
    df$color <- defaultColor
    df$group <- ""
    
  }
  
  
  dp <- ggplot(df[idx,], aes(x = x, y = y)) + 
    coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
    xlab(xlabel) + ylab(ylabel) + 
    ggtitle(title)
    # .geom_point_rast2(
    #   size = size, raster.dpi = dpi, alpha = alpha, 
    #   raster.width = min(par('fin')),
    #   raster.height = (ratioYX * min(par('fin')))
    # ) 
  if (pointBorder) {
    dp <- dp + geom_point_interactive(
        aes(fill = group, tooltip = group, data_id = group), 
        size = size, alpha = alpha, pch = 21
      ) + scale_fill_manual(values = pallete)
  } else {
    dp <- dp + geom_point_interactive(
      aes(color = group, tooltip = group, data_id = group), 
      size = size, alpha = alpha
    ) + scale_color_manual(values = pallete) 
  }
  dp <- dp +
    theme_classic() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    axis.title = element_text(size = labelAxisSize)
  ) 
  
  if (showLegend==FALSE) {
    dp <- dp + theme(legend.position = "none")
  } else {
    dp <- dp + theme(
      legend.position = legendPosition,
      legend.title = element_blank(), 
      legend.text = element_text(size = labelLegendSize),
      legend.direction = legendDirection,
    ) +
      scale_color_manual(
        values = pallete, 
        labels = function(x) ifelse(nchar(x)>labelLegendTrim, paste0(substr(x,1,labelLegendTrim),"..."), x)
      ) +
      scale_fill_manual(
        values = pallete, 
        labels = function(x) ifelse(nchar(x)>labelLegendTrim, paste0(substr(x,1,labelLegendTrim),"..."), x)
      ) +
      guides(
        colour = guide_legend(nrow=legendRows, override.aes = list(size = legendSize, alpha = 1)),
        fill = guide_legend(nrow=legendRows, override.aes = list(size = legendSize, alpha = 1))
      )
  }
  if (showLabels==TRUE) {
    dt <- as.data.table(df)
    dtlab <- dt[,lapply(.SD,median),by=group,.SDcols=c("x","y")]
    dp <- dp + geom_text_repel(
      data=dtlab, 
      aes_string("x", "y", label="group"), 
      box.padding = 0.5, 
      max.overlaps = Inf,
      bg.color = "white", 
      bg.r = 0.15 
    )
  }
  dp

  
}

#' Plot gene accessibility on UMAP
#' 
#' @param df data.frame with two columns with UMAP coordinates, rownames should be cell names
#' 
plot_umap_gene <- function(
  df, 
  sc_values,
  color_scale = c("gray95","lightyellow","khaki1","orange","orangered2","#520c52"),
  group_values = NULL,
  labelMeans = TRUE,  
  size = 1,
  alpha = 1, 
  printBorder = FALSE,
  xlim = NULL, 
  ylim = NULL, 
  extend = 0.05, 
  xlabel = "UMAP1", 
  ylabel = "UMAP2",
  title = "", 
  downsample = TRUE,
  legend.position = "bottom", # c(0.9, 0.1)
  randomize = FALSE, 
  seed = 1,
  sc_min = NULL,
  sc_max = NULL,
  labelAxisSize = 8,
  labelLegendSize = 8,
  ratioYX = 1, 
  rastr = FALSE, 
  interactive = TRUE,
  dpi = 300,
  ...
  
) {
  
  set.seed(seed)
  
  x=df[,1]
  y=df[,2]
  
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }

  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(sc_values))
      sc_values <- sc_values[include]
    if (!is.null(group_values)) 
      group_values <- group_values[include]
      
  }
  
  # add groups to df
  if (!is.null(group_values)) {
    df$group <- group_values
  } else {
    df$group <- ""
  }
  
  # limits for sc values
  if(is.null(xlim))
    xlim <- range(df$x) %>% extendrange(f = extend)
  if(is.null(ylim))
    ylim <- range(df$y) %>% extendrange(f = extend)
  
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  
  # sc values
  if (is.null(sc_min))
    sc_min = min(sc_values[!is.infinite(sc_values)&!is.na(sc_values)], na.rm = TRUE)
  message(sprintf("min: %s",sc_min))
  if (is.null(sc_max))
    sc_max = max(sc_values[!is.infinite(sc_values)&!is.na(sc_values)], na.rm = TRUE)
  message(sprintf("max: %s",sc_max))
    
  
  # apply min/max to sc vector
  sc_values [ sc_values < sc_min ] = sc_min
  sc_values [ sc_values > sc_max ] = sc_max
  
  # add sc values to df
  df$accessibility <- sc_values
  
  # order for plotting
  if (randomize){
    set.seed(seed)
    idx <- sample(seq_along(sc_values), length(sc_values))
  } else {
    idx <- order(sc_values, decreasing=FALSE)
  }
  
  # optionally downsample
  if(downsample==TRUE) {
    idx <- rev(rev(idx)[1:pmin(1e4,0.2*nrow(df))])
    # df <- df[ds_id,]
    # sc_values <- sc_values[ds_id]
    # group_values <- group_values[ds_id]
  }
  
  dp <- ggplot(df[idx,], aes(x = x, y = y))
  if (printBorder) {
    dp <- dp + geom_point_interactive(
      aes(fill = accessibility, tooltip = group, data_id = group), 
      size = size, alpha = alpha, pch = 21, color = "black"
    ) + scale_fill_gradientn(colours = color_scale) 
  } else {
    dp <- dp + geom_point_interactive(
      aes(color = accessibility, tooltip = group, data_id = group), 
      size = size, alpha = alpha
    ) + scale_color_gradientn(colours = color_scale) 
  }
  dp <- dp +
    coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) + 
    labs(x = xlabel, y = ylabel) + 
    ggtitle(title) +
    theme_classic() + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.title = element_text(size = labelAxisSize),
      legend.title = element_text(size=labelLegendSize, vjust=0.5,hjust=0.5),
      legend.text = element_text(size=labelLegendSize),
      legend.direction = "horizontal",
      legend.position = legend.position
    )
  dp
}

# Accessibility barplot
plot_gene_access_bar <- function(
    df, 
    x_order = NULL,
    y_log = FALSE,
    group_column = NULL, 
    group_order = NULL,
    col_pal = NULL,
    barPosition = "stack",
    showLabels = TRUE,
    showCounts = FALSE,
    legendPosition = "none",
    legendRows = 3,
    legendSize = 4
) {
  
  dt <- setDT(copy(df))
  if (!is.null(x_order))
    dt[[1]] <- factor(dt[[1]], levels = x_order)
  
  if (!is.null(group_column)) {
    if (!is.null(group_order)) {
      dt[[group_column]] <- factor(dt[[group_column]], levels=group_order)
      counts_dt <- dt[,.(N=.N),c(colnames(dt)[1],group_column)]
      
    } else {
      dt[[group_column]] <- factor(dt[[group_column]], levels=unique(dt[[group_column]]))
      counts_dt <- dt[,.(N=.N),c(colnames(dt)[1])]
    }
    dt <- merge.data.table(dt, counts_dt, by = intersect(colnames(dt), colnames(counts_dt)))
    
    if (is.null(x_order)) {
      setorderv(dt, group_column)
      dt[[1]] <- factor(dt[[1]], levels = unique(dt[[1]]))
    }
    as <- aes_string(x = colnames(dt)[1], fill = group_column)
  } else {
    as <- aes_string(x = colnames(dt)[1])
  }
  
  gb <- ggplot(dt, as) + 
    geom_bar_interactive(aes(tooltip = get(group_column), data_id = get(group_column)), position = barPosition)
  # if (y_log) {
  #   gb <- gb + scale_y_log10(expand = expansion(mult = c(0,0.1)))
  # } else {
    gb <- gb + scale_y_continuous(expand = expansion(mult = c(0,0.4)))
  # }
  if (showCounts) {
    counts_dt <- unique(dt)[,.(N=sum(N)),c(colnames(dt)[1])]
    as_lab <- aes_string(x = colnames(counts_dt)[1], y = "N", label = "N")
    gb <- gb + geom_text(data = counts_dt, as_lab, inherit.aes = FALSE, angle = 45, vjust = 0, hjust = 0.1)
  }
  gb <- gb +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "lightgray"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = legendPosition
    ) +
    guides(
      fill = guide_legend(nrow=legendRows, override.aes = list(size = legendSize, alpha = 1))
    ) +
    labs(y = switch(barPosition, "stack" = "count", "fill" = "fraction", "count"))
  if (showLabels) {
    gb <- gb + theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    )
  } else {
    gb <- gb + theme(
      axis.text.x = element_blank()
    ) 
  }
  if (!is.null(col_pal))
    gb <- gb + 
      scale_fill_manual(values = col_pal)
  
}

# Accessibility boxplot
plot_gene_access_box <- function(
    df, 
    group_column = NULL, 
    col_pal = NULL,
    showLabels = TRUE,
    legendPosition = "none"
) {
  
  dt <- setDT(copy(df))
  
  if (!is.null(group_column)) {
    dt[,group:=factor(dt[[group_column]], levels=unique(dt[[group_column]]))]
    setorderv(dt, group_column)
    dt[[1]] <- factor(dt[[1]], levels = unique(dt[[1]]))
    as <- aes_string(x = colnames(dt)[1], y = colnames(dt)[2], fill = group_column)
  } else {
    as <- aes_string(x = colnames(dt)[1], y = colnames(dt)[2])
  }
  gb <- ggplot(dt, as) + 
    geom_violin_interactive(aes(tooltip = as.character(group), data_id = group), scale = "width", alpha = 0.8, color = "black") +
    geom_boxplot_interactive(aes(tooltip = as.character(group), data_id = group), width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black") + 
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    theme_classic() + 
    theme(
      panel.grid.major.y = element_line(size = 0.5, color = "lightgray"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = legendPosition
    )
  if (showLabels) {
    gb <- gb + theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    )
  } else {
    gb <- gb + theme(
      axis.text.x = element_blank()
    )
  }
  if (!is.null(col_pal))
    gb <- gb + 
    scale_fill_manual(values = col_pal)
  
}

# Plot TF activity vs expression
#' @param chr_exp_ct data.table with TF activity and expression
#' @param gn gene id
#' @param ct_cols cell type colors
plot_tf_act_exp <- function(chr_exp_ct, gn, ct_cols) {

  #chr_exp_ct <- fread("Nematostella_scATAC/Results/GRN/metacell/gene_expression_fc_chromVAR_genes_exp_FC2_acc_FC4_spearman.tsv.gz")
  #chr_exp_ct <- chr_exp_ct[insilico_ChIP_threshold == 0.1]
  #setnames(chr_exp_ct, c("motif", "zscore"), c("archetype_name", "motif_deviation"))
  #chr_exp_ct <- chr_exp_ct[, .(
  #  expression = mean(expression),
  #  motif_deviation = mean(motif_deviation)
  #), .(gene, gene_name, common_name, pfam, og, archetype_name, cell_type, stage)]
  # aggregate per cell type
  #chr_exp_ct <- chr_exp_ct[, .(
  #  expression = mean(expression),
  #  motif_deviation = mean(motif_deviation)
  #), .(gene, gene_name, common_name, pfam, og, archetype_name, cell_type, stage)]
  #saveRDS(chr_exp_ct, "Data/TF_activity_expression.rds")

  tf_dt <- chr_exp_ct[gene==gn]
  if (tf_dt$common_name[1] != "") gn <- sprintf("%s (%s)", gn, tf_dt$common_name[1])
  pf <- tf_dt$pfam[1]
  og <- tf_dt$og[1]
  tf_tt <- sprintf("%s\nPFAM: %s\nOG: %s", gn, pf, og)
  tf_gp <- ggplot(tf_dt, aes(
          expression, motif_deviation, 
          label = cell_type,
          fill = cell_type, 
          color = cell_type,
          shape = stage
      )
    ) +
    geom_point(size = 4) +
    #geom_point_interactive(
    #  aes(tooltip = cell_type, data_id = cell_type), 
    #  size = 4, alpha = alpha, pch = 21
    #) +
    ggrepel::geom_text_repel(size = 4, alpha = 0.6) +
    scale_fill_manual(values = ct_cols) +
    scale_color_manual(
      values = c(
        structure(
          colorspace::darken(ct_cols, 0.5),
          names = names(ct_cols)
        ),
        structure(
          colorspace::lighten(ct_cols, 0.5),
          names = paste0(names(ct_cols), "_sc")
        )
      )
    ) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_shape_manual(
        values = c("adult" = 21, "gastrula" = 24)
    ) +
    labs(
        x = "TF expression",
        y = "TF activity",
        title = tf_tt
    ) +
    theme_classic() + 
    theme(
      # strip.text = element_text(size = 8),
      # axis.text = element_text(size = 8),
      # axis.title = element_text(size = 8),
      # title = element_text(size = 8),
      legend.position = "none"
    )

    tf_gp
}

#' Select genes for heatmap
genes_select_dt <- function(sterm, nmat, annt) {
  sterms <- strsplit(sterm, ',')[[1]]
  gsfid <- unique(unlist(lapply(
    sterms, function(sterm) {
      sterm = stringr::str_trim(sterm,side="both")
      grep1 <- grep(sterm,annt[[1]],ignore.case=TRUE)
      grep2 <- grep(sterm,annt[[2]],ignore.case=TRUE)
      grep3 <- grep(sterm,annt[[3]],ignore.case=TRUE)
      rids <- sort(unique(c(grep1, grep2, grep3)))
      gs <- annt[rids][[1]]
      gsf <- gs[gs %in% rownames(nmat)]
      match(gsf,annt[[1]])
    }), use.names = FALSE))
  annt[gsfid,1:3]
}
genes_upload_dt <- function(gs,annt) {
  gsfid <- match(gs,annt[[1]])
  annt[gsfid,1:3]
}
genes_select_names <- function(dt,rid) {
  dt[rid,][[1]]
}

#' Plot multigene heatmap to show expression of selected genes
mgenes_hmap_height <- function(nmat, gids, annt, min_expression_fc = NULL){
  tryCatch({
    gs <- intersect(gids,rownames(nmat))
    if (!is.null(min_expression_fc)) {
      flt <- apply(nmat[gs,], 1, function(x) !(sort(x,decreasing=TRUE,na.last=TRUE)[1]<min_expression_fc))
      gs <- gs[flt]
    }
    message("Length: ",length(gs))
    return(length(gs))
  }, error=function(e) {
    message(e)
    message("Using fixed length")
    return(2)
  })
}
mgenes_hmap <- function(
    nmat, annt, gids,
    min_expression_fc = NULL,  # gene filtering
    scale_expression_fc = 4, # trim expression values, i.e. set anything > scale_expression_fc to this value
    heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
    ct_table, cell_type_palette = NULL, cluster_genes = TRUE,
    gene_font_size = 12,mcid_font_size = 12, mc_annotaion_height = unit(2, "mm")
){
  
  # selected genes
  message("Selected genes: ",length(gids))
  gs <- intersect(gids,rownames(nmat))
  rid <- match(gs, annt[[1]])
  gns <- annt[rid][[2]]
  badgns <- gns=="" | is.na(gns)
  gns[badgns] <- gs[badgns]
  row_labels <- gns
  names(row_labels) <- gs
  
  message("Genes in heatmap: ",nrow(nmat[gs,]))
  hm <- as.matrix(nmat[gs,,drop=FALSE])
  
  # filter genes
  if (!is.null(min_expression_fc)) {
    message("Filtering genes by min FC ", min_expression_fc)
    flt <- apply(hm, 1, function(x) !(sort(x,decreasing=TRUE,na.last=TRUE)[1]<min_expression_fc))
    hm <- hm[flt,]
  }
  
  # expression matrix
  hm <- pmin(hm,scale_expression_fc)
  hm[is.na(hm)] <- 0
  hm[is.nan(hm)] <- 0
  hm[is.infinite(hm)] <- 0
  
  # order genes
  if (cluster_genes == TRUE) {
    message("Clustering genes")
    gord <- order(apply(hm, 1,function(x) which.max(rollmean(x,1))))
    hm <- hm[gord,]
  }
  
  # matrix
  hmat <- cbind(hm)
  
  # heatmap colors
  if (is.null(min_expression_fc)) min_expression_fc=0; max_expression_fc=5
  max_expression_fc <- pmin(scale_expression_fc, max_expression_fc)
  message("Scaling colors for heatmap between ", min_expression_fc , " and ", max_expression_fc)
  col_fun = circlize::colorRamp2(
    breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)),
    colors = heatmap_colors
  )
  
  # cell type colours
  message("Cell types colors")
  if (is.null(cell_type_palette)) {
    ctpalette_dt <- unique(ct_table[,.(cell_type,color)])
    cell_type_palette <- structure(ctpalette_dt$color, names = ctpalette_dt$cell_type)
  }
  # all cell types and colors
  cell_types <- ct_table[match(colnames(hmat),mc)]$cell_type
  cell_colours <- cell_type_palette[cell_types]
  
  # cell type annotiation bar
  message("Cell types annotations")
  ncts_ann <- length(cell_types)
  # print(sprintf("length col ann %s", ncts_ann))
  ncts_hm <- ncol(hmat)
  # print(sprintf("ncol hmat %s", ncts_hm))
  if (ncts_hm == ncts_ann) {
    ct_ann <- ComplexHeatmap::columnAnnotation(
      ct = cell_types, col = list(ct=cell_colours), height = mc_annotaion_height,
      gp = gpar(fontsize = mcid_font_size), border = FALSE,
      show_annotation_name = FALSE, show_legend = FALSE
    )
  } else {
    ct_ann <- ComplexHeatmap::columnAnnotation(
      anno_empty(which = "column")
    )
  }
  
  # gene annotation
  message("Gene annotations")
  gs <- rownames(hm)
  annid <- match(gs,annt[[1]])
  ganns <- annt[annid][[3]]
  ganns[is.na(ganns)] <- ""
  ganns_tr <- unlist(lapply(ganns, function(x){
    if(nchar(x)>30) {
      paste0(substr(x,1,27),"...")
    } else {
      x
    }
  }))
  names(ganns_tr) <- gs
  if (any(ganns_tr!="")) {
    gs_ann <- ComplexHeatmap::HeatmapAnnotation(
      which = "row", gn = anno_text(ganns_tr,which="row", gp = gpar(fontsize = gene_font_size)),
      show_annotation_name = FALSE, show_legend = FALSE,
      gp = gpar(fontsize = gene_font_size)
    )
  } else {
    gs_ann <- ComplexHeatmap::HeatmapAnnotation(
      which = "row", gn = anno_empty(which="row", border = FALSE),
      show_annotation_name = FALSE, show_legend = FALSE,
    )
  }
  
  # expression heatmap
  message("Building heatmap")
  ht1 <- ComplexHeatmap::Heatmap(
    hmat, name = "expression FC", show_heatmap_legend = TRUE,
    cluster_columns = FALSE, cluster_rows = FALSE,
    #row_labels = row_labels[rownames(hm)],
    show_column_dend = FALSE, show_row_dend = FALSE,
    show_column_names = TRUE, show_row_names = TRUE,
    row_names_side = "left", column_names_side = "bottom",
    column_names_gp = gpar(fontsize = mcid_font_size),
    row_names_gp = gpar(fontsize = gene_font_size),
    border = TRUE,
    col = col_fun, rect_gp = gpar(col = "gray88", lwd = 0.1),
    bottom_annotation = ct_ann, top_annotation = ct_ann,
    right_annotation = gs_ann,
    heatmap_legend_param = list(
      legend_direction = "horizontal",
      legend_width = unit(5, "cm"),
      border = TRUE
    )
  )
  draw(ht1, heatmap_legend_side = "bottom")
}

# Adapted from 
# https://github.com/tidyverse/ggplot2/blob/660aad2db2b3495ae0d8040915a40d247133ffc0/R/geom-point.r
# from https://github.com/VPetukhov/ggrastr/blob/master/R/geom-point-rast.R
# This funciton now handles issues with Cairo installation that can lead to plot errors
.geom_point_rast2 <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  raster.width = min(par('fin')), 
  raster.height = min(par('fin')), 
  raster.dpi = 300
){
  
  GeomPointRast <- tryCatch({
    
    if(!.checkCairo()){
      stop()
    }
    
    #Try to create a geom rast for points if not then just use normal geom_point
    ggplot2::ggproto(
      "GeomPointRast",
      ggplot2::GeomPoint,
      required_aes = c("x", "y"),
      non_missing_aes = c("size", "shape", "colour"),
      default_aes = aes(
        shape = 19, colour = "black", size = 1.5, fill = NA,
        alpha = NA, stroke = 0.5
      ),
      
      draw_panel = function(data, panel_params, coord, na.rm = FALSE, 
                            raster.width=min(par('fin')), raster.height=min(par('fin')), raster.dpi=300){
        
        #From ggrastr  
        prevDevID <- dev.cur()
        
        p <- ggplot2::GeomPoint$draw_panel(data, panel_params, coord)
        
        devID <- Cairo::Cairo(
          type='raster', 
          width=raster.width*raster.dpi, 
          height=raster.height*raster.dpi, 
          dpi=raster.dpi, 
          units='px', 
          bg="transparent"
        )[1]
        
        grid::pushViewport(grid::viewport(width=1, height=1))
        
        grid::grid.points(
          x=p$x, 
          y=p$y, 
          pch = p$pch, 
          size = p$size,
          name = p$name, 
          gp = p$gp, 
          vp = p$vp, 
          draw = TRUE
        )
        
        grid::popViewport()
        gridCapture <- grid::grid.cap()
        
        dev.off(devID)
        
        dev.set(prevDevID)
        
        grid::rasterGrob(
          gridCapture, 
          x=0, 
          y=0, 
          width = 1,
          height = 1,
          default.units = "native",
          just = c("left","bottom")
        )
        
      }
      
    )
    
  }, error = function(e){
    
    if(.checkCairo()){
      message("WARNING: Error found with trying to rasterize geom. Continuing without rasterization.")
    }else{
      message("WARNING: Error found with Cairo installation. Continuing without rasterization.")
    }
    
    #Default geom_point
    ggplot2::ggproto(
      "GeomPoint", 
      ggplot2::GeomPoint,
      required_aes = c("x", "y"),
      non_missing_aes = c("size", "shape", "colour"),
      default_aes = aes(
        shape = 19, colour = "black", size = 1.5, fill = NA,
        alpha = NA, stroke = 0.5
      ),
      
      draw_panel = function(data, panel_params, coord, na.rm = FALSE, 
                            raster.width=min(par('fin')), raster.height=min(par('fin')), raster.dpi=300){
        if (is.character(data$shape)) {
          data$shape <- ggplot2:::translate_shape_string(data$shape) #Hidden ggplot2
        }
        
        coords <- coord$transform(data, panel_params)
        
        pGrob <- grid::pointsGrob(
          x = coords$x, 
          y = coords$y,
          pch = coords$shape,
          gp = grid::gpar(
            col = scales::alpha(coords$colour, coords$alpha),
            fill = scales::alpha(coords$fill, coords$alpha),
            # Stroke is added around the outside of the point
            fontsize = coords$size * .pt + coords$stroke * .stroke / 2,
            lwd = coords$stroke * .stroke / 2
          )
        )
        
        pGrob
        
      },
      
      draw_key = ggplot2::draw_key_point
    )
    
    
  })
  
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomPointRast,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      raster.width=raster.width,
      raster.height=raster.height,
      raster.dpi=raster.dpi,
      ...
    )
  )
  
}

.checkCairo <- function(){
  tryCatch({
    tmp <- dev.cur()
    Cairo::Cairo(type='raster')
    dev.off()
    dev.set(tmp)
    TRUE
  }, error = function(e){
    FALSE
  })
}

# (3) calculations ----------------------------------------------------------------------------------

# generate table for highly expressed genes in a group of metacells (mcs)
# method="absolute" - all mcs (- lky percentage) must have gene fc above threshold
# method="median" - median in all mcs must be above threshold
mc_gene_summary <- function(
  mc_ids, fc=2, method=c("absolute","median"), methodbg=method, lky=0, 
  usefcbg=FALSE, fcbg=fc, lkybg=lky, 
  mc_fp, mc_counts, annt, tfannt
){
  mc_counts <- as.data.frame(mc_counts)
  mcs <- colnames(mc_fp)
  mc_ids <- as.character(mc_ids)
  mc_others <- mcs[!(mcs %in% mc_ids)]
  if (length(method)>1) method=method[1]
  if (!method %in% c("absolute","median")) stop("Method must be either 'absolute' or 'median'!")
  
  # horizontal UMI fraction
  #tot_umis <- rowSums(mc_counts)
  if (length(mc_ids)==1) {
    umisum_selected <- mc_counts[[mc_ids]]
    names(umisum_selected) <- rownames(mc_counts)
  } else {
    umisum_selected <- rowSums(cbind(mc_counts[,mc_ids]))
  }
  umisum_others <- rowSums(mc_counts[mc_others])
  umi_frac <- umisum_selected/(umisum_selected+umisum_others)
  umi_frac[is.nan(umi_frac)] <- 0
  
  # median gene FC
  median_fc <- apply(cbind(mc_fp[,mc_ids]),1,median,na.rm=TRUE)
  names(median_fc) <- rownames(mc_fp)
  
  # "gap genes" - specifically expressed in selected mcs
  if (method == "absolute") {
    
    ntarget <- length(mc_ids) - length(mc_ids)*lky
    abovefc <- rowSums(cbind(mc_fp[,mc_ids]) > fc) 
    gap_genes_abovefc <- abovefc > ntarget | abovefc == ntarget
    if (usefcbg==TRUE) {
      if (methodbg == "absolute") {
        ntargetbg <- length(mc_others) - lkybg*length(mc_others)
        belowfc <- rowSums(cbind(mc_fp[,mc_others]) < fcbg)
        gap_genes_belowfc <- belowfc > ntargetbg | belowfc== ntargetbg
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      } else if ( methodbg == "median") {
        gap_genes_belowfc <-  apply(cbind(mc_fp[,mc_others]),1,function(x) median(x,na.rm=TRUE)<fcbg)
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      }
    } else {
      fcbg <- NULL
      gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc==TRUE)]
    }
    gap_genes_annot <- intersect(gap_genes,annt$gene_id)
    anntids <- match(gap_genes_annot,annt$gene_id)
    
  } else if (method == "median") {
    
    gap_genes_abovefc <- apply(cbind(mc_fp[,mc_ids]),1,function(x) median(x,na.rm=TRUE)>fc)
    if (usefcbg==TRUE) {
      if (methodbg == "absolute") {
        ntargetbg <- length(mc_others) - lkybg*length(mc_others)
        belowfc <- rowSums(cbind(mc_fp[,mc_others]) < fcbg)
        gap_genes_belowfc <- belowfc > ntargetbg | belowfc== ntargetbg
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      } else if ( methodbg == "median") {
        gap_genes_belowfc <-  apply(cbind(mc_fp[,mc_others]),1,function(x) median(x,na.rm=TRUE)<fcbg)
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      }
    } else {
      fcbg <- NULL
      gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc==TRUE)]
    }
    
    
    gap_genes_annot <- intersect(gap_genes,annt$gene_id)
    anntids <- match(gap_genes_annot,annt$gene_id)
    
  }
  
  # gene table
  gap_dt <- annt[anntids,1:3]
  gap_dt[,':='(
    tot_umi = umisum_selected[gap_genes_annot],
    umi_frac = round(umi_frac[gap_genes_annot]*100,2),
    median_fc = round(median_fc[gap_genes_annot],2),
    tf = ifelse(gene_id %in% tfannt$gene,"yes","no")
  )]
  gap_dt[,tf:=factor(tf,levels=c("yes","no"))]
  setorder(gap_dt,-median_fc)
  setnames(gap_dt,c("gene_id","tot_umi","umi_frac","median_fc"),c("gene ID","total UMIs","% UMIs","median fc"))
  
  # summary of the search and results
  search_dt <- data.table(
    `selected` = paste( as.character(mc_ids), collapse=", "),
    `minimum fold change in selected metacells` = fc, 
    `maximum fold change in non-selected metacells` = fcbg,
    `number of genes` = length(gap_genes_annot)
  )
  search_tdt <- data.table::transpose(search_dt, keep.names = "V0")
  setnames(search_tdt,c(""," "))
  # output
  list(gene_summary=gap_dt, summary=search_tdt)
}

# same for clusters in scatac
atac_gene_summary <- function(
  mc_ids, fc=2, method=c("absolute","median"), methodbg=method, lky=0, 
  usefcbg=FALSE, fcbg=fc, lkybg=lky, 
  mc_fp, mc_counts, annt, tfannt
){
  annt <- setDT(annt[,1:3])
  setnames(annt, c("gene_id","PFAM_domain","gene_name"))
  
  tfannt <- setDT(tfannt[,1:3])
  setnames(tfannt, c("gene_id","PFAM_domain","gene_name"))
  
  mc_counts <- as.data.frame(mc_counts)
  mcs <- colnames(mc_fp)
  mc_ids <- as.character(mc_ids)
  mc_others <- mcs[!(mcs %in% mc_ids)]
  if (length(method)>1) method=method[1]
  if (!method %in% c("absolute","median")) stop("Method must be either 'absolute' or 'median'!")
  
  # horizontal UMI fraction
  if (length(mc_ids)==1) {
    umisum_selected <- mc_counts[,mc_ids]
    names(umisum_selected) <- rownames(mc_counts)
  } else {
    umisum_selected <- rowSums(cbind(mc_counts[,mc_ids]))
  }
  umisum_others <- rowSums(mc_counts[mc_others])
  umi_frac <- umisum_selected/(umisum_selected+umisum_others)
  umi_frac[is.nan(umi_frac)] <- 0
  
  # median gene FC
  median_fc <- apply(cbind(mc_fp[,mc_ids]),1,median,na.rm=TRUE)
  names(median_fc) <- rownames(mc_fp)
  
  # "gap genes" - specific for selected groups
  if (method == "absolute") {
    
    ntarget <- length(mc_ids) - length(mc_ids)*lky
    abovefc <- rowSums(cbind(mc_fp[,mc_ids]) > fc) 
    gap_genes_abovefc <- abovefc > ntarget | abovefc == ntarget
    if (usefcbg==TRUE) {
      if (!methodbg %in% c("absolute","median")) stop("Method must be either 'absolute' or 'median'!")
      message(sprintf("Getting %s accessible genes with %s background",method,methodbg))
      if (methodbg == "absolute") {
        ntargetbg <- length(mc_others) - lkybg*length(mc_others)
        belowfc <- rowSums(cbind(mc_fp[,mc_others]) < fcbg)
        gap_genes_belowfc <- belowfc > ntargetbg | belowfc== ntargetbg
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      } else if ( methodbg == "median") {
        gap_genes_belowfc <-  apply(cbind(mc_fp[,mc_others]),1,function(x) median(x,na.rm=TRUE)<fcbg)
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      }
    } else {
      message(sprintf("Getting %s accessible genes",method))
      fcbg <- NULL
      gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc==TRUE)]
    }
    gap_genes_annot <- intersect(gap_genes,annt$gene_id)
    anntids <- match(gap_genes_annot,annt$gene_id)
    
  } else if (method == "median") {
    
    gap_genes_abovefc <- apply(cbind(mc_fp[,mc_ids]),1,function(x) median(x,na.rm=TRUE)>fc)
    if (usefcbg==TRUE) {
      if (methodbg == "absolute") {
        ntargetbg <- length(mc_others) - lkybg*length(mc_others)
        belowfc <- rowSums(cbind(mc_fp[,mc_others]) < fcbg)
        gap_genes_belowfc <- belowfc > ntargetbg | belowfc== ntargetbg
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      } else if ( methodbg == "median") {
        gap_genes_belowfc <-  apply(cbind(mc_fp[,mc_others]),1,function(x) median(x,na.rm=TRUE)<fcbg)
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      }
    } else {
      fcbg <- NULL
      gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc==TRUE)]
    }
    gap_genes_annot <- intersect(gap_genes,annt[[1]])
    anntids <- match(gap_genes_annot,annt[[1]])
    
  }
  
  # gene table
  gap_dt <- annt[anntids,1:3]
  gap_dt[,':='(
    access_sum = round(umisum_selected[gap_genes_annot],2),
    access_frac = round(umi_frac[gap_genes_annot]*100,2),
    median_fc = round(median_fc[gap_genes_annot],2),
    tf = ifelse(gene_id %in% tfannt[[1]],"yes","no")
  )]
  gap_dt[,tf:=factor(tf,levels=c("yes","no"))]
  setorder(gap_dt,-median_fc)
  gap_dt[nchar(gene_name)>60, gene_name:=sprintf("%s...",substr(gene_name,1,57))]
  gap_dt[nchar(PFAM_domain)>60, PFAM_domain:=sprintf("%s...",substr(PFAM_domain,1,57))]
  setnames(
    gap_dt,
    c("gene_id","PFAM_domain","gene_name","access_sum","access_frac","median_fc"),
    c("gene ID","PFAM domain","gene name","acessibility sum","acessibility fraction","median FC")
  )
  
  # summary of the search and results
  search_dt <- data.table(
    `selected metacells` = paste( as.character(mc_ids), collapse=", "),
    `minimum fold change in selected metacells` = fc, 
    `maximum fold change in non-selected metacells` = fcbg,
    `number of genes` = length(gap_genes_annot)
  )
  search_tdt <- data.table::transpose(search_dt, keep.names = "V0")
  setnames(search_tdt,c(""," "))
  # output
  list(gene_summary=gap_dt, summary=search_tdt)
}

# generate table for highly expressed genes in selected cell type(s)
cell_type_gene_summary <- function() {
  # TBA
}

# imputation
impute_mat <- function(mat, weightList) {
  
  imputeMat <- lapply(seq_along(weightList), function(x){
    h5df <- h5ls(weightList[[x]])
    blocks <- gtools::mixedsort(grep("block",h5df$name,value=TRUE))
    matx <- lapply(seq_along(blocks), function(y){
      #Read In Weights and Names
      bn <- h5read(weightList[[x]], paste0(blocks[y], "/Names"))
      by <- h5read(weightList[[x]], paste0(blocks[y], "/Weights"))
      colnames(by) <- bn
      rownames(by) <- bn
      #Multiply
      Matrix::t(by %*% Matrix::t(mat[, paste0(bn), drop = FALSE]))
    }) %>% Reduce("cbind", .)
    matx[, colnames(mat)]
  }) %>% Reduce("+", .)
  
  imputeMat / length(weightList)
}

#' Function to get a set of features in overlap between multiple cell types,
#' which are also not present in other cell types
#' 
#' @param ct_fg character cell types for which to find overlap of features 
#'   (i.e. foreground cell types)
#' @param ct_bg character cell types to use as background
#' @param feats_dt data.frame with features annotations, it should contain 
#'   at least the following columns: feat_column, val_column and cell_type
#' @param feat_column character, feature for which to get overlaps
#' @param val_column character column with values to use for selecting features
#' @param min_fg numeric minimum value in selected cell types
#' @param max_bg numeric maximum value in other cell types
#' @param level_dt data.frame with cell_type and broad_cell_type mapping
#'   e.g. data.table(cell_type = c("muscle_1", "muscle_2"), broad_cell_type = c("muscle", "muscle"))
#' @param feats_blacklist
#' 
#' @return data.table with features in overlap
#' 
feats_ct_ovl <- function(
    ct_fg, ct_bg = NULL,
    feats_dt, feat_column = "peak", val_fg_column = "Log2FC", val_bg_column = "Log2FC",
    min_fg = 1, max_bg = 1, level_dt = NULL, feats_blacklist = NULL
) {
  
  # we need these columns
  feats_wdt <- copy(feats_dt)
  setDT(feats_wdt)
  feats_wdt[, feature := feats_dt[[feat_column]]]
  feats_wdt[, value_fg := feats_dt[[val_fg_column]]]
  feats_wdt[, value_bg := feats_dt[[val_bg_column]]]
  feats_cols <- c("feature", "cell_type", "value_fg", "value_bg")
  feats_wdt <- unique(feats_wdt[, ..feats_cols])
  setorder(feats_wdt, value_fg, value_bg)
  feats_wdt <- feats_wdt[, .SD[1], .(feature, cell_type)]
  
  # if background cell types are not given, use all other cell types
  if (is.null(ct_bg)) {
    ct_bg <- setdiff(unique(feats_wdt$cell_type), ct_fg)
  }
  if (!any(ct_bg %in% feats_wdt$cell_type)) {
    warning("Some background cell types not present in the feature data: ", paste(setdiff(ct_bg, feats_wdt$cell_type), collapse = ", "))
  }
  if (!any(ct_fg %in% feats_wdt$cell_type)) {
    stop("Some foreground cell types not present in the feature data: ", paste(setdiff(ct_fg, feats_wdt$cell_type), collapse = ", "))
  }
  if (!any(ct_fg %in% ct_bg)) {
    # stop("Some cell types both in foreground and background: ", paste(ct_fg[ct_fg %in% ct_bg], collapse = ", "))
    ct_bg <- setdiff(ct_bg, ct_fg)
  }
  
  # if fg includes broad cell types, remove corresponding cell types from bg
  if (!is.null(level_dt)) {
    ct_bg <- setdiff(ct_bg, level_dt[broad_cell_type %in% ct_fg]$cell_type)
  }
  
  message(sprintf("Foreground: %s\nBackground: %s", paste(ct_fg, collapse = ", "), paste(ct_bg, collapse = ", ")))
  
  # exclude blacklisted peaks 
  if (!is.null(feats_blacklist)) {
    feats_wdt <- feats_wdt[!feature %in% feats_blacklist]
  }
  
  # make sure all features have a value for every cell type (fg and bg)
  feats_wdt <- feats_wdt[order(feature, cell_type)][cell_type %in% c(ct_fg, ct_bg)]
  feats_ndt <- CJ(unique(feats_wdt$feature), c(ct_fg, ct_bg))
  setnames(feats_ndt, c("feature", "cell_type"))
  feats_wdt <- merge.data.table(feats_ndt, feats_wdt, by = c("feature", "cell_type"), all.x = TRUE, sort = FALSE)
  feats_wdt[is.na(value_fg), value_fg := 0]
  feats_wdt[is.na(value_bg), value_bg := 0]
  
  # select peaks based on value in foreground in any of selected cell types
  feats_fg_dt <- feats_wdt[cell_type %in% ct_fg & value_fg > min_fg]
  feats_fg_dt[, n_selected_cell_types := .N, feature]
  feats_fg <- unique(feats_fg_dt$feature)
  message(sprintf(
    "Found %s features with %s > %s in any of the selected foreground cell types.",
    length(feats_fg), val_fg_column, min_fg
  ))
  message(sprintf(
    "Of those, %s features are in all of %s selected cell types",
    length(unique(feats_fg_dt[n_selected_cell_types==length(ct_fg)]$feature)),
    length(ct_fg)
  ))
  
  # subset peaks based on value in bg
  feats_bg_dt <- feats_wdt[feature %in% feats_fg & cell_type %in% ct_bg][order(-value_bg)][, .SD[1], feature]
  feats_fg_bg <- unique(feats_bg_dt[value_bg < max_bg]$feature)
  feats_fg_bg_dt <- feats_fg_dt[feature %in% feats_fg_bg]
  feats_fg_bg_dt[, n_selected_cell_types := .N, feature]
  message(sprintf(
    "Found %s features with %s > %s in selected foreground cell types and %s < %s in background cell types",
    length(feats_fg_bg),
    val_fg_column, min_fg,
    val_bg_column, max_bg
  ))
  message(sprintf(
    "Of those, %s features are in all of %s selected cell types",
    length(unique(feats_fg_bg_dt[n_selected_cell_types==length(ct_fg)]$feature)),
    length(ct_fg)
  ))
  
  # sort data
  feats_fg_bg_dt[, c("value_fg", "value_bg") := NULL]
  setorder(feats_fg_bg_dt, -n_selected_cell_types, feature)
  
  # return subset of input data with info 
  setnames(feats_fg_bg_dt, "feature", feat_column)  
  feats_res_dt <- merge.data.table(
    feats_dt, feats_fg_bg_dt,
    by = intersect(colnames(feats_dt), colnames(feats_fg_bg_dt)),
    all.y = TRUE, sort = FALSE
  )
  feats_res_dt[["feature"]] <- feats_res_dt[[feat_column]]
  feats_res_dt[, selected_cell_types := paste(sort(unique(.SD$cell_type)),collapse=","), feature]
  setorderv(feats_res_dt, c("n_selected_cell_types", feat_column), order = c(-1L, 1L))
  setcolorder(feats_res_dt, colnames(feats_dt))
  feats_res_dt[, feature := NULL]
  return(feats_res_dt)
}

# (4) trees ----------------------------------------------------------------------------------

#' @param fp_mat footprint matrix
#' @param fc_thr threshold for peak selection
parse_fps <- function(fp_mat, fc_thr = 2) {
  
  # subset variable peaks
  feats <- names(which(apply(
    fp_mat, 1, function(x) {
      max(x) > fc_thr
    }
  )))
  
  # selected features
  feats_select <- intersect(feats, rownames(fp_mat))
  
  # subset matrices
  table_fp <- fp_mat[feats_select, ]
  
  # distance matrix
  dmatrix <- dist(t(table_fp))
  
  return(dmatrix)
}

#' Peaks overlap
#' 
peaks_ovl <- function(fp_mat, fc_thr = 1.5) {
  jacc_mt <- matrix(
    NA,
    nrow = ncol(fp_mat),
    ncol = ncol(fp_mat),
    dimnames = list(
      colnames(fp_mat), colnames(fp_mat)
    )
  )
  for (ct_1 in colnames(fp_mat)) {
    for (ct_2 in colnames(fp_mat)) {
      pk1 <- names(which(fp_mat[, ct_1] > fc_thr))
      pk2 <- names(which(fp_mat[, ct_2] > fc_thr))
      jac <- length(intersect(pk1, pk2)) / length(union(pk1, pk2))
      jacc_mt[ct_1, ct_2] <- jac
    }
  }
  
  return(jacc_mt)
  
}

#' Cluster similar motifs
#' @param sim_mat motif similarity matrix
#' @param hclust clustering of motifs
#' @param hclust_method optionaly, cluster motifs on the go
siml_mot <- function(sim_mat, hclust = NULL, hclust_method = NULL, k = 1250) {
  ord <- rownames(sim_mat)
  if (!is.null(hclust_method)) {
    hclust <- hclust(dist(sim_mat), method = hclust_method)
  }
  cutree(hclust, k = k)
}

#' Parse motifs
#' @mot_dt motifs enrichment data frame
#' @param cluster_by either "fc" or "padj"
parse_mot <- function(mot_dt, ctr, cluster_by = "fc", padj_thrs = 0.05, fc_thrs = 1.5) {
  
  # select one motif per cluster
  mot_dt[, cluster := ctr[motif]]
  setorder(mot_dt, cluster, padj)
  if (cluster_by == "fc") {
    select_motifs <- mot_dt[padj < 0.05][fc != Inf][order(-fc)][, .SD[1], cluster]$motif
  } else if (cluster_by == "padj") {
    select_motifs <- mot_dt[fc != Inf][order(padj)][, .SD[1], cluster]$motif
  }
  mot_dt <- mot_dt[motif %in% select_motifs]
  
  # select significant motifs
  sg_mot <- mot_dt[padj < padj_thrs & fc > fc_thrs]$motif
  mot_dt <- mot_dt[motif %in% sg_mot]
  
  # transform data for clustering
  mot_dc <- dcast.data.table(
    mot_dt, motif ~ cell_type_label,
    value.var = cluster_by
  )
  mot_mt <- data.matrix(mot_dc[, -1])
  rownames(mot_mt) <- mot_dc$motif
  
  # distance matrix
  dmatrix <- tgs_dist(t(mot_mt))
  
  return(dmatrix)
}

#' Plots tree
#' @param dmatrix distance matrix
trees_funct <- function(dmatrix, groups, colors, rename_pattern = NULL) {
  tree <- nj(dmatrix)
  if (!is.null(rename_pattern))
    tree$tip.label <- str_replace_all(tree$tip.label, rename_pattern)
  tips <- tree$tip.label
  tree <- groupOTU(
    tree,
    sapply(
      groups,
      function(ct) grep(ct, tips),
      USE.NAMES = TRUE,
      simplify = FALSE
    )
  )
  gt <- ggtree(tree, aes(color = group), ladderize = FALSE) +
    scale_color_manual(values = colors) +
    theme_tree2() +
    geom_tiplab()
  xlim <- max(gt$data$x)
  if (xlim < 1) {
    xlim <- xlim + 2 * xlim
  } else {
    xlim <- xlim + 0.8 * xlim
  }
  gt <- gt +
    scale_x_continuous(limits = c(NA, xlim)) +
    theme(legend.position = "none")
  return(gt)
}


