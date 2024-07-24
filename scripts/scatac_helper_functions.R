require(data.table)
require(ggplot2)
require(ggrepel)
require(ggiraph)

accessibility_summary <- function(sc_matrix, sc_labels, fun=sum, nbins=10L) {
    tryCatch(
      t(apply(sc_matrix, 1,  function(x) tapply(x, sc_labels, function(y) fun(y)))),
      error = function(e) {
        warning(e)
        message("Calculating summary in ",nbins," bins")
        sum_list <- vector("list",nbins)
        sl <- split(1:nrow(sc_matrix),cut(1:nrow(sc_matrix),nbins))
        for (i in seq_along(sl)) {
          message(i," / ", nbins)
          sum_sub <- sc_matrix[sl[[i]],]
          sum_mat <- t(apply(sum_sub, 1,  function(x) tapply(x, sc_labels, function(y) fun(y))))
          sum_list[[i]] <- sum_mat 
        }
        do.call(rbind, sum_list)
      }
    )
}

accessibility_footprint <- function(sc_matrix, sc_labels, nbins=10L) {
  ct_geomean = tryCatch(
    t(apply(sc_matrix, 1,  function(x) tapply(x, sc_labels, function(y) exp(mean(log(1 + y))) - 1))),
    error = function(e) {
      warning(e)
      message("Calculating summary in ",nbins," bins")
      sum_list <- vector("list",nbins)
      sl <- split(1:nrow(sc_matrix),cut(1:nrow(sc_matrix),nbins))
      for (i in seq_along(sl)) {
        message(i," / ", nbins)
        sum_sub <- sc_matrix[sl[[i]],]
        sum_mat <- t(apply(sum_sub, 1,  function(x) tapply(x, sc_labels, function(y) exp(mean(log(1 + y))) - 1)))
        sum_list[[i]] <- sum_mat 
      }
      do.call(rbind, sum_list)
    }
  )
  ct_meansize = tapply(colSums(sc_matrix), sc_labels, mean)
  ideal_cell_size = pmin(1000,median(ct_meansize))
  g_fp = t(ideal_cell_size * t(ct_geomean) / as.vector(ct_meansize))
  fp_reg = 0.05
  g_fp_n = (fp_reg + g_fp) / apply(fp_reg + g_fp, 1, median)
  return(g_fp_n)
}

#' Read bed or bed-like formated file
#' @param fn character, filename 
#' @param format character, either "bed" or "narrowPeak"
#' @param column_names character, specify column names explicitly (otherwise they are inferred from `format`)
#' @param as_Granges logical, whether to return GRanges object (default: TURE), otherwise returns data.table
#' 
read_bed <- function(
  filename, 
  format = "bed",
  column_names = NULL, 
  as_GRanges = TRUE
) {
  if (is.null(column_names))
    column_names <- switch (
      format,
      "bed" = c("seqnames","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"),
      "narrowPeak" = c("seqnames","start","end","name","score","strand","signalValue","pValue","qValue","peak")
    )
  dt <- fread(filename, sep="\t")
  setnames(dt, column_names[1:ncol(dt)])
  if (as_GRanges==TRUE)
    dt <- makeGRangesFromDataFrame(as.data.frame(dt), keep.extra.columns=TRUE)
  return(dt)
}

read_gtf <- function(filename, as_GRanges = TRUE) {
  dt <- fread(filename)
  gtf_cols <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  setnames(dt, gtf_cols)
  if (as_GRanges==TRUE)
    dt <- makeGRangesFromDataFrame(dt, keep.extra.columns=TRUE)
  return(dt)
}

write_bed <- function(x, filename) {
  x <- as.data.frame(x, keep.extra.columns = TRUE)
  bed_cols <- c("seqnames", "start", "end", "score", "name", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
  bed_cols <- bed_cols[bed_cols %in% colnames(x)]
  x <- x[,bed_cols]
  fwrite(x, filename, sep="\t", quote=FALSE, col.names = FALSE)
}

write_gtf <- function(x, filename) {
  gtf_cols <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  x <- as.data.frame(x)[,gtf_cols]
  fwrite(x, filename, sep="\t", quote=FALSE, col.names = FALSE)
}

#' Plot scATAC 2D projection
#' 
#' @param df data.frame with two columns with cell coordinates, rownames should be cell names
#' @param color named vector of colors of same length as `nrow(df)`, names are groups (e.g. cell types)
#' 
plot_2d_proj <- function(
  df, 
  color=NULL,
  defaultColor = "#C4C4C4",
  highlightPoints = NULL,
  size = 1,
  xlim = NULL, 
  ylim = NULL, 
  extend = 0.05, 
  xlabel = "UMAP1", 
  ylabel = "UMAP2",
  title = "", 
  randomize = FALSE, 
  seed = 1,
  alpha = 1, 
  baseSize = 12, 
  showLabels = TRUE,
  showLegend = FALSE,
  legendTextSize = 12,
  legendSize = 4,
  legendPosition = "bottom",
  legendDirection = "vertical",
  legendRows = 3,
  ratioYX = 1, 
  bgWidth = 1,
  labelSize = 3,
  panelBorder = TRUE,
  pointBorder = FALSE,
  pointBorderColor = "black",
  labelBorderColor = "white",
  rastr = FALSE, 
  dpi = 300,
  interactive = FALSE,
  ...
) {
  
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

  if (!is.null(rownames(df))) {
    dnames <- rownames(df)
  } else {
    dnames <- ""
  }
  df <- data.frame(x = x, y = y, name = dnames)

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
      if(length(highlightPoints) < length(color)){
        color[-highlightPoints] <- defaultColor
        idx <- c(idx[-highlightPoints], idx[highlightPoints])
      }
    }
    
    df$color <- color
    if (is.null(names(color))) {
      names(color) <- color
      showLabels <- FALSE
    }
    
    df$group <- names(color)
    df$tooltip <- paste(df$group, df$name)
    group_values <- unique(names(color))
    pallete <- color[match(group_values, names(color))]
    names(pallete) <- group_values
    pallete <- pallete[sort(group_values)]
    df$group <- factor(df$group, levels = sort(group_values))   
    
    
  } else {
    
    df$color <- defaultColor
    df$group <- ""
    
  }
  
  # set up plot
  dp <- ggplot(df[idx,], aes(x = x, y = y, color=group)) + 
    coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
    xlab(xlabel) + ylab(ylabel) + 
    ggtitle(title) 
  
  # plot points
  # if (interactive == TRUE) {
    
    if (pointBorder == TRUE) {
      dp <- dp + geom_point_interactive(
        data = df[idx,], aes(x = x, y = y, fill = group, tooltip = tooltip, data_id = group), 
        color = pointBorderColor, pch = 21,
        size = size, alpha = alpha
      )
    } else {
      dp <- dp + geom_point_interactive(
        data = df[idx,], aes(x = x, y = y, color = group, tooltip = tooltip, data_id = group), 
        size = size, alpha = alpha
      )
    }

  # } else {
  #   
  #   if (pointBorder == TRUE) {
  #     dp <- dp + .geom_point_rast2(
  #       data = df[idx,], aes(x = x, y = y, fill = group), 
  #       color = pointBorderColor, pch = 21,
  #       size = size, raster.dpi = dpi, alpha = alpha, 
  #       raster.width = min(par('fin')),
  #       raster.height = (ratioYX * min(par('fin')))
  #     )
  #   } else {
  #     dp <- dp + .geom_point_rast2(
  #       size = size, raster.dpi = dpi, alpha = alpha, 
  #       raster.width = min(par('fin')),
  #       raster.height = (ratioYX * min(par('fin')))
  #     )
  #   }
  #   
  # }
  
  if (!is.null(highlightPoints)) {
    dp <- dp + .geom_point_rast2(data = df[highlightPoints,],
                                 size = size, raster.dpi = dpi, alpha = alpha, 
                                 raster.width = min(par('fin')),
                                 raster.height = (ratioYX * min(par('fin')))
    )
  }
  
  dp <- dp + 
    scale_color_manual(values = pallete) + 
    scale_fill_manual(values = pallete) + 
    theme_classic() + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(size=1, fill=NA)
    ) 
  
  if (panelBorder == FALSE) {
    dp <- dp + theme(
      panel.border = element_blank()
    )
  }
  
  if (showLegend==FALSE) {
    dp <- dp + theme(
      legend.position = "none"
    )
  } else {
    dp <- dp + theme(
      legend.position = legendPosition,
      legend.title = element_blank(), 
      legend.text = element_text(size = legendTextSize),
      legend.direction = legendDirection
    ) + guides(colour = guide_legend(nrow=legendRows, override.aes = list(size = legendSize, alpha = 1)))
  }
  
  if (showLabels==TRUE) {
    dt <- as.data.table(df)
    dtlab <- dt[,lapply(.SD,median),by=group,.SDcols=c("x","y")]
    dp <- dp + geom_text_repel(
      data=dtlab, 
      aes_string("x", "y", label="group"), 
      box.padding = 0.5, 
      max.overlaps = Inf,
      bg.color = labelBorderColor, 
      bg.r = 0.15,
      segment.color = NA
    )
  }
  
  if (interactive == TRUE) {
    girafe(
      ggobj = dp,
      width_svg = 8, height_svg = 8,
      options = list(
        opts_selection(type = "none", only_shiny = FALSE),
        opts_hover_inv(css = "opacity:0.1;"),
        opts_hover(css = "opacity:1;")
      )
    )
  } else {
    dp
  }
}

#' Plot gene accessibility on UMAP
#' 
#' @param df data.frame with two columns with UMAP coordinates, rownames should be cell names
#' 
plot_2d_proj_gene <- function(
  df, 
  sc_values,
  color_scale = c("gray95","lightyellow","khaki1","orange","orangered2","#520c52"),
  name = "accessibility",
  labelMeans = TRUE,
  highlightPoints = NULL,
  size = 1,
  xlim = NULL,
  ylim = NULL,
  extend = 0.05,
  xlabel = "UMAP1",
  ylabel = "UMAP2",
  title = "",
  legend.position = c(0.9, 0.1),
  randomize = FALSE,
  seed = 1,
  sc_min = NULL,
  sc_max = NULL,
  alpha = 1,
  pointBorder = FALSE,
  pointBorderColor = "black",
  panelBorder = TRUE,
  baseSize = 10,
  legendSize = 3,
  ratioYX = 1,
  rastr = FALSE,
  interactive = TRUE,
  dpi = 300,
  na_color = "white",
  ...
) {
  
  x=df[,1]
  y=df[,2]
  
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }
  
  name <- str_replace_all(name, " ", "_")

  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(sc_values)){
      sc_values <- sc_values[include]
    }
  }
  
  if(is.null(xlim)){
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  
  if(is.null(ylim)){
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  
  # sc values
  if (is.null(sc_min))
    sc_min = quantile(sc_values[!is.infinite(sc_values)&!is.na(sc_values)], 0.03, na.rm = TRUE)
  message(sprintf("min: %s",sc_min))
  if (is.null(sc_max))
    sc_max = quantile(sc_values[!is.infinite(sc_values)&!is.na(sc_values)], 0.97, na.rm = TRUE)
  message(sprintf("max: %s",sc_max))
  
  # apply min/max to sc vector
  sc_values [ sc_values < sc_min ] = sc_min
  sc_values [ sc_values > sc_max ] = sc_max
  
  df[, name] <- sc_values
  
  # order for plotting
  if (randomize) {
    set.seed(seed)
    idx <- sample(seq_along(sc_values), length(sc_values))
  } else {
    idx <- order(sc_values,decreasing=FALSE)
  }
  
  if (pointBorder) { 
    gp <- ggplot(df[idx,], aes_string(x = "x", y = "y", fill=name)) +
      scale_fill_gradientn(colours = color_scale, na.value = na_color) 
    pch <- 21
  } else {
    gp <- ggplot(df[idx,], aes_string(x = "x", y = "y", color=name)) +
      scale_color_gradientn(colours = color_scale, na.value = na_color)
    pch <- 16
  }
  gp <- gp + 
    coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) + 
    xlab(xlabel) + ylab(ylabel) + 
    ggtitle(title) +
    theme_classic() + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      legend.title = element_text(vjust=0.5,hjust=0.5),
      legend.direction = "horizontal",
      legend.position = legend.position,
      panel.border = element_rect(size=1, fill=NA)
    )
  
  if (!panelBorder) {
    gp <- gp + theme(panel.border = element_blank())
  }

  if (rastr) {
    
    gp + .geom_point_rast2(
      color = pointBorderColor, pch = pch,
      size = size, raster.dpi = dpi, alpha = alpha, 
      raster.width = min(par('fin')), 
      raster.height = (ratioYX * min(par('fin')))
    )
    
    
  } else {
    
    if (interactive) {
      
      gp <- gp + geom_point_interactive(aes(tooltip = accessibility), pch = pch)
      girafe(ggobj = gp)
      
    } else {
      
      gp + geom_point(size = size, alpha = alpha, pch = pch, )
      
    }
    
  }
}


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

# GO analysis

require(topGO)

#' @param list_interest character, gene names
#' @param gomap named list of GO annotations for genes, names of list should be gene names
#' @param output_name character, prefix for output file names
#' @param name_geneset character, used to construct file name for output files
#' @param ontology_set character(s) indicating GO ontology to use, `c("BP","CC","MF")` 
#' @param tg_test character, which test to use, one of `c("fisher","t")`, see `statistic` in `?topGO::runTest`
#' @param tg_algorithm character, which algorithm to use, see `algorithm` in `?topGO::runTest`
#' @param printfile logical, whether to save plot and table
#' @param p_adj_method character, multiple correction method to use
topgofun  <- function(
  list_interest, gomap, output_name, name_geneset, ontology_set, tg_test="fisher", tg_algorithm="classic", 
  topnum=20, nodesize=10, printfile=TRUE, p_adj_method="BH", firstSigNodes=10
) {
  
  library(topGO)
  
  # Input 
  list_interest = unique(list_interest)
  genom = names(gomap)
  gesel = factor(as.integer(genom %in% list_interest))
  names(gesel) = genom
  
  # shortened go mappings without empty transcripts
  gomap_nonempty = gomap[lapply(gomap,length)>0]
  
  namepref <- paste0(output_name,".",name_geneset,".topgo",".",tg_test,tg_algorithm)
  # if(printfile){
  #   pdf(file=paste0(namepref,".pdf"),height=4.5,width=4)
  # }
  par(mar=c(5,12,5,2))
  
  topgo_tau_tot = data.frame()
  
  if (length(list_interest[list_interest %in% names(gomap_nonempty)])>1) {
    
    for (ontology_seti in ontology_set) {
      # topGO setup 
      
      GOdata = new(
        "topGOdata", ontology=ontology_seti, allGenes=gesel,
        annot=annFUN.gene2GO, gene2GO=gomap
      )
      
      num_interest_feasible = sum(GOdata@feasible & genom %in% list_interest)
      
      # topGO analysis
      topgo_res = runTest(GOdata, algorithm = tg_algorithm, statistic = tg_test)
      topgo_tau = GenTable(
        GOdata, pval_test = topgo_res, orderBy = "pval_test", 
        topNodes = length(usedGO(object = GOdata))
      )
      topGO::printGraph(
        GOdata, result=topgo_res, firstSigNodes=firstSigNodes, # all the nodes in the graph: length(usedGO(object = GOdata)) -- a mess
        useInfo="all", fn.prefix=paste(namepref,ontology_seti,sep="."), pdfSW=TRUE
      )
      topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
      topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
      topgo_tau$ontology = ontology_seti
      topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
      
      # Output 
      # ploti=barplot(height = rev(head(log(topgo_tau$pval_test,10),topnum)),
      #               names.arg = rev(head(paste(topgo_tau$Term,topgo_tau$GO.ID),topnum)),
      #               xlim=c(0,-5),horiz=T,las=1,col="slategray3",border=NA,
      #               cex.names=0.35,cex.axis=0.6,cex.lab=0.6,cex.sub=0.6,cex.main=0.6,
      #               main=paste(name_geneset,"top GO:",ontology_seti,tg_test,tg_algorithm),
      #               sub =paste("n=",num_interest_feasible,"/",length(list_interest), sep=""),
      #               xlab="log(p)")
      # abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
      # abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
      # text(x=0,ploti,labels = paste("p =",signif(rev(head(topgo_tau$pval_test,topnum)),3)),
      #      col="red",pos=4,cex=0.35)
    }
    
  }else {
    print("skip, no annotations in interest list!")
  }
  
  if(printfile){
    write.table(
      topgo_tau_tot,
      file=paste(output_name,".",name_geneset,".topgo",".",tg_test,tg_algorithm,".txt",sep=""),
      sep="\t", quote=F, col.names=T, row.names=F, append = F)
    dev.off()
  }
  
  return(topgo_tau_tot)
}

# Revigo

plot_revigo <- function(
  df,
  legend_position = "bottom",
  legend_box = "vertical"
) {

    pdt <- copy(df)
    pdt[, plot_X := as.numeric(as.character(pdt$PC_0))]
    pdt[, plot_Y := as.numeric(as.character(pdt$PC_1))]
    pdt <- pdt[!is.na(plot_X) & !is.na(plot_Y)]
    pdt[, num_annots := as.numeric(as.character(pdt$LogSize))]
    pdt[, log10padj := as.numeric(as.character(pdt$Value)) * -1]
    pdt[, frequency := as.numeric(as.character(pdt$Frequency))]
    pdt[, uniqueness := as.numeric(as.character(pdt$Uniqueness))]
    pdt[, dispensability := as.numeric(as.character(pdt$Dispensability))]
    class(pdt) <- "data.frame"
    
    ex <- pdt[pdt$dispensability < 0.15, ]
    x_range <- max(pdt$plot_X) - min(pdt$plot_X)
    y_range <- max(pdt$plot_Y) - min(pdt$plot_Y)
    p_range <- max(pdt$log10padj)
    
    p1 <- ggplot(pdt) +
        geom_point(
            aes(plot_X, plot_Y, fill = log10padj, size = num_annots),
            shape = 21,
            alpha = 1
        ) +
        scale_fill_gradientn(
          name = "log10\nadjusted pvalue",
          colours = c(RColorBrewer::brewer.pal(4, "Blues")[-1], "#012d66", "#01153d"),
          limits = c(0, p_range)
        ) +
        scale_size(
          name = "number of\nannotations",
          range = c(2, 10)
        ) +
        geom_text_repel(
            data = ex,
            aes(plot_X, plot_Y, label = Name),
            colour = I(alpha("black", 0.85)),
            size = 4
        ) +
        labs(y = "MDS2", x = "MDS1") +
        coord_fixed() +
        xlim(
            min(pdt$plot_X) - x_range / 10,
            max(pdt$plot_X) + x_range / 10
        ) +
        ylim(
            min(pdt$plot_Y) - y_range / 10,
            max(pdt$plot_Y) + y_range / 10
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = legend_position,
            legend.box = legend_box
        )

    p1
}

