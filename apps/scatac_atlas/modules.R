atac2DUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box(
        title="Atlas", width=12, height=1500, solidHeader = TRUE,
        selectInput(
          inputId = ns("colorby"), 
          label = "Color 2D projection by", 
          choices = list("Cell type" = "cell_type", "Broad cell type" = "broad_cell_type"),
          selected = "cell_type",
          selectize = TRUE
        ),
        tryCatch(girafeOutput(ns("plot_2d_proj"), height=800), error=function(e) warning(e)),
        tryCatch(girafeOutput(ns("plot_bar"), height=400), error=function(e) warning(e))
        
      ),
      box(
        title="Gene accessibility", width=12, height=1400, solidHeader = TRUE, 
        selectInput(
          inputId=ns("searchid"), 
          label="Search gene by name or id", 
          choices=NULL,
          selectize = TRUE
        ),
        tryCatch(girafeOutput(ns("plot_2d_gene_proj"), height=800), error=function(e) warning(e)),
        tryCatch(girafeOutput(ns("plot_box_gene"), height=350), error=function(e) warning(e)),
        br(),
        tableOutput(ns("single_gene"))
        # column(width = 3, sliderInput(ns("sc_min"), label = "Scale color to minimum value:", min=0, max=5, value=0, step=0.1, width="100%")),
        # column(width = 3, sliderInput(ns("sc_max"), label = "Scale color to maximum value:", min=0, max=5, value=2.5, step=0.1, width="100%"))
      ),
      box(
        title = "Gene summary", width = 12, height = 800, solidHeader = TRUE,
        fluidRow(column(width = 12, 
                        "The table below summarizes accessibility scores for selected gene in cells belonging to the same annotation group (cell type, broad cell type, or sample).",
                        "For each group, ncells is te total number of cells, ncells_access is the number of cells and percent_acess is the fraction of cells in which a gene is acessible.",
                        "A gene is considered acessible if its score is above the selected threshold.")),
        br(),
        fluidRow(
          column(width = 3, selectInput(ns("acces_group"), label = "Group cells by", list("Cell type" = "cell_type", "Broad cell type" = "broad_cell_type"), selected = "cell_type", selectize = TRUE)),
          column(width = 3, sliderInput(ns("acces_thrs"), label = "Accessibility threshold", min=0, max=5, value=1, step=0.5))
        ),
        fluidRow(column(width=12, dataTableOutput(ns("single_gene_summary"))))
      )
    )
  )
  
}

atac2DServer <- function(id, config_file="config.yaml", config_id) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir
      
      # load metadata
      sc_METADATA <- fread(file.path(INPUT_DIR, conf[[config_id]]$meta_file_sc))
      class(sc_METADATA) <- "data.frame"
      rownames(sc_METADATA) <- sc_METADATA$Cell
      cells <- rownames(sc_METADATA)
      
      mc_METADATA <- fread(file.path(INPUT_DIR, conf[[config_id]]$meta_file_mc))
      class(mc_METADATA) <- "data.frame"
      rownames(mc_METADATA) <- mc_METADATA$SEACell
      SEACells <- rownames(mc_METADATA)
      
      # load 2D projections
      mc_UMAP <- fread(file.path(INPUT_DIR, conf[[config_id]]$umap_file_mc))
      mc_UMAP <- as.data.frame(mc_UMAP)
      rownames(mc_UMAP) <- mc_UMAP$SEACell
      mc_UMAP <- mc_UMAP[,2:3]
      colnames(mc_UMAP) <- c("x","y")
      
      # load gene scores 
      genescore_mat <- readRDS(file.path(INPUT_DIR, conf[[config_id]]$genescore_file_mc))
      SEACells <- SEACells[SEACells %in% colnames(genescore_mat)]
      mc_METADATA <- mc_METADATA[SEACells, ]
      mc_UMAP <- mc_UMAP[SEACells,]
      
      # plot 2D projection
      plot_umap_f <- reactive({
        showLeg <- TRUE
        showLab <- FALSE
        groups <- as.character(unique(sc_METADATA[,input$colorby]))
        legRows <- round(length(groups))
        
        col_pal_mc <- structure(mc_METADATA[,paste0(input$colorby,"_color")], names=mc_METADATA[,input$colorby])
        groups_names <- as.character(mc_METADATA[SEACells,input$colorby])
        groups_colors <- col_pal_mc[groups_names]
        p_mc <- plot_umap(
          df = as.data.frame(mc_UMAP), 
          color = groups_colors, 
          groups = groups,
          downsample = FALSE,
          pointBorder = TRUE,
          size = 4, 
          alpha = 1,
          xlabel = "UMAP1", 
          ylabel = "UMAP2",
          showLegend = showLeg,
          showLabels = showLab,
          legendRows = legRows
        ) + theme(legend.position = "right")
        
        # combine plots
        girafe(
          code = print(p_mc),
          width_svg = 12, height_svg = 9,
          options = list(
            opts_selection(type = "none", only_shiny = FALSE),
            opts_hover_inv(css = "opacity:0.1;"),
            opts_hover(css = "opacity:1;")
          )
        )
      })
      
      output$plot_2d_proj <- renderGirafe(plot_umap_f())
      
      plot_bar_f <- reactive({
        
        req(input$colorby)
        
        # barplot per sample
        groups <- as.character(unique(sc_METADATA[,input$colorby]))
        legRows <- ceiling(length(groups)/10)
        sample_dt <- sc_METADATA[,c(input$colorby, "Sample")]
        sample_vc <- intersect(samples, sample_dt$Sample)
        
        p_bar_count <- plot_gene_access_bar(
          df = sample_dt, 
          x_order = groups,
          group_column = "Sample", 
          group_order = sample_vc,
          col_pal = orig_cols[sample_vc], 
          showLabels = TRUE,
          showCounts = TRUE,
          barPosition = "stack",
          legendPosition = "bottom", 
          legendRows = legRows
        ) + labs(y = "cells")
        
        # p_bar_frac <- plot_gene_access_bar(
        #   df = sample_dt, 
        #   x_order = groups,
        #   group_column = "Sample", 
        #   group_order = sample_vc,
        #   col_pal = orig_cols[sample_vc], 
        #   showLabels = TRUE,
        #   barPosition = "fill",
        #   legendPosition = "none", 
        #   legendRows = 2
        # )
        
        # combine plots
        # combined_plot <- patchwork::wrap_plots(list(p_bar_count, p_bar_frac), nrow = 2)
        
        girafe(
          code = print(p_bar_count),
          width_svg = 12, height_svg = 5,
          options = list(
            opts_selection(type = "none", only_shiny = FALSE),
            opts_hover_inv(css = "opacity:0.1;"),
            opts_hover(css = "opacity:1;")
          )
        )
        
      })
      
      output$plot_bar <- renderGirafe(plot_bar_f())
      
      # gene search - subset only genes with accessibility
      genes <- rownames(genescore_mat)
      GENE_ANNT_ATAC <- GENE_ANNT[gene_id %in% genes]
      
      updateSelectizeInput(
        session, "searchid",
        choices = GENE_ANNT_ATAC$search_id,
        selected = grep("Pou4 ", GENE_ANNT_ATAC$search_id, value = TRUE)[1],
        server = TRUE
      )
      
      # UMAP of gene accessibility
      plot_gene_umap_f <- reactive({
        selected_gene <- GENE_ANNT_ATAC[search_id==input$searchid,gene_id]
        mc_vals <- as.numeric(genescore_mat[match(selected_gene,genes),SEACells])
        names(mc_vals) <- SEACells
        p_mc <- tryCatch(
          plot_umap_gene(
            df = mc_UMAP,
            group_values = sc_METADATA[SEACells,input$colorby],
            color_scale = proj_colors,
            sc_values = mc_vals,
            downsample = FALSE,
            printBorder = TRUE,
            size = 4,
            alpha = 1,
            xlabel = "UMAP1",
            ylabel = "UMAP2",
            legend.position = "right"
          ),
          error = function(e) warning(e)
        )
        girafe(
          code = print(p_mc),
          width_svg = 12, height_svg = 9,
          options = list(
            opts_selection(type = "none", only_shiny = FALSE),
            opts_hover_inv(css = "opacity:0.1;"),
            opts_hover(css = "opacity:1;")
          )
        )
      })
      
      # boxplot of gene accessibility
      plot_gene_box_f <- reactive({
        selected_gene <- GENE_ANNT_ATAC[search_id==input$searchid,gene_id]
        mc_vals <- as.numeric(genescore_mat[match(selected_gene,genes),SEACells])
        names(mc_vals) <- SEACells
        access_dt <- data.table(group = mc_METADATA[SEACells,input$colorby], accessibility = mc_vals)
        col_pal_mc <- structure(mc_METADATA[, paste0(input$colorby, "_color")], names = mc_METADATA[,input$colorby])
        p_box <- plot_gene_access_box(
          access_dt, 
          group_column = "group", 
          col_pal = col_pal_mc, 
          showLabels = TRUE
        )
        girafe(
          code = print(p_box),
          width_svg = 10, height_svg = 4,
          options = list(
            opts_selection(type = "none", only_shiny = FALSE),
            opts_hover_inv(css = "opacity:0.1;"),
            opts_hover(css = "opacity:1;")
          )
        )
      })
      
      # combine plots
      # combined_gene_plot <- p_mc / p_box + plot_layout(nrow = 2, heights = c(4,1)) 
      # combined_gene_plot <- patchwork::wrap_plots(list(p_mc, p_box), nrow = 2, heights = c(4,1))
      
      output$plot_2d_gene_proj <- renderGirafe({
        shiny::req(input$searchid)
        plot_gene_umap_f()
      })
      
      output$plot_box_gene <- renderGirafe({
        shiny::req(input$searchid)
        plot_gene_box_f()
      })
      
      # mini table showing gene name and id
      output$single_gene <- renderTable(
        GENE_ANNT_ATAC[search_id==input$searchid, 1:(ncol(GENE_ANNT_ATAC)-1)],
        width = "100%"
      )
      
      # gene summary table
      atacgsm <- reactive({
        req(input$searchid)
        selected_group <- input$acces_group
        selected_gene <- GENE_ANNT_ATAC[search_id==input$searchid,gene_id]
        mc_vals <- as.numeric(genescore_mat[match(selected_gene,genes),SEACells])
        names(mc_vals) <- SEACells
        wdt <- copy(mc_METADATA)
        setDT(wdt)
        wdt[,acessibility:=mc_vals]
        wdt[,nSEACells:=.N, by=selected_group]
        wdt[,access:=mean(acessibility), by=selected_group]
        wdts <- wdt[,.(nSEACells_access = sum(acessibility > input$acces_thrs)), .(get(selected_group),nSEACells,access)]
        wdts[,`percent_access`:=nSEACells_access/nSEACells]
        setnames(wdts, "get", "group")
        setcolorder(wdts, c("group","access","nSEACells","nSEACells_access","percent_access"))
        setorder(wdts,-access)
        wdts
      })
      
      output$single_gene_summary <- DT::renderDataTable(
        datatable(atacgsm()),
        rownames = FALSE,
        options = list(
          dom = 'tp', scrollX = TRUE, ordering = TRUE, pageLength = 20,
          columnDefs = list(list(className = 'dt-center', targets = 0:2))
        ),
        selection = list(mode = 'none')
      )
      
    }
  ) 
  
}

multiGeneUI <- function(id) {
  ns <- NS(id)
  
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title="Gene selection", width=12, height=NULL, solidHeader=TRUE,
        shiny::h5("Choose genes for which to plot the expression across metacells.
              Either search for genes by typing (part of) the name or gene id in the search bar,
              or upload a text file with genes."),
        shiny::h5("Select genes to show on the heatmap either by clicking on the individual rows
              of the search results table followed by 'Add selected genes', or just by clicking
              'Add all genes' to include all the search results."
        ), br(),
        radioButtons(
          ns("geneselecttype"), label = "Choose genes",
          choices = list("search for genes" = "search",  "upload list of genes" = "upload"),
          selected = "search"
        ),
        conditionalPanel(
          condition = "input.geneselecttype == 'search'", ns = ns,
          searchInput(
            inputId = ns("free_genes"),
            label = "Search for gene:",
            placeholder = "Type gene name or ID",
            btnSearch = icon("search"),
            btnReset = icon("remove"),
            width = "50%"
          )
        ),
        conditionalPanel(
          condition = "input.geneselecttype == 'upload'", ns = ns,
          shiny::h5("Upload a text text file with genes. Each gene should be on the new line."),
          fileInput(
            ns("genefile"), "Choose gene file",
            multiple = FALSE,
            accept = c("text/tsv","text/csv","text/tab-separated-values,text/plain","text/comma-separated-values,text/plain")
          )
        ),
        DTOutput(ns("geneSelectDT")),
        hr(),
        actionButton(ns("add_genes"), "Add selected genes"),
        actionButton(ns("add_all_genes"), "Add all genes")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="Choosen genes:", width=12, solidHeader=TRUE,
        DTOutput(ns("selected_genes_dt")), hr(),
        actionButton(ns("clear_genes"), "Clear selection"),
        downloadButton(ns("download_genes_data"),"Download table")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="Heatmap", width=12, solidHeader=TRUE,
        selectInput(
          ns("mcselecttype"), label = "Choose annotations from:",
          choices = list("Cell type" = "cell_type", "Broad cell type" = "broad_cell_type"), 
          selected = "cell_type"
        ),
        tags$style(".btn-custom {background-color: #b8b8b8; color: #FFF;}"),
        dropdownButton(
          prettySwitch(ns("clustergenes"), "Cluster genes", value=TRUE, inline=FALSE),
          sliderInput(
            inputId = ns("min_expression_fc"), label="Filter genes by min FC:",
            min=0, max=5, step=0.1, round=FALSE, value=0, width = "100%"
          ),
          sliderInput(
            inputId = ns("scale_expression_fc"), label="Scale to max FC value:",
            min=1, max=5, step=0.1, round=FALSE, value=5, width = "100%"
          ),
          sliderInput(
            inputId=ns("mcid_font_size"),
            label = "Metacell lables size",
            min = 0, max = 12, step = 1, value = 10
          ),
          sliderInput(
            inputId=ns("gene_font_size"),
            label = "Gene lables size",
            min = 0, max = 12, step = 1, value = 10
          ),
          circle = TRUE, status = "custom", icon = icon("gear"), width = "300px",
          tooltip = tooltipOptions(title = "Click to customize")
        ),
        withSpinner(
          uiOutput(ns("ui_genes_heatmap")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        ),
        br(),
        downloadButton(ns("download_genes_hmap"),"Download heatmap")
      )
    )
  )
  
}

multiGeneServer <- function(id, config_file="config.yaml", config_id) {
  
  shiny::moduleServer(
    id,
    
    function(input, output, session) {
      
      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir
      
      # gene annotation
      gene_annotation <-  conf[['default']]$gene_annot_file
      GENE_ANNT <- fread_gene_annotation(
        file = file.path(INPUT_DIR, gene_annotation),
        cols = c(1:3),
        colnames = c("gene_id","gene name","PFAM domain")
      )
      
      # load metadata
      sc_METADATA <- fread(file.path(INPUT_DIR, conf[[config_id]]$meta_file_sc))
      class(sc_METADATA) <- "data.frame"
      rownames(sc_METADATA) <- sc_METADATA$Cell
      cells <- rownames(sc_METADATA)
      
      mc_METADATA <- fread(file.path(INPUT_DIR, conf[[config_id]]$meta_file_mc))
      class(mc_METADATA) <- "data.frame"
      rownames(mc_METADATA) <- mc_METADATA$SEACell
      SEACells <- rownames(mc_METADATA)
      
      # gene search - subset only genes with accessibility
      genes <- lfile_mc$get.attribute.df(MARGIN = 1, attribute.names = "gene")[,1]
      GENE_ANNT <- GENE_ANNT[gene_id %in% genes]
      
      # Select genes
      genes_dt <- reactive({
        
        if (input$mcselecttype == "cell_type") {
          fpfile <- conf[[config_id]]$footprint_file_genes
        } else if (input$mcselecttype == "broad_cell_type") {
          fpfile <- conf[[config_id]]$footprint_file_genes_broad
        } 
        print(fpfile)
        FP <-readRDS(file.path(INPUT_DIR,fpfile))
        
        if (input$geneselecttype == "search") {
          df <- genes_select_dt(sterm = input$free_genes, nmat = FP, annt = GENE_ANNT)
        } else if (input$geneselecttype == "upload") {
          req(input$genefile)
          gs <- fread(input$genefile$datapath, header=FALSE)[[1]]
          df <- genes_upload_dt(gs=gs, annt=GENE_ANNT)
        }
        message(nrow(x),"genes")
        df
      })
      
      output$geneSelectDT <- renderDT(
        genes_dt(), rownames = FALSE,
        options = list(
          dom = 'tp', scrollX = TRUE, ordering = FALSE, pageLength = 10,
          columnDefs = list(list(className = 'dt-center', targets = 0:2))
        )
      )
      selected_genes <- reactiveValues()
      selected_genes$values <- c()
      observeEvent(
        input$add_genes,
        tryCatch({
          g <- genes_select_names(dt=genes_dt(), rid=input$geneSelectDT_rows_selected)
          selected_genes$values <- c(selected_genes$values, g)
        }, error = function(e) message("Select at least two genes!"))
      )
      observeEvent(
        input$add_all_genes,
        tryCatch({
          g <- genes_dt()$gene_id
          selected_genes$values <- c(selected_genes$values, g)
        }, error = function(e) message("There should be at least two genes!"))
      )
      observeEvent(input$clear_genes, {
        selected_genes$values <- c()
      })
      genes_selected <- reactive({
        rid <- unique(match(selected_genes$values,GENE_ANNT$gene_id))
        tryCatch(GENE_ANNT[rid, 1:3], error=function(e) NULL)
      })
      output$selected_genes_dt <- renderDT(
        if (!is.null(selected_genes$values))
          genes_selected(),
        rownames = FALSE,
        options = list(
          dom = 'tp', scrollX = TRUE, ordering = FALSE, pageLength = 10,
          columnDefs = list(list(className = 'dt-center', targets = 0:2))
        ),
        selection = list(mode = 'none')
      )
      output$download_genes_data <- downloadHandler(
        filename = function() {
          "selected_genes.tsv"
        },
        content = function(file) {
          write.table(genes_selected(), file, sep = "\t", quote = FALSE)
        }
      )
      
      # Heatmap of selected genes
      plotting_f <- function() {
        message("Plotting function")
        tryCatch(mgenes_hmap(
          nmat=FP, annt=GENE_ANNT, gids=selected_genes$values,
          min_expression_fc=input$min_expression_fc,
          scale_expression_fc=pmax(input$scale_expression_fc,input$min_expression_fc),
          cluster_genes=input$clustergenes,
          heatmap_colors=heatmap_colors,
          ct_table=CELL_ANNT,
          mcid_font_size=input$mcid_font_size,
          gene_font_size=input$gene_font_size
        ), error = function(e) message("Select at least two genes!"))
      }
      output$genes_hmap <- shiny::renderPlot({
        if (!is.null(selected_genes$values) & length(selected_genes$values)>1) {
          tryCatch(print(plotting_f()), error=function(e) NULL)
        }
      })
      hmh <- reactiveValues()
      hmap_height <- function() {
        hmh$ng <- mgenes_hmap_height(
          nmat = MCFP, gids = selected_genes$values, annt=GENE_ANNT,
          min_expression_fc=input$min_expression_fc
        )
        if (hmh$ng<5) {
          sf <- 80
        } else if (hmh$ng<15) {
          sf <- 50
        } else if (hmh$ng<25) {
          sf <- 20
        } else {
          sf <- 15
        }
        return(hmh$ng*sf)
      }
      output$ui_genes_heatmap <- renderUI({
        ns <- session$ns
        if (!is.null(selected_genes$values))
          shiny::plotOutput(ns("genes_hmap"), height = hmap_height(), width = "100%")
      })
      output$download_genes_hmap <- downloadHandler(
        filename = "selected_genes_heatmap.png",
        content = function(file) {
          png(file,width=2400,height=hmap_height(),res=110)
          print(plotting_f())
          dev.off()
        }
      )
      
      
    })
}

geneSummaryUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box(
        title="1. Annotation selection", width=4, solidHeader = TRUE,
        h5("Find genes that are specifically acessible in a group of cells."),
        h5("First select level of clustering (Clusters or Annotation), then select one or more groups of cells
           for which the gee acessibility will be calculated."),
        br(),
        selectInput(
          ns("mcselecttype"), label = "Choose annotations from:",
          choices = list("Cell type" = "cell_type", "Broad cell type" = "broad_cell_type"), 
          selected = "cell_type"
        ),
        selectInput(
          inputId=ns("ids_mcs"), 
          label="Select one or more cell groups", 
          choices=NULL, multiple=TRUE, selected=list(1),
          selectize = TRUE
        )
      ),
      box(
        title="2. Fold change (FC) selection", width=4, solidHeader = TRUE,
        h5("Select minimum FC threshold to be used for filtering genes. Genes that have FC 
               above this value in all selected groups (if using absolute method) or median FC 
               above this value in all selected groups (if using median method) will be shown in the 
               summary table."),
        radioButtons(
          ns("fc_sum_method"),"Select method for summarizing gene acessibility in selected groups of cells:", 
          choices = c("absolute","median"), selected = "median", inline = TRUE
        ),
        sliderInput(
          inputId = ns("fc_selection"), label="Choose minimum FC for selected groups of cells:", 
          min=1.0, max=5, step=0.1, round=FALSE, value=2
        ),
        h5("It is also possible to set a maximum FC threshold for non-selected (background) groups of cells. 
               In this case, the summary will only include genes that have FC above minimum threshold 
               in selected groups of cells, and below maximum threshold in all other groups."),
        checkboxInput(ns("fcbg"), "Choose maximum background FC threshold", value = FALSE),
        conditionalPanel(
          condition = "input.fcbg", ns=ns,
          radioButtons(
            ns("fcbg_sum_method"),"Select method for summarizing gene acessibility in non-selected groups of cells:", 
            choices = c("absolute","median"), selected = "absolute", inline = TRUE
          ),
          sliderInput(
            inputId = ns("fcbg_selection"), label="Choose maximum FC for non-selected groups of cells", 
            min=1.0, max=5, step=0.1, round=FALSE, value=2
          )
        )
      ),
      box(
        title="3. Generate summary", width=4, solidHeader=TRUE,
        h5("Click to generate summary using specified parameters."),
        actionButton(ns("dosummary"), "Go!"),
        h4("Your search: "),
        tableOutput(ns("summary"))
      )
    ),
    fluidRow(
      box(
        title="Gene table", width=12, solidHeader=TRUE,
        DT::DTOutput(ns("genes_summary_table")), shiny::br(),
        h5("Acessibility sum is the total summed acessibility of a gene in selected groups of cells
           (of note, it will be correlated with the number of cells in selected groups)."),
        h5("Acessibility fraction is the percentage of summed acessibility for a gene that come from the selected groups of cells."),
        h5("Median FC is the median enrichment of acessibility for a gene in selected groups of cells."),
        shiny::br(),
        downloadButton(outputId=ns("download_table"), label="Download gene table")
      )
    )
  )
}

geneSummaryServer <- function(id, config_file="config.yaml", config_id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir
      
      # load metadata
      sc_METADATA <- fread(file.path(INPUT_DIR, conf[[config_id]]$meta_file_sc))
      class(sc_METADATA) <- "data.frame"
      rownames(sc_METADATA) <- sc_METADATA$Cell
      cells <- rownames(sc_METADATA)
      
      mc_METADATA <- fread(file.path(INPUT_DIR, conf[[config_id]]$meta_file_mc))
      class(mc_METADATA) <- "data.frame"
      rownames(mc_METADATA) <- mc_METADATA$SEACell
      SEACells <- rownames(mc_METADATA)
      
      # update scluster selection options
      annts <- reactiveValues()
      observeEvent(input$mcselecttype, {
        annts$opts <- sort(unique(sc_METADATA[,input$mcselecttype]))
        updateSelectizeInput(
          session, 
          "ids_mcs",
          choices = as.character(unique(annts$opts)),
          server = TRUE
        )
      })
      
      # select footprint
      FP <- reactive({
        if (input$mcselecttype=="cell_type") {
          fpfile <- conf[[config_id]]$footprint_file_genes
        } else if (input$mcselecttype=="broad_cell_type") {
          fpfile <- conf[[config_id]]$footprint_file_genes_broad
        } 
        mat <- readRDS(file.path(INPUT_DIR,fpfile))
        colnames(mat)[colnames(mat) == "neuronal_gastrula"] <- "neuronal"
        mat
      })
      ACCESS <- reactive({
        if (input$mcselecttype=="cell_type") {
          sfile <- conf[[config_id]]$sum_file_genes
        } else if (input$mcselecttype=="broad_cell_type") {
          sfile <- conf[[config_id]]$sum_file_genes_broad
        }
        mat <- readRDS(file.path(INPUT_DIR,sfile))
        colnames(mat)[colnames(mat) == "neuronal_gastrula"] <- "neuronal"
        mat
      })
      
      # calculate genes summary
      gsmcs <- eventReactive(
        eventExpr = input$dosummary,
        valueExpr = {
          ns <- session$ns
          # selected cells
          cells <- rownames(sc_METADATA[input$mcselecttype %in% input$ids_mcs])
          message(sprintf("Summarizing %s cells; %s method; fc>%s", length(cells), input$fc_sum_method, input$fc_selection))
          # acessibility values
          atac_gene_summary(
            mc_ids=input$ids_mcs, method=input$fc_sum_method, fc=input$fc_selection,
            usefcbg=input$fcbg, fcbg=input$fcbg_selection, methodbg=input$fcbg_sum_method,
            mc_fp=FP(), mc_counts=ACCESS(), annt=GENE_ANNT, tfannt=TF_ANNT
          )
        }
      )
      
      # table summarizing gap genes in queried metacells
      output$genes_summary_table <- DT::renderDataTable(
        datatable(gsmcs()$gene_summary) %>% formatStyle('tf', target='row', backgroundColor = styleEqual(c("yes","no"), c("AntiqueWhite","white"))),
        rownames = FALSE,
        extensions = 'FixedColumns',
        options = list(dom = 'tp', scrollX = TRUE, fixedColumns = TRUE, ordering = TRUE, pageLength = 20),
        selection = 'none'
      )
      
      # table summarizing query and results
      output$summary  <- renderTable(
        gsmcs()$summary
      )
      
      # download multiple metacell query table
      output$download_table <- downloadHandler(
        filename <- function(){
          mc_ids_names <- red_mc_vector(selected_mcs(),range_sep="-")
          fcn <- paste0("fc",input$fc_selection)
          if (input$fcbg==TRUE)
            fcn <- paste0(fcn,"_fcbg",input$fcbg_selection)
          sprintf("mc_summary_%s_%s.tsv",mc_ids_names,fcn)
        },
        content <- function(file){
          write.table(
            gsmcs()$gene_summary,
            file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t"
          )
        }
      )
    }
  )  
  
}

peaksUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box(
        title="1. Annotation selection", width = 4, solidHeader = TRUE,
        h5("Find peaks that are acessible in all selected cell types (foreground), optionally excluding peaks accessible in other cell types (background)."),
        h5("If nothing is selected as background, all cell types not in foreground will be used."),
        h5("If you don't want to use background filtering, set accessibility in background to a very high value (e.g. 1e5)."),
        br(),
        selectInput(
          inputId=ns("ids_mcs"), 
          label="Select one or more cell types (foreground)", 
          choices=NULL, multiple=TRUE, selected=list(1),
          selectize = TRUE
        ),
        selectInput(
          inputId=ns("bgs_mcs"), 
          label="Select one or more cell types to compare against (background)",,
          choices=NULL, multiple=TRUE, selected=list(1),
          selectize = TRUE
        ),
        shiny::br()
      ),
      box(
        title="2. Accessibility threshold selection", width = 4, solidHeader = TRUE,
        numericInput(
          inputId = ns("access_selection"), label = "Choose minimum peak accessibility:", 
          min = -1e5, max = 1e5, value = 1
        ),
        numericInput(
          inputId = ns("access_bg_selection"), label = "Choose maximum peak accessibility in background:", 
          min = -1e5, max = 1e5, value = 100
        ),
        shiny::br(),
        plotOutput(ns("access_dist"), width = "100%", height = 300),
        checkboxInput(
          inputId = ns("log_scale_density"), label = "Show accessibility distribution on log10 scale", value = TRUE
        ),
        checkboxInput(
          inputId = ns("exclude_cp"), label = "Exclude CP genes", value = TRUE
        ),
      ),
      box(
        title="3. Generate summary", width = 4, solidHeader=TRUE,
        h5("Click to generate summary using specified parameters."),
        actionButton(ns("dosummary"), "Go!"),
        h4("Your search: "),
        tableOutput(ns("summary")),
        shiny::br()
      ),
      box(
        title="Master peak table", width=12, height=700, solidHeader=TRUE,
        DT::dataTableOutput(ns("peaks_summary_table"), height = 560),
        downloadButton(ns('download_peaks_table'), 'Download Peaks Table'),
        shiny::br()
      ),
      box(
        title = "Peaks overlap", width=12, height=1500, solidHeader=TRUE,
        plotOutput(ns("peak_ATAC_euler"), width = "50%"),
        shiny::br(),
        tableOutput(ns("peak_ATAC_fit")),
        shiny::br()
      )
    )
    
  )
  
}

peaksServer <- function(id, config_file="config.yaml", config_id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir
      
      # scATAC peaks
      peaks_file <- reactive({
        p_fn <- conf[['default']]$peaks_file
        peaks_file <- readRDS(file.path(INPUT_DIR, p_fn))
        peaks_file[, accessibility  := score_avg_qnorm]
        unique(peaks_file[, .(cell_type, peak, promoter, promoter_adult, promoter_gastrula, accessibility, gene, expression_fc)])
      })
      
      output$access_dist <- renderPlot({
        req(peaks_file())
        peaks_dt <- unique(peaks_file()[, .(peak,promoter,cell_type,accessibility)])
        peaks_dt <- peaks_dt[cell_type %in% c(input$ids_mcs, input$bgs_mcs)]
        if (input$exclude_cp == TRUE)
          peaks_dt <- peaks_dt[promoter != "CP"]
        gp_dens <- ggplot(peaks_dt, aes(x=accessibility)) + 
          geom_density(fill="grey") + 
          labs(subtitle = "Accessibility of peaks in selected cell types") +
          theme_minimal() +
          theme(plot.subtitle = element_text(face = "bold"))
        if (input$log_scale_density == TRUE)
          gp_dens <- gp_dens + scale_x_log10()
        gp_dens + 
          geom_vline(xintercept = input$access_bg_selection, linetype="dashed", color="red") + 
          geom_vline(xintercept = input$access_selection, linetype="dashed", color="blue")
      })
      
      # update id_mcs with values from cell_type column of peaks_file
      observe({
        pf <- peaks_file()
        updateSelectInput(
          session, "ids_mcs",
          choices = unique(pf$cell_type)
        )
        updateSelectInput(
          session, "bgs_mcs",
          choices = unique(pf$cell_type)
        )
      })
      
      # subset scATAC peaks table 
      peak_ATAC_summary <- eventReactive(
        eventExpr = input$dosummary,
        valueExpr = {
          # selected cells
          pt_mc <- peaks_file()
          
          # blacklist CP
          if (input$exclude_cp == TRUE) {
            peaks_blacklist <- pt_mc[promoter == "CP" | promoter_adult == "CP" | promoter_gastrula == "CP"]$peak
          } else {
            peaks_blacklist <- NULL
          }
          
          # find overlaps
          pt_mc_out <- feats_ct_ovl(
            ct_fg = input$ids_mcs, ct_bg = input$bgs_mcs, 
            unique(pt_mc), feat_column = "peak", 
            val_fg_column = "accessibility", val_bg_column = "accessibility",
            min_fg = input$access_selection, max_bg = input$access_bg_selection,
            level_dt = NULL, feats_blacklist = peaks_blacklist
          )
          
          # small summary table
          summary_dt <- unique(pt_mc_out[,.(peak,n_selected_cell_types)])[,.N,n_selected_cell_types]
          setnames(summary_dt, c("n_selected_cell_types","N"), c("selected_cell_types", "peaks"))
          
          # add gene info
          gene_ann <- copy(GENE_ANNT)
          gene_ann[, search_id := NULL]
          setnames(gene_ann, c("gene_id","PFAM domain", "gene name"), c("gene","pfam","name"), skip_absent = TRUE)
          gene_ann <- gene_ann[, .SD[1], gene]
          pt_mc_out <- merge.data.table(
            pt_mc_out, gene_ann, by = "gene", all.x = TRUE, sort = FALSE
          )
          pt_mc_out[is.na(name), name := ""]
          pt_mc_out[is.na(pfam), pfam := ""]
          if ("promoter_adult" %in% colnames(pt_mc_out)) {
            pt_mc_out[, promoter_adult := NULL]
          }
          if ("promoter_gastrula" %in% colnames(pt_mc_out)) {
            pt_mc_out[, promoter_gastrula := NULL]
          }
          pt_mc_out[, TF := "no"]
          pt_mc_out[gene %in% TF_ANNT[[1]], TF := "yes"]
          setcolorder(pt_mc_out, c("selected_cell_types", "cell_type", "peak", "accessibility", "gene", "pfam", "name", "TF", "expression_fc"))
          
          # eulerr fit
          eul_dt <- unique(pt_mc_out[,.(peak, cell_type)])[, val := TRUE]
          eul_dt <- dcast.data.table(
            eul_dt, peak ~ cell_type, value.var = "val", fill = FALSE
          )
          eul_dt[, c("peak") := NULL]
          class(eul_dt) <- "data.frame"
          
          # return
          list("peaks_table" = pt_mc_out, "euler_fit" = eul_dt, "summary" = summary_dt)
        }
      )
      
      # table summarizing query and results
      output$summary  <- renderTable(
        peak_ATAC_summary()$summary
      )
      
      # output peaks master table
      output$peaks_summary_table <- DT::renderDataTable(
        datatable(
          peak_ATAC_summary()$peaks_table, 
          fillContainer = TRUE, 
          selection = "none", 
          filter = list(position = 'top', clear = FALSE),
          rownames = FALSE,
          #extensions = 'Buttons', options = list(
          #  dom = 'Bfrtip'
            #buttons = 
            #  list(list(
            #    extend = 'collection',
            #    buttons = c('csv', 'excel'),
            #    text = 'Download'
            #  ))
          #)
        ) %>% formatRound(columns = c(
          "accessibility","expression_fc"
        ))
      )

      output$download_peaks_table <- downloadHandler(
        filename = function() {
          paste("peaks_summary", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(peak_ATAC_summary()$peaks_table, file, row.names = FALSE)
        }
      )

      # output venn diagram
      output$peak_ATAC_euler <- renderPlot({
        eul_dt <- peak_ATAC_summary()$euler_fit
        if (!is.null(eul_dt) && ncol(eul_dt) > 1) {
          fit <- euler(eul_dt)
          plot(
            fit,
            quantities = TRUE,
            fill = ct_cols[colnames(eul_dt)]
          )
        } else {
          print("Not enough data to plot euler diagram")
          plot_spacer()
        }
      })
      
      output$peak_ATAC_fit <- renderTable({
        eul_dt <- peak_ATAC_summary()$euler_fit
        fit <- euler(eul_dt)
        print(fit)
        fit_dt <- data.table(
          combination = names(fit$original.values),
          original = fit$original.values
        )
        fit_dt[, fitted := fit$fitted.values[combination]]
        fit_dt[, residuals := fit$residuals[combination]]
        fit_dt
      })
      
    }
    
  )
}

tfUI <- function(id) {

  ns <- NS(id)
  
  tagList(

    fluidRow(
      box(
        # TF activity
        title="TF activity", width=12, height=800, solidHeader = TRUE, 
        selectInput(
          inputId=ns("searchid"), 
          label="Search gene by name or id", 
          choices=NULL,
          selectize = TRUE
        ),
        tryCatch(girafeOutput(ns("plot_gene_act"), height=600), error=function(e) warning(e))
      ),
    ),
    fluidRow(
      box(
        title="1. Annotation selection", width = 4, solidHeader = TRUE,
        h5("Find TF genes that are active in all selected cell types (foreground), optionally excluding genes active in other cell types (background)."),
        h5("If nothing is selected as background, all cell types not in foreground will be used."),
        h5("If you don't want to use background filtering, set threshold value in background to a very high value (e.g. 1e5)."),
        br(),
        selectInput(
          inputId=ns("ids_mcs"), 
          label="Select one or more cell types (foreground)", 
          choices=NULL, multiple=TRUE, selected=list(1),
          selectize = TRUE
        ),
        selectInput(
          inputId=ns("bgs_mcs"), 
          label="Select one or more cell types to compare against (background)",,
          choices=NULL, multiple=TRUE, selected=list(1),
          selectize = TRUE
        ),
        shiny::br()
      ),
      box(
        title="2. Activity threshold selection", width = 4, solidHeader = TRUE,
        numericInput(
          inputId = ns("act_selection"), label = "Choose minimum value in selected foreground cell types:", 
          min = -1e5, max = 1e5, value = 1
        ),
        numericInput(
          inputId = ns("act_bg_selection"), label = "Choose maximum value in background:", 
          min = -1e5, max = 1e5, value = 100
        ),
        shiny::br(),
        plotOutput(ns("act_dist"), width = "100%", height = 300),
        checkboxInput(
          inputId = ns("log_scale_density"), label = "Show distribution on log10 scale", value = TRUE
        )
      ),
      box(
        title="3. Generate summary", width = 4, solidHeader=TRUE,
        h5("Click to generate summary using specified parameters."),
        actionButton(ns("dosummary"), "Go!"),
        h4("Your search: "),
        tableOutput(ns("summary")),
        shiny::br()
      ),
      box(
        title="Master TF table", width=12, height=700, solidHeader=TRUE,
        DT::dataTableOutput(ns("tfs_summary_table"), height = 560),
        downloadButton(ns('download_tfs_table'), 'Download TFs Table'),
        shiny::br()
      ),
      box(
        title = "TFs overlap", width=12, height=1500, solidHeader=TRUE,
        plotOutput(ns("tfs_act_euler"), width = "50%"),
        shiny::br(),
        tableOutput(ns("tfs_act_fit")),
        shiny::br()
      )
    )
  )
}

tfServer <- function(id, config_file="config.yaml", config_id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir
      
      # TF chrom var z score (activity) and expression
      chr_exp_fn <- file.path(INPUT_DIR, conf[['default']]$act_exp_file)
      chr_exp_ct <- readRDS(chr_exp_fn)

      # awful fix
      chr_exp_ct[cell_type == "neuronal_gastrula", cell_type := "neuronal"]
      # chr_exp_ct[, cell_type := str_replace_all(cell_type, cell_type_rename)]

      # gene search - subset only TFs with accessibility
      genes <- chr_exp_ct[['gene']]
      genes <- genes[genes %in% TF_ANNT[[1]]]
      GENE_ANNT_ATAC <- GENE_ANNT[gene_id %in% genes]
      
      updateSelectizeInput(
        session, "searchid",
        choices = GENE_ANNT_ATAC$search_id,
        selected = grep("Pou4 ", GENE_ANNT_ATAC$search_id, value = TRUE)[1],
        server = TRUE
      )
      
      # TF activity vs expression scatterplot
      plot_gene_act_f <- reactive({
        selected_gene <- GENE_ANNT_ATAC[search_id==input$searchid,gene_id]
        p_mc <- tryCatch(
          plot_tf_act_exp(chr_exp_ct, gn = selected_gene, ct_cols = ct_cols),
          error = function(e) warning(e)
        )
        girafe(
          code = print(p_mc),
          width_svg = 8, height_svg = 8,
          options = list(
            opts_selection(type = "none", only_shiny = FALSE),
            opts_hover_inv(css = "opacity:0.1;"),
            opts_hover(css = "opacity:1;")
          )
        )
      })

      output$plot_gene_act <- renderGirafe({
        shiny::req(input$searchid)
        plot_gene_act_f()
      })


      # TF actvity
      tfs_file <- reactive({
        dt <- unique(chr_exp_ct[, .(cell_type, stage, gene, gene_name, common_name, og, pfam, expression, motif_deviation)])
        cs <- unique(sort(dt$cell_type))
        if (!all(cs %in% cell_types)) {
          warning("Unknown cell types in TF activity table: ", paste(cs[!cs %in% cell_types]))
        }
        dt[, cell_type := factor(cell_type, levels = cell_types)]
        setorder(dt, cell_type, gene)
        dt
      })
      
      output$act_dist <- renderPlot({
        req(tfs_file())
        peaks_dt <- tfs_file()
        peaks_dt <- peaks_dt[cell_type %in% c(input$ids_mcs, input$bgs_mcs)]
        gp_dens <- ggplot(peaks_dt, aes(x = motif_deviation)) + 
          geom_density(fill = "grey") + 
          labs(subtitle = "Activity of TFs in selected cell types") +
          theme_minimal() +
          theme(plot.subtitle = element_text(face = "bold"))
        if (input$log_scale_density == TRUE)
          gp_dens <- gp_dens + scale_x_log10()
        gp_dens + 
          geom_vline(xintercept = input$act_bg_selection, linetype = "dashed", color = "red") + 
          geom_vline(xintercept = input$act_selection, linetype = "dashed", color  ="blue")
      })
      
      # update id_mcs with values from cell_type column of tfs_file
      observe({
        pf <- tfs_file()
        updateSelectInput(
          session, "ids_mcs",
          choices = unique(pf$cell_type)
        )
        updateSelectInput(
          session, "bgs_mcs",
          choices = unique(pf$cell_type)
        )
      })
      
      # subset scATAC peaks table 
      tfs_act_summary <- eventReactive(
        eventExpr = input$dosummary,
        valueExpr = {
          # selected cells
          pt_mc <- tfs_file()
          
          # find overlaps
          pt_mc_out <- feats_ct_ovl(
            ct_fg = input$ids_mcs, ct_bg = input$bgs_mcs, 
            unique(pt_mc), feat_column = "gene", 
            val_fg_column = "motif_deviation", val_bg_column = "motif_deviation",
            min_fg = input$act_selection, max_bg = input$act_bg_selection,
            level_dt = NULL, feats_blacklist = NULL
          )
          
          # small summary table
          summary_dt <- unique(pt_mc_out[,.(gene,n_selected_cell_types)])[,.N,n_selected_cell_types]
          setnames(summary_dt, c("n_selected_cell_types","N"), c("selected_cell_types", "TFs"))
          
          # eulerr fit
          eul_dt <- unique(pt_mc_out[,.(gene, cell_type)])[, val := TRUE]
          eul_dt <- dcast.data.table(
            eul_dt, gene ~ cell_type, value.var = "val", fill = FALSE
          )
          eul_dt[, c("gene") := NULL]
          class(eul_dt) <- "data.frame"
          
          # return
          list("tfs_table" = pt_mc_out, "euler_fit" = eul_dt, "summary" = summary_dt)
        }
      )
      
      # table summarizing query and results
      output$summary <- renderTable(
        tfs_act_summary()$summary
      )
      
      # output TFs master table
      output$tfs_summary_table <- DT::renderDataTable(
        datatable(
          tfs_act_summary()$tfs_table, 
          fillContainer = TRUE, 
          selection = "none", 
          filter = list(position = 'top', clear = FALSE),
          rownames = FALSE,
          #extensions = 'Buttons', options = list(
          #  dom = 'Bfrtip'
            #buttons = 
            #  list(list(
            #    extend = 'collection',
            #    buttons = c('csv', 'excel'),
            #    text = 'Download'
            #  ))
          #)
        ) %>% formatRound(columns = c(
          "motif_deviation","expression"
        ))
      )

      output$download_tfs_table <- downloadHandler(
        filename = function() {
          paste("tfs_summary", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(peak_act_summary()$peaks_table, file, row.names = FALSE)
        }
      )

      # output venn diagram
      output$tfs_act_euler <- renderPlot({
        eul_dt <- tfs_act_summary()$euler_fit
        if (!is.null(eul_dt) && ncol(eul_dt) > 1) {
          fit <- euler(eul_dt)
          plot(
            fit,
            quantities = TRUE,
            fill = ct_cols[colnames(eul_dt)]
          )
        } else {
          print("Not enough data to plot euler diagram")
          plot_spacer()
        }
      })
      
      output$tfs_act_fit <- renderTable({
        eul_dt <- tfs_act_summary()$euler_fit
        fit <- euler(eul_dt)
        print(fit)
        fit_dt <- data.table(
          combination = names(fit$original.values),
          original = fit$original.values
        )
        fit_dt[, fitted := fit$fitted.values[combination]]
        fit_dt[, residuals := fit$residuals[combination]]
        fit_dt
      })
      

    }
  )
}

cellTreeUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      # genes scores
      box(
        title = "Gene accessibility", width = 4, height=1500, solidHeader = TRUE,
        selectInput(
          inputId = ns("subset_cts_genes"), 
          label = "Cell types subset", 
          choices = list("All" = "scatac", "Adult" = "adult", "Gastrula" = "gastrula"),
          selected = "All",
          selectize = TRUE
        ),
        sliderInput(
          inputId = ns("genes_thr"),
          label = "Gene score fold change threshold",
          min = 1, max = 5, step = 0.1, value = 2
        ),
        plotOutput(ns("genes_tree"))
      ),
      # peaks overlaps
      box(
        title =  "Peaks overlap", width = 4, height=1500, solidHeader = TRUE,
        selectInput(
          inputId = ns("subset_cts_peaks"), 
          label = "Cell types subset", 
          choices = list("All" = "scatac", "Adult" = "adult", "Gastrula" = "gastrula"),
          selected = "All",
          selectize = TRUE
        ),
        sliderInput(
          inputId = ns("peaks_ovl_thr"),
          label = "Peak accessibility fold change threshold",
          min = 0, max = 5, step = 0.1, value = 2
        ),
        plotOutput(ns("peaks_tree"))
      ),
      # sequence features
      box(
        title = "Squence features", width = 4, height=1500, solidHeader = TRUE,
        selectInput(
          inputId = ns("subset_cts_seqs"), 
          label = "Cell types subset", 
          choices = list("All" = "scatac", "Adult" = "adult", "Gastrula" = "gastrula"),
          selected = "All",
          selectize = TRUE
        ),
        radioButtons(
          inputId = ns("model"),
          label = "Select sequence model",
          choices = list("gkm-SVM" = "gkmSVM"),
          inline = TRUE, 
          selected = "gkmSVM"
        ),
        plotOutput(ns("gksvm_tree"))
      )
    )
  )
}

cellTreeServer <- function(id, config_file="config.yaml", config_id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir
      
      # # # # # # # # # # # # 
      #         GENES       #
      # # # # # # # # # # # # 
      
      # data
      genes_mt <- reactive({
        genes_fn <- conf[[input$subset_cts_genes]]$footprint_file_genes
        genes_fp <- readRDS(file.path(INPUT_DIR, genes_fn))
        dmatrix <- parse_fps(genes_fp, fc_thr = input$genes_thr)
        as.matrix(dmatrix)
      })
      
      # tree
      output$genes_tree <- renderPlot({
        trees_funct(genes_mt(), groups = cell_types, colors = ct_cols, rename_pattern = cell_type_rename)
      }, height = 600)

      # # # # # # # # # # # # 
      #         PEAKS       #
      # # # # # # # # # # # # 
      
      # data
      peaks_ovl_mt <- reactive({
        peaks_fn <- conf[[input$subset_cts_peaks]]$footprint_file_peaks
        peaks_fp <- readRDS(file.path(INPUT_DIR, peaks_fn))
        mt <- peaks_ovl(peaks_fp, input$peaks_ovl_thr)
        as.matrix(dist(mt))
      })
      
      # tree
      output$peaks_tree <- renderPlot({
        trees_funct(peaks_ovl_mt(), groups = cell_types, colors = ct_cols, rename_pattern = cell_type_rename)
      }, height = 600)
     
      # # # # # # # # # # # # 
      #       sequence      #
      # # # # # # # # # # # # 
      
      # data
      gksvm_mt <- reactive({
        req(input$model)
        auc_fn <- switch(
          input$model,
          "gkmSVM" = conf[['default']]$auc_file,
          NULL
        )
        auc_mt <- readRDS(file.path(INPUT_DIR, auc_fn))
        if (input$subset_cts_seqs != "scatac") {
          auc_mt <- auc_mt[grep(input$subset_cts_seqs, rownames(auc_mt)), grep(input$subset_cts_seqs, colnames(auc_mt))]
        }
        as.matrix(dist(auc_mt))
      })
      
      # tree
      output$gksvm_tree <-  renderPlot({
        trees_funct(gksvm_mt(), groups = cell_types, colors = ct_cols, rename_pattern = cell_type_rename)
      }, height = 600)
      
    }
  )
}

UI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box()
    )
  )
  
}

Server <- function(id, config_file="config.yaml", config_id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
    }
    
  )
}
