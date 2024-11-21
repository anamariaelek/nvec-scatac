library(shiny)
library(data.table)
library(ggplot2)
library(ggseqlogo)
library(plotly)
library(stringr)
library(universalmotif)
library(htmlwidgets)
library(knitr)
library(kableExtra)

# Cell types
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
cts <- names(ct_cols)

# UI
ui <- fluidPage(
  titlePanel("Motif Co-occurrence Enrichment"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        HTML("This app allows to explore motif pairs co-occurrence in cell type specific peaks."),
        br(),
        HTML("It might take a few minutes to start, because loading counts for all the motif pairs in all peaks is slow. I am terribly sorry about that."),
        br(),
        HTML("Below you can set the thresholds to filter peaks in which to summarize the motif pairs occurences."),
        br(),
        HTML("Press 'Load peaks' to apply the filters. This might take a few moments."),
        br(),
        HTML("The table with number of peaks and motifs per cell type will be displayed below."),
        br(),
        HTML("Then you can check enriched motif pairs on the right."),
        br()
      ),
      # Filter input peaks by cell type specificity
      br(),
      fluidRow(
        numericInput(
          "peak_log2fc", "Peak cell type Log2FC threshold",
          min = 0, max = 10, value = 1, step = 0.1
        ),
        HTML("To select all accessible peaks, set this to 0."),
        HTML("To select differential cell type specific peaks, set this to a value greater than 1."),
        br(), br(),
        numericInput(
          "peak_fdr", "Peak cell type FDR threshold",
          min = 0, max = 1, value = 1, step = 0.1
        ),
        HTML("To select all accessible peaks, set this to 1."),
        HTML("To select differential cell type specific peaks, set this to a lower value."),
        br(), br(),
        # Filter input motifs by cell type specificity
        numericInput(
          "motif_log2fc", "Motif cell type Log2FC threshold",
          min = 0, max = 10, value = 1, step = 0.1
        ),
        HTML("To select all motifs, set this to 0."),
        HTML("To select cell type enriched motifs, set this to a value greater than 1."),
        br(), br(),
        numericInput(
          "motif_padj", "Motif cell type FDR threshold",
          min = 0, max = 10, value = 0.05, step = 0.1
        ),
        HTML("To select all motifs, set this to 1."),
        HTML("To select cell type enriched motifs, set this to a lower value."),
        br(), br(),
        actionButton("reload_peaks", "Load peaks"),
        br(),
      ),
      fluidRow(
        uiOutput("peaks_info")
      )
    ),
    mainPanel = mainPanel(
      tabsetPanel(
        tabPanel(
          "Heatmap",
          br(),
          fluidRow(
            column(2,
              # Select cell types
              selectInput(
                "cell_type",
                "Cell type enriched motif pairs",
                choices = cts,
                selected = "cnidocyte"
              )
            ),
            column(2,
              # Filter peaks by number of motif pairs
              numericInput(
                "pair_count",
                "Number of peaks with motif pair",
                min = 0, max = 1e10, value = 100, step = 1)
            ),
            column(2,
              # Filter peaks by fraction of motif pairs
              numericInput(
                "pair_frac",
                "Fraction of peaks with motif pair",
                min = 0, max = 1, value = 0.05, step = 0.01)
            )
          ),
          fluidRow(
            column(8,
              HTML("Press 'Load heatmap' to filter and cluster motif pairs in selected cell type."),
              HTML("This might take a few moments."),
              br(), br(),
              actionButton("reload_button_heatmap", "Load heatmap"),
              br(), br()
            )
          ),
          fluidRow(
            column(6,
              h4("Motif pair co-occurrence heatmap:"),
              HTML("Click on a heatmap cell to see motif pair details."),
              plotlyOutput("enrich_heatmap", height = "800px"),
              h4("Cell type co-occurence:"),
              htmlOutput("motif_pair_info")
            ),
            column(6,
              br(),
              plotOutput("motif1_logo", height = "100px"),
              uiOutput("motif1_info"),
              plotOutput("motif2_logo", height = "100px"),
              uiOutput("motif2_info"),
              plotOutput("motif_pair_bar", height = "500px")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # Paths and constants
  q <- 0.95
  arc_id <- "PPM-PCC-0.8-IC0.5-5bp"
  dat_dir <- "data"

  # Load per cell type peaks
  pks_dt <- readRDS(file.path(
    dat_dir, "peaks-per-cell-type-annotation.rds"
  ))

  # Load motifs cell type enrichment
  mta_enr_dt <- readRDS(file.path(
    dat_dir, sprintf(
      "motif-enrichment-cell-type-%s-mona-q-%s-enr.rds",
      arc_id, q
    )
  ))

  # Load data for co-occurence of motifs assigned to TFs
  mta_pair_novl_gen <- readRDS(file.path(
    dat_dir, sprintf(
      "motif-co-occurrences-nonovl-%s-mona-q-%s-TFs.rds",
      arc_id, q
    )
  ))

  # Reactive variables
  rv_mt <- reactiveValues(
    enr_mts = NULL, # motifs enriched in cell types
    pks_sig = NULL, # peaks enriched in cell types
    cts_num = NULL, # peak and motif counts in cell types
    mtd_cnt = NULL # motif pairs counts in cell types
  )

  observeEvent(input$reload_peaks, {
    print(sprintf(
      "Filtering peaks: Log2FC >= %.2f & FDR <= %.2f",
      input$peak_log2fc, input$peak_fdr
    ))
    
    # Filter enriched peaks
    pks_sig <- pks_dt[
      Log2FC >= input$peak_log2fc &
      FDR <= input$peak_fdr
    ]
    rv_mt$pks_sig <- pks_sig

    # Count peaks
    pks_num <- pks_sig[, .(n_peaks = .N), cell_type]
    if (nrow(pks_sig) > 0) {

      # Filter enriched motifs
      mta_enr_sig <-  mta_enr_dt[
        padj < input$motif_padj &
        log2(fc) > input$motif_log2fc
      ]
      rv_mt$enr_mts <- mta_enr_sig

      # Count number of motifs
      mta_num <- mta_enr_sig[, .(n_motifs = .N), cell_type]
      
      # Per cell type counts
      cts_num <- merge.data.table(pks_num, mta_num, by = "cell_type")
      setorder(cts_num, -n_peaks)
      rv_mt$cts_num <- cts_num

      # Subset motif occurence data for enriched motifs
      mta_pair_novl_enr <- mta_pair_novl_gen[
        motif1 %in% mta_enr_sig$motif |
        motif2 %in% mta_enr_sig$motif
      ]

      # Subset motif occurence data for hits in significant peaks
      mta_pair_novl_pks <- merge.data.table(
        mta_pair_novl_enr, pks_sig, by = "peak",
        all = FALSE, allow.cartesian = TRUE
      )

      # Count motif pairs in significan peaks per cell type
      mta_pair_novl_sig <- mta_pair_novl_pks[, .(count_peak_pair = .N), .(
        cell_type,
        motif1, motif2,
        gene_motif1, gene_motif2,
        gene_name_motif1, gene_name_motif2,
        tf_family_motif1, tf_family_motif2
      )]

      # Add number of significant peaks and motifs
      mdt_cnt <- merge.data.table(
        mta_pair_novl_sig, cts_num, by = "cell_type"
      )
      mdt_cnt[, frac_peak_pair := count_peak_pair / n_peaks, by = 1:nrow(mdt_cnt)]

      # add motif pair column
      mdt_cnt[, motif_pair := sprintf(
        "%s + %s", motif1, motif2
      ), by = 1:nrow(mdt_cnt)]
      mdt_cnt[, gene_pair := sprintf(
        "%s + %s", gene_motif1, gene_motif2
      ), by = 1:nrow(mdt_cnt)]
      mdt_cnt[, gene_name_pair := sprintf(
        "%s + %s", gene_name_motif1, gene_name_motif2
      ), by = 1:nrow(mdt_cnt)]

      # factor cell types
      mdt_cnt[, cell_type := factor(cell_type, levels = cts)]
      rv_mt$mtd_cnt <- unique(mdt_cnt)

    }

  })

  # Table with peaks info
  output$peaks_info <- renderUI({
    req(!is.null(rv_mt$cts_num))
    cts_dt <- unique(rv_mt$cts_num)
    setnames(
      cts_dt,
      c("cell_type", "n_peaks", "n_motifs"),
      c("Cell type", "Number of peaks", "Number of motifs")
    )
    if (nrow(cts_dt) == 0) return(NULL)
    cts_dt <- kable(cts_dt, format = "html") %>%
      kable_styling("striped", full_width = TRUE)
    HTML(cts_dt)
  })

  # Load motif PWM data
  mta_pwm <- read_meme(file.path(dat_dir, "motif-archetypes-all.meme"))
  names(mta_pwm) <- sapply(mta_pwm, function(m) m@name)

  # # # # ## # # # #
  #     Heatmap    #
  # # # # ## # # # #

  # Filter motifs pairs to plot as heatmap
  rv_hm <- reactiveValues(
    enr_mts = NULL, # vector of selected enriched motif pairs
    mtd_cnt = NULL # table with enriched motif pairs
  )
  observeEvent(
    ignoreInit = TRUE,
    list(input$reload_peaks, input$reload_button_heatmap),
    {
      rv_hm$enr_mts <- NULL
      rv_hm$mtd_cnt <- NULL
      req(!is.null(rv_mt$mtd_cnt))

      print(sprintf(
        "Selecting motif pairs in %s: fraction >= %.2f & count <= %.2f",
        input$cell_type, input$pair_frac, input$pair_count
      ))

      # Select motif pairs
      mts <- rv_mt$enr_mts[cell_type == input$cell_type]$motif
      enr_mts <- unique(
        rv_mt$mtd_cnt[
          cell_type == input$cell_type &
          frac_peak_pair > input$pair_frac &
          count_peak_pair > input$pair_count &
          (motif1 %in% mts | motif2 %in% mts)
        ]$motif_pair
      )
      rv_hm$enr_mts <- enr_mts

      # Subset data for selected motifs
      mtd_cnt <- rv_mt$mtd_cnt[motif_pair %in% enr_mts]
          
      # Order cell types
      stopifnot(all(mtd_cnt$cell_type %in% cts))
      mtd_cnt[, cell_type := factor(
        cell_type, levels = intersect(cts, unique(mtd_cnt$cell_type))
      )]
      print(sprintf("There are %d motif pairs", length(unique(mtd_cnt$motif_pair))))

      # Cluster motif pairs
      if (nrow(mtd_cnt) > 2) {
        print(sprintf("Clustering motif pairs"))
        tryCatch({
          mdt_dc <- dcast.data.table(
            unique(mtd_cnt[, .(motif_pair, cell_type, frac_peak_pair)]),
            motif_pair ~ cell_type, value.var = "frac_peak_pair"
          )
          mdt_mt <- as.matrix(mdt_dc[, -1])
          rownames(mdt_mt) <- mdt_dc[[1]]
          hclust_res <- stats::hclust(dist(
            mdt_mt, method = "manhattan"
          ), method = "ward.D2")
          motif_pair_lvl <- rownames(mdt_mt)[hclust_res$order]
          mtd_cnt[, motif_pair := factor(motif_pair, levels = motif_pair_lvl)]
        }, error = function(e) {
          message(e)
          mtd_cnt[, motif_pair := factor(motif_pair)]
        })
        mtd_cnt[, motif_pair_idx := as.integer(motif_pair)]
        setorder(mtd_cnt, motif_pair_idx)
        #print(head(unique(mtd_cnt[, .(gene_name_pair, motif_pair_idx)])))
        #print("...")
        #print(tail(unique(mtd_cnt[, .(gene_name_pair, motif_pair_idx)])))
        rv_hm$mtd_cnt <- mtd_cnt
        print(sprintf("Selected %d motif pairs", length(unique(mtd_cnt$motif_pair))))
      }
    }
  )


  # Plot heatmap
  output$enrich_heatmap <- renderPlotly({
    req(!is.null(rv_hm$mtd_cnt))
    print(sprintf("Ploting heatmap"))
    hm_gp <- ggplot(
        rv_hm$mtd_cnt,
        aes(
          cell_type, motif_pair,
          fill = frac_peak_pair, text = gene_name_pair
        )
      ) +
      geom_tile() +
      scale_fill_gradientn(
        colours = c("#edf8b1", "#7fcdbb", "#2c7fb8", "#225ea8", "#081d58", "#000b29"), 
        na.value = "#edf8b1"
      ) +
      labs(
        x = "cell type peaks",
        y = sprintf(
          "%s motif pairs present in > %.0f%% peaks", 
          length(unique(rv_hm$mtd_cnt$motif_pair)),
          input$pair_frac * 100
        ),
        fill = "fraction of \ncell type peaks\nwith both motifs"
      ) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right"
      )
    
    # Make interactive plot
    p <- ggplotly(hm_gp, tooltip = "text")

    # Register the click event
    event_register(p, "plotly_click")
    p
  })

  
  # Render motif logos
  output$motif1_logo <- renderPlot({
    event_data <- event_data("plotly_click")
    req(event_data, rv_hm$mtd_cnt)
    selected_row <- as.integer(event_data$pointNumber[[1]][1] + 1)
    selected <- rv_hm$mtd_cnt[motif_pair_idx == selected_row][1]
    print(sprintf("selected_row: %s", selected_row))
    tryCatch({
      ggseqlogo(mta_pwm[[selected$motif1]]@motif) + theme_void()
    }, error = function(e) {
      message(e)
      ggplot() + theme_void()
    })
  })
  
  output$motif2_logo <- renderPlot({
    event_data <- event_data("plotly_click")
    req(event_data, rv_hm$mtd_cnt)
    selected_row <- as.integer(event_data$pointNumber[[1]][1] + 1)
    selected <- rv_hm$mtd_cnt[motif_pair_idx == selected_row][1]
    tryCatch({
      ggseqlogo(mta_pwm[[selected$motif2]]@motif) + theme_void()
    }, error = function(e) {
      message(e)
      ggplot() + theme_void()
    })
  })
  
  # Render detailed tables
  output$motif1_info <- renderUI({
    event_data <- event_data("plotly_click")
    req(event_data, rv_hm$mtd_cnt)
    selected_row <- as.integer(event_data$pointNumber[[1]][1] + 1)
    selected <- rv_hm$mtd_cnt[motif_pair_idx == selected_row][1]
    motif1_info <- rv_hm$mtd_cnt[
      motif1 == selected$motif1 & motif2 == selected$motif2 & cell_type == input$cell_type, .(
        motif1, gene_motif1, gene_name_motif1, tf_family_motif1
      )
    ]
    if (nrow(motif1_info) == 0) {
      return(NULL)
    } else {
      print(motif1_info)
      motif1_table <- kable(t(motif1_info), format = "html", col.names = NULL) %>%
        kable_styling("striped", full_width = TRUE)
      HTML(motif1_table)
    }
  })
  
  output$motif2_info <- renderUI({
    event_data <- event_data("plotly_click")
    req(event_data, rv_hm$mtd_cnt)
    selected_row <- as.integer(event_data$pointNumber[[1]][1] + 1)
    selected <- rv_hm$mtd_cnt[motif_pair_idx == selected_row][1]
    motif2_info <- rv_hm$mtd_cnt[
      motif1 == selected$motif1 & motif2 == selected$motif2 & cell_type == input$cell_type, .(
        motif2, gene_motif2, gene_name_motif2, tf_family_motif2
      )
    ]
    if (nrow(motif2_info) == 0) {
      return(NULL)
    } else {
      print(motif2_info)
      motif2_table <- kable(t(motif2_info), format = "html", col.names = NULL) %>%
        kable_styling("striped", full_width = TRUE)
      HTML(motif2_table)
    }
  })

  pair_info <- reactive({
    event_data <- event_data("plotly_click")
    req(event_data, rv_hm$mtd_cnt)
    selected_row <- as.integer(event_data$pointNumber[[1]][1] + 1)
    selected <- rv_hm$mtd_cnt[motif_pair_idx == selected_row][1]
    pair_info <- rv_hm$mtd_cnt[
      motif1 == selected$motif1 & motif2 == selected$motif2, .(
        cell_type, count_peak_pair, frac_peak_pair
      )
    ]
    setorder(pair_info, -frac_peak_pair)
    pair_info
  })

  
  output$motif_pair_bar <- renderPlot({
    req(event_data("plotly_click"))
    dt <- pair_info()
    dt[, cell_type := factor(
      cell_type, levels = intersect(cts, unique(dt$cell_type))
    )]
    dt <- melt.data.table(
      dt, 
      measure.vars = c("count_peak_pair", "frac_peak_pair"),
      variable.name = "count_type", value.name = "count"
    )
    ggplot(dt, aes(cell_type, count, fill = cell_type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = ct_cols) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      facet_grid(
        count_type ~ ., 
        scales = "free_y", 
        switch = "y",
        labeller = labeller(count_type = c(
          "count_peak_pair" = "Peak count",
          "frac_peak_pair" = "Peak fraction"
        ))
      ) +
      labs(y = "Count of peaks\nwith motif pair") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
      )
  })
  
  output$motif_pair_info <- renderUI({
    req(event_data("plotly_click"))
    if (nrow(pair_info()) == 0) {
      return(NULL)
    } else {
      pair_table <- kable(pair_info(), format = "html") %>%
        kable_styling("striped", full_width = TRUE)
      tagList(
        HTML(pair_table)
      )
    }
  })

}

shinyApp(ui = ui, server = server)
