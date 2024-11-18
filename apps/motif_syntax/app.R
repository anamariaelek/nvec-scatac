library(shiny)
library(data.table)
library(ggplot2)
library(ggseqlogo)
library(plotly)
library(universalmotif)
library(htmlwidgets)
library(knitr)
library(kableExtra)


bcts <- c(
    "cnidocyte", "cnidocyte_gastrula", "digestive_filaments",
    "ecto", "EMS", "epidermis", "gastro_circular_muscle",
    "gastro_parietal_muscle", "gastro", "gland_mucin", "gland",
    "muscle", "neuronal", "neuron_GATA_Islet", "neuron_Pou4_FoxL2",
    "NPC", "precursors"
)

bct_cols <- c(
    "cnidocyte"                 = "#ff42ff",
    "cnidocyte_gastrula"        = "#f7abf7",
    "ecto"                      = "#51a0be",
    "EMS"                       = "#bdf5bd",
    "gastro_circular_muscle"    = "#73b009",
    "gastro_parietal_muscle"    = "#8ceb10",
    "gastro"                    = "#85c90e",
    "muscle"                    = "#ffd700",
    "digestive_filaments"       = "#e33d3d",
    "precursors"                = "#bebebe",
    "NPC"                       = "#808d91",
    "epidermis"                 = "#04ccd4",
    "neuron_GATA_Islet"         = "#1175f0",
    "neuron_Pou4_FoxL2"         = "#101cde",
    "neuronal"                  = "#063cb9",
    "gland"                     = "#ff6f08",
    "gland_mucin"               = "#ff8f12"
)
bcts <- names(bct_cols)

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

ui <- fluidPage(
  titlePanel("Motif Co-occurrence Enrichment"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons(
        "lvl",
        "Select peaks for which to show enrichmnet:",
        choices = c(
          "Cell type peaks" = "cell_type", 
          "Differential cell type peaks" = "cell_type_differential"
        ),
        selected = "cell_type"
      ),
      HTML("<b>Suggested parameters:</b><br><br>"),
      HTML("<em>Cell type peaks</em><br>"),
      HTML("Number of peaks with motif pair: 100<br>"),
      HTML("Fraction of peaks with motif pair: 0.02<br>"),
      HTML("Pair co-occurence enrichment: 1<br><br>"),
      HTML("<em>Differential cell type peaks</em><br>"),
      HTML("Number of peaks with motif pair: 10<br>"),
      HTML("Fraction of peaks with motif pair: 0.02<br>"),
      HTML("Pair co-occurence enrichment: 1.5<br><br>"),
      HTML("<b>Number of peaks:</b><br><br>"),
      uiOutput("peaks_info")
    ),
    mainPanel = mainPanel(
      tabsetPanel(
        tabPanel(
          "Heatmap",
          br(),
          fluidRow(
            column(2,
              numericInput(
                "pair_count",
                "Number of peaks with motif pair",
                min = 0, max = 1e10, value = 100, step = 1)
            ),
            column(2,
              numericInput(
                "pair_frac",
                "Fraction of peaks with motif pair",
                min = 0, max = 1, value = 0.02, step = 0.01)
            ), 
            column(2,
              numericInput(
                "pair_enr",
                "Pair co-occurence enrichment",
                min = 0, max = 10, value = 1, step = 0.01)
            )
          ),
          fluidRow(
            column(2,
              actionButton("reload_button_heatmap", "Load heatmap"),
              br(), br()
            )
          ),
          fluidRow(
            column(6,
              h4("Motif pair co-occurrence heatmap:"),
              HTML("Click on a heatmap cell to see motif pair details."),
              plotlyOutput("enrich_heatmap", height = "1000px")
            ),
            column(6,
              h4("Selected motif pair:"),
              plotOutput("heatmap_barplot", height = "400px"),
              br(),
              plotOutput("motif1_hm_logo", height = "100px"),
              uiOutput("motif1_hm_info"),
              plotOutput("motif2_hm_logo", height = "100px"),
              uiOutput("motif2_hm_info"),
              h4("Cell type co-occurence enrichment:"),
              uiOutput("heatmap_row_table")
            )
          )
        ),
        tabPanel(
          "Cell type motif pairs",
          br(),
          fluidRow(
            column(7,
              selectInput("cell_type", "Select cell type:", choices = cts),
              actionButton("reload_button", "Load data"),
              br(), br(),
              HTML("Click on a point to see motif pair details."),
              plotlyOutput("plot", height = "800px"),
              uiOutput("motif_pair_info")
            ),
            column(5,
              br(),
              plotOutput("motif1_logo", height = "100px"),
              uiOutput("motif1_info"),
              plotOutput("motif2_logo", height = "100px"),
              uiOutput("motif2_info"),
              br(),
              plotOutput("enr_bars", height = "500px")
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
  
  # Load peaks
  pks_dt <- reactive({
    fread(file.path(dat_dir, input$lvl, sprintf("Peaks_per_%s_mapped.tsv.gz", input$lvl)))
  })
  # Load enrichment data across cell types
  mta_enr_cts <- reactive({
    mta_enr_cts <- readRDS(
      file.path(dat_dir, input$lvl, "motif-co-occurrences-enrichment.rds")
    )
    if (input$lvl %in% c("cell_type", "cell_type_differential")) {
      mta_enr_cts[, cell_type := factor(cell_type, levels = cts)]
    } else if (input$lvl == "broad_cell_type") {
      mta_enr_cts[, cell_type := factor(cell_type, levels = bcts)]
    }
  })

  # Table with peaks info
  output$peaks_info <- renderUI({
    pks_dt <- unique(pks_dt()[, .(peak, cell_type)])
    pks_dt <- pks_dt[, .("number of peaks" = .N), cell_type]
    setorder(pks_dt, -`number of peaks`)
    if (nrow(pks_dt) == 0) return(NULL)
    pks_dt <- kable(pks_dt, format = "html", col.names = NULL) %>%
      kable_styling("striped", full_width = TRUE)
    HTML(pks_dt)
  })

  # Load motif PWM data
  mta_pwm <- read_meme(file.path(dat_dir, "motif-archetypes-all.meme"))
  names(mta_pwm) <- sapply(mta_pwm, function(m) m@name)

  # # # # # # # # # # #
  #     Enrichment    #
  # # # # # # # # # # #

  # Reactive value to store data
  rv <- reactiveValues(
    mta_enr_gp = NULL,
    mta_pwm = NULL,
    mta_enr_gen = NULL
  )
  
  # Reload data when button is clicked
  observeEvent(input$reload_button, {
    req(input$cell_type)

    mta_enr_gen <- mta_enr_cts()[cell_type == input$cell_type]

    # Filter and store data
    rv$mta_enr_gen <- mta_enr_gen
    rv$mta_enr_gp <- mta_enr_gen[frac_peak_pair > 0.01]
  })
  
  # Plot limits helper
  plot_limits_center <- function(x, center = 0) {
    xmax <- max(abs(x) - center)
    c(-xmax, xmax)
  }
  
  # Render Plotly plot
  output$plot <- renderPlotly({
    req(rv$mta_enr_gp)
    gp_enr <- ggplot(rv$mta_enr_gp, aes(count_peak_pair, log_enrichment)) +
      geom_point(aes(color = motifs_similarity, text = gene_pair)) +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(limits = plot_limits_center) +
      scale_color_viridis_c() +
      labs(
        x = "Co-occurrence peaks",
        y = "Co-occurrence enrichment (log2)",
        color = "Motifs\nsimilarity"
      ) +
      theme_minimal()
    ggplotly(gp_enr, tooltip = "text")
  })
  
  # Render motif logos
  output$motif1_logo <- renderPlot({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    ggseqlogo(mta_pwm[[selected$motif1]]@motif) + theme_void()
  })
  
  output$motif2_logo <- renderPlot({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    ggseqlogo(mta_pwm[[selected$motif2]]@motif) + theme_void()
  })
  
  # Render detailed tables
  output$motif1_info <- renderUI({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp, rv$mta_enr_gen)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    motif1_info <- rv$mta_enr_gen[
      motif1 == selected$motif1 & motif2 == selected$motif2, .(
        motif1, count_peak_motif1, frac_peak_motif1,
        gene_motif1, gene_name_motif1, tf_family_motif1
      )
    ]
    if (nrow(motif1_info) == 0) return(NULL)
    motif1_table <- kable(t(motif1_info), format = "html", col.names = NULL) %>%
      kable_styling("striped", full_width = TRUE)
    HTML(motif1_table)
  })

  output$motif2_info <- renderUI({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp, rv$mta_enr_gen)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    motif2_info <- rv$mta_enr_gen[
      motif1 == selected$motif1 & motif2 == selected$motif2, .(
        motif2, count_peak_motif2, frac_peak_motif2,
        gene_motif2, gene_name_motif2, tf_family_motif2
      )
    ]
    if (nrow(motif2_info) == 0) return(NULL)
    motif2_table <- kable(t(motif2_info), format = "html", col.names = NULL) %>%
      kable_styling("striped", full_width = TRUE)
    HTML(motif2_table)
  })

  output$motif_pair_info <- renderUI({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp, rv$mta_enr_gen)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    pair_info <- rv$mta_enr_gen[
      motif1 == selected$motif1 & motif2 == selected$motif2, .(
        motif1, motif2, gene_pair, count_peak_pair, frac_peak_pair, log_enrichment, motifs_similarity
      )
    ]
    if (nrow(pair_info) == 0) return(NULL)
    pair_table <- kable(t(pair_info), format = "html") %>%
      kable_styling("striped", full_width = TRUE)
    tagList(
      h4("Selected Motif Pair"),
      HTML(pair_table)
    )
  })
  
  # Render enrichment accross cell types
  output$enr_bars <- renderPlot({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp, rv$mta_enr_gen)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    sel_mta_enr_cts <- mta_enr_cts()[motif1 == selected$motif1 & motif2 == selected$motif2]
    ggplot(sel_mta_enr_cts, aes(cell_type, log_enrichment, fill = cell_type)) + 
      geom_bar(stat = "identity", color = "black") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = c(ct_cols, bct_cols)) +
      labs(
        x = "Cell type",
        y = "Enrichment (log2)",
        title = "Motif pair enrichment across cell types"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none"
      )
  })

  # # # # ## # # # #
  #     Heatmap    #
  # # # # ## # # # #

  # Filter motifs pairs to plot as heatmap
  rv_hm <- reactiveValues(
    enr_mts = NULL,
    mta_enr_gp = NULL
  )
  observeEvent(input$reload_button_heatmap, {
    
    # Select motif pairs
    enr_mts <- unique(
      mta_enr_cts()[
        frac_peak_pair > input$pair_frac &
        enrichment > input$pair_enr &
        count_peak_pair > input$pair_count
      ]$motif_pair
    )
    rv_hm$enr_mts <- enr_mts

    # Subset data
    mta_enr_gp <- mta_enr_cts()[motif_pair %in% enr_mts]
        
    # Order pairs
    stopifnot(all(mta_enr_gp$cell_type %in% cts))
    mta_enr_gp[, cell_type := factor(cell_type, levels = cts)]
    mta_enr_gp[order(motif_pair, -enrichment), max_enr_ct := as.integer(.SD[1]$cell_type), motif_pair]
    setorder(mta_enr_gp, max_enr_ct, -enrichment, motif_pair)
    mta_enr_gp[, motif_pair := factor(motif_pair, levels = unique(mta_enr_gp$motif_pair))]
    mta_enr_gp[, gene_pair := factor(gene_pair, levels = unique(mta_enr_gp$gene_pair))]

    # Ranges for plotting
    mta_enr_gp[, log_enrichment_trimmed := log_enrichment]
    mta_enr_gp[, log_enrichment_trimmed := pmax(log_enrichment_trimmed, 0)]
    mta_enr_gp[, log_enrichment_trimmed := pmin(log_enrichment_trimmed, 2)]

    rv_hm$mta_enr_gp <- mta_enr_gp
    
  })


  # Plot heatmap
  col_vec <- c("#4d4d4d", "#e0e0e0", "#b2182b")
  output$enrich_heatmap <- renderPlot({
    req(!is.null(rv_hm$mta_enr_gp))
    ggplot(rv_hm$mta_enr_gp, aes(cell_type, gene_pair, fill = log_enrichment_trimmed)) +
      geom_tile() +
      scale_fill_gradient2(
        low = col_vec[1], mid = col_vec[2], high = col_vec[3], midpoint = 1
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      labs(
        y = sprintf("%s motif pairs", length(unique(rv_hm$mta_enr_gp$gene_pair))),
      )
  })
output$enrich_heatmap <- renderPlotly({
  req(!is.null(rv_hm$mta_enr_gp))

  # color vector and corresponding positions
  col_vec <- c("#4d4d4d", "#e0e0e0", "#b2182b")
  positions <- c(0, 0.5, 1)

  # Normalize the actual values to the range [0, 1]
  actual_values <- rv_hm$mta_enr_gp$log_enrichment_trimmed
  z_min <- min(actual_values)
  z_max <- max(actual_values)
  normalized_positions <- (actual_values - z_min) / (z_max - z_min)
  colorscale <- lapply(seq_along(normalized_positions), function(i) list(normalized_positions[i], col_vec[i]))

  # Prepare data for plotly heatmap
  heatmap_data <- rv_hm$mta_enr_gp
  heatmap_data[, cell_type := factor(cell_type, levels = cts)]
  heatmap_data[, gene_pair := factor(gene_pair, levels = unique(heatmap_data$gene_pair))]
  p <- plot_ly(
    data = heatmap_data,
    x = ~cell_type,
    y = ~gene_pair,
    z = ~log_enrichment_trimmed,
    type = "heatmap",
    colorscale = colorscale, # Map colors to actual values
    zmin = z_min, # Ensure Plotly knows the actual z-range
    zmax = z_max,
    #colorscale = list(c(0, 0.5, 1), col_vec),
    colorbar = list(title = "Log2\nco-occurence\nenrichment"),
    source = "heatmap_click"
  ) %>%
    layout(
      xaxis = list(
        title = sprintf("%s cell types", length(unique(heatmap_data$cell_type))),
        tickangle = -90
      ),
      yaxis = list(
        title = sprintf("%s motif pairs", length(unique(heatmap_data$gene_pair))),
        showticklabels = FALSE
      )
    )
    
    # Register the click event
    event_register(p, "plotly_click")
  })

  # Add reactive value for selected heatmap row
  rv_hm_selected <- reactiveValues(selected_row = NULL)

  # Capture heatmap click event
  observe({
    heatmap_event <- event_data("plotly_click", source = "heatmap_click")
    req(heatmap_event, rv_hm$mta_enr_gp)
    
    # Extract selected row (gene_pair)
    clicked_gene_pair <- rv_hm$mta_enr_gp$gene_pair[heatmap_event$pointNumber + 1]
    rv_hm_selected$selected_row <- clicked_gene_pair
  })
  
  # Render barplot for selected row
  output$heatmap_barplot <- renderPlot({
    req(rv_hm_selected$selected_row, rv_hm$mta_enr_gp)
    
    # Selected data
    selected_data <- rv_hm$mta_enr_gp[gene_pair == rv_hm_selected$selected_row]
    
    # Add 0s for missing cell types
    all_combinations <- CJ(
      gene_pair = unique(selected_data$gene_pair),
      cell_type = unique(selected_data$cell_type)
    )

    # Join with the original data.table
    selected_data <- merge(
      selected_data,
      all_combinations,
      by = c("gene_pair", "cell_type"),
      all = TRUE
    )
    selected_data[is.na(enrichment), enrichment := 0]

    # Plot
    ggplot(selected_data, aes(cell_type, enrichment, fill = cell_type)) +
      geom_bar(stat = "identity", color = "black") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = c(ct_cols, bct_cols)) +
      labs(
        x = "Cell type",
        y = "Enrichment (log2)",
        title = rv_hm_selected$selected_row
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none"
      )
  })

  # Render motif logos
  output$motif1_hm_logo <- renderPlot({
    event_data <- event_data("plotly_click", source = "heatmap_click")
    req(rv_hm_selected$selected_row, rv_hm$mta_enr_gp)
    selected <- rv_hm$mta_enr_gp[event_data$pointNumber + 1]
    ggseqlogo(mta_pwm[[selected$motif1]]@motif) + theme_void()
  })
  
  output$motif2_hm_logo <- renderPlot({
    event_data <- event_data("plotly_click", source = "heatmap_click")
    req(rv_hm_selected$selected_row, rv_hm$mta_enr_gp)
    selected <- rv_hm$mta_enr_gp[event_data$pointNumber + 1]
    ggseqlogo(mta_pwm[[selected$motif2]]@motif) + theme_void()
  })

  # Render motif tables
  output$motif1_hm_info <- renderUI({
    event_data <- event_data("plotly_click", source = "heatmap_click")
    req(event_data, rv_hm$mta_enr_gp)
    selected <- rv_hm$mta_enr_gp[event_data$pointNumber + 1]
    motif1_info <- unique(selected[, .(
      motif1, gene_motif1, gene_name_motif1, tf_family_motif1
    )])
    if (nrow(motif1_info) == 0) return(NULL)
    motif1_table <- kable(t(motif1_info), format = "html", col.names = NULL) %>%
      kable_styling("striped", full_width = TRUE)
    HTML(motif1_table)
  })

  output$motif2_hm_info <- renderUI({
    event_data <- event_data("plotly_click", source = "heatmap_click")
    req(event_data, rv_hm$mta_enr_gp)
    selected <- rv_hm$mta_enr_gp[event_data$pointNumber + 1]
    motif2_info <- unique(selected[, .(
      motif2, gene_motif2, gene_name_motif2, tf_family_motif2
    )])
    if (nrow(motif2_info) == 0) return(NULL)
    motif2_table <- kable(t(motif2_info), format = "html", col.names = NULL) %>%
      kable_styling("striped", full_width = TRUE)
    HTML(motif2_table)
  })

  # Render table for selected row
  output$heatmap_row_table <- renderUI({
    req(rv_hm_selected$selected_row, rv_hm$mta_enr_gp)
    
    selected_table <- rv_hm$mta_enr_gp[gene_pair == rv_hm_selected$selected_row, .(
      cell_type, enrichment, frac_peak_pair, count_peak_pair, motifs_similarity
    )]
    setorder(selected_table, -enrichment)
    if (nrow(selected_table) == 0) return(NULL)
    
    HTML(
      kable(selected_table, format = "html", row.names = FALSE) %>%
        kable_styling("striped", full_width = TRUE)
    )
  })


}

shinyApp(ui = ui, server = server)
