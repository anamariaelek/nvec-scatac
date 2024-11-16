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

ui <- fluidPage(
  titlePanel("Motif Co-occurrence Enrichment"),
  fluidRow(
    column(8, 
      selectInput("cell_type", "Select cell type:", choices = bcts),
      actionButton("reload_button", "Load data"),
      plotlyOutput("plot", height = "800px"),
      uiOutput("motif_pair_info")
    ),
    column(4,
      br(),
      plotOutput("motif1_logo", height = "150px"),
      uiOutput("motif1_info"),
      plotOutput("motif2_logo", height = "150px"),
      uiOutput("motif2_info"),
      br(),
      plotOutput("enr_bars", height = "500px")
    )
  )
)

server <- function(input, output, session) {
  # Paths and constants
  lvl <- "broad_cell_type"
  lvl <- "cell_type"
  q <- 0.95
  arc_id <- "PPM-PCC-0.8-IC0.5-5bp"
  dat_dir <- "data"
  
  # Load enrichment data across cell types
  mta_enr_cts <- readRDS(
    file.path(dat_dir, lvl, "motif-co-occurrences-enrichment.rds")
  )
  mta_enr_cts[, cell_type := factor(cell_type, levels = bcts)]
  

  # Reactive value to store data
  rv <- reactiveValues(
    mta_enr_gp = NULL,
    mta_pwm = NULL,
    mta_enr_gen = NULL
  )
  
  # Reload data when button is clicked
  observeEvent(input$reload_button, {
    req(input$cell_type)

    # Load enrichment data for selected cell type
    mta_enr <- fread(
      file.path(dat_dir, lvl, sprintf(
        "motif-co-occurrences-enrichment-%s-%s-mona-q-%s.tsv.gz",
        input$cell_type, arc_id, q
      ))
    )
    mta_enr[, log_enrichment := log2(enrichment)]
    mta_enr[, motif_pair := sprintf("%s + %s", motif1, motif2)]
    
    # Load motif PWM data
    rv$mta_pwm <- read_meme(file.path(dat_dir, "motif-archetypes-all.meme"))
    names(rv$mta_pwm) <- sapply(rv$mta_pwm, function(m) m@name)
    
    # Load assignment data
    assign_dt <- fread(file.path(
      dat_dir,
      sprintf("motif-assignment-archetypes-%s.tsv.gz", arc_id)
    ))
    assign_dt <- assign_dt[, .(archetype_name, gene, gene_name, og, common_name, tf_family)]
    assign_dt[common_name != "", gene_name := common_name][, common_name := NULL]
    assign_dt[gene_name == "", gene_name := ifelse(nchar(og) > 40, paste0(substr(og, 1, 37), "..."), og)][, og := NULL]
    assign_dt <- assign_dt[, lapply(.SD, paste, collapse = ";"), .SDcols = c("gene", "gene_name", "tf_family"), by = archetype_name]
    
    # Process enrichment data
    mta_enr_gen <- copy(mta_enr)
    setnames(assign_dt, c("motif1", "gene_motif1", "gene_name_motif1", "tf_family_motif1"))
    mta_enr_gen <- merge.data.table(mta_enr_gen, assign_dt, by = "motif1", allow.cartesian = TRUE)
    setnames(assign_dt, c("motif2", "gene_motif2", "gene_name_motif2", "tf_family_motif2"))
    mta_enr_gen <- merge.data.table(mta_enr_gen, assign_dt, by = "motif2", allow.cartesian = TRUE)
    mta_enr_gen[, gene_pair := sprintf("%s + %s", gene_name_motif1, gene_name_motif2)]
    
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
    ggseqlogo(rv$mta_pwm[[selected$motif1]]@motif) + theme_void()
  })
  
  output$motif2_logo <- renderPlot({
    event_data <- event_data("plotly_click")
    req(event_data, rv$mta_enr_gp)
    selected <- rv$mta_enr_gp[event_data$pointNumber + 1]
    ggseqlogo(rv$mta_pwm[[selected$motif2]]@motif) + theme_void()
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
    sel_mta_enr_cts <- mta_enr_cts[motif1 == selected$motif1 & motif2 == selected$motif2]
    ggplot(sel_mta_enr_cts, aes(cell_type, enrichment, fill = cell_type)) + 
      geom_bar(stat = "identity", color = "black") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = ct_cols) +
      labs(
        x = "Cell type",
        y = "Enrichment (log2)",
        title = "Motif pair enrichment across cell types"
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none"
      )
  })
  
}

shinyApp(ui = ui, server = server)
