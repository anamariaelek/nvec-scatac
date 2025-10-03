
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

shiny_env = new.env()
shinyServer(function(input, output, session) {

  atac2DServer(id = "scatac", config_id = "scatac")
  atac2DServer(id = "adult", config_id = "adult")
  atac2DServer(id = "gastrula", config_id = "gastrula")
  geneSummaryServer(id = "gene_summary", config_id="scatac")
  peaksServer(id = "peaks", config_id="scatac")
  tfServer(id = "tfs", config_id="scatac")
  cellTreeServer(id = "cell_type_trees", config_id = "scatac")
  
})

