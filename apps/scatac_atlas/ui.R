
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(dashboardPage(
  
  dashboardHeader(title = "scATAC"),
  
  dashboardSidebar(
    disable = FALSE,
    sidebarMenu(
      id="tabs",
      menuItem("Joint atlas", tabName = "scatac"),
      menuItem("Adult", tabName = "adult"),
      menuItem("Gastrula", tabName = "gastrula"),
      menuItem("Genes summary", tabName = "gene_summary"),
      menuItem("Peaks summary", tabName = "peaks"),
      menuItem("TF summary", tabName = "tfs"),
      menuItem("Cell type trees", tabName = "cell_type_trees")
    )
  ), 
  
  dashboardBody(
    
    tabItems(
      tabItem(tabName = "scatac", atac2DUI(id = "scatac")),
      tabItem(tabName = "adult", atac2DUI(id = "adult")),
      tabItem(tabName = "gastrula", atac2DUI(id = "gastrula")),
      tabItem(tabName = "gene_summary", geneSummaryUI(id = "gene_summary")),
      tabItem(tabName = "peaks", peaksUI(id = "peaks")),
      tabItem(tabName = "tfs", tfUI(id = "tfs")),
      tabItem(tabName = "cell_type_trees", cellTreeUI(id = "cell_type_trees"))
    )
  )
  
))


