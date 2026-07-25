library(shiny)
library(ggplot2)
library(umap)
library(scales)
library(RomicsProcessor)

# Define UI for app
ui <- fluidPage(
  # App title
  titlePanel("RomicsExplorer"),

  # Tabs
  tabsetPanel(
    # First tab for loading data
    tabPanel("Load Data",
             fileInput("datafile", "Choose a romics_object to be loaded",
                       accept = c("text/csv", "text/comma-separated-values,text/plain")),
             br(),
             actionButton("load", "Load Data"),
             textOutput("load_status")
    ),

    # Second tab for grouping figure
    tabPanel("Grouping Figure",
             selectInput("group_var", "Select a grouping variable:",
                         choices = c("", "Group1", "Group2", "Group3")),
             br(),
             plotOutput("group_plot")
    ),

    # Third tab for boxplots
    tabPanel("Boxplots",
             selectInput("gene", "Select a gene:",
                         choices = c("", "Gene1", "Gene2", "Gene3")),
             br(),
             plotOutput("boxplot")
    ),

    # Fourth tab for volcano plot
    tabPanel("Volcano Plot",
             br(),
             plotOutput("volcano")
    )
  )
)

# Define server function
server <- function(input, output) {

  # Load data function
  data <- reactive({
    req(input$datafile)
    read.csv(input$datafile$datapath)
  })

  # Load status output
  output$load_status <- renderText({
    if (input$load == 0) {
      "Waiting for data to be loaded..."
    } else {
      "Data loaded successfully!"
    }
  })

  # Grouping plot output
  output$group_plot <- renderPlot({
    req(data())
    ggplot(data(), aes(x = umap_x, y = umap_y, color = input$group_var)) +
      geom_point() +
      theme_bw() +
      scale_color_discrete(name = "Grouping Variable") +
      labs(x = "UMAP1", y = "UMAP2", title = "UMAP Plot")
  })

  # Boxplot output
  output$boxplot <- renderPlot({
    req(data(), input$gene)
    ggplot(data(), aes(x = input$gene, y = expression_value, fill = input$group_var)) +
      geom_boxplot() +
      theme_bw() +
      scale_fill_discrete(name = "Grouping Variable") +
      labs(x = "Gene", y = "Expression Value", title = "Boxplot")
  })

  # Volcano plot output
  output$volcano <- renderPlot({
    req(data())
    ggplot(data(), aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(abs(log2FoldChange) > 1 & pvalue < 0.05, "Significant", "Not Significant")), size = 1) +
      theme_bw() +
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
      labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot")
  })

}

# Run the app
shinyApp(ui = ui, server = server)
