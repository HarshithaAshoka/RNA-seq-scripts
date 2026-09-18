library(shiny)
library(readxl)
library(dplyr)
library(VennDiagram)
library(openxlsx)
library(grid)

# Increase maximum upload size to 50 MB
options(shiny.maxRequestSize = 50*1024^2)  

ui <- fluidPage(
  titlePanel("DESeq2 Contrast Venn Diagram (Up/Downregulated)"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload Contrast 1 Excel"),
      textInput("name1", "Contrast 1 Name", "Contrast1"),
      fileInput("file2", "Upload Contrast 2 Excel"),
      textInput("name2", "Contrast 2 Name", "Contrast2"),
      
      numericInput("padj_thresh", "padj threshold", 0.05, min=0, step=0.01),
      numericInput("logfc_thresh", "log2FC threshold (absolute)", 0.585, min=0, step=0.01),
      
      numericInput("label_pos1", "Label position 1", 0, step=1),
      numericInput("label_pos2", "Label position 2", 0, step=1),
      numericInput("label_dist1", "Label distance 1", 0.05, step=0.01),
      numericInput("label_dist2", "Label distance 2", 0.05, step=0.01),
      
      radioButtons("regulation", "Regulation type", choices = c("All", "Up", "Down"), selected = "All"),
      
      downloadButton("download_venn", "Download Venn Diagram"),
      downloadButton("download_excel", "Download DEG Excel Sheets")
    ),
    
    mainPanel(
      plotOutput("vennPlot")
    )
  )
)

server <- function(input, output, session) {
  
  # helper to filter DEGs and keep only desired columns
  get_DEGs <- function(file, padj_thresh, logfc_thresh, regulation="All"){
    req(file)
    df <- read_excel(file$datapath)
    
    # keep only relevant columns if they exist
    keep_cols <- c("ensembl", "entrezgene", "gene_biotype", 
                   "chromosome_name", "gene_symbol", 
                   "log2FoldChange", "pvalue", "padj")
    df <- df %>% select(any_of(keep_cols))
    
    # filter DEGs
    df <- df %>%
      filter(!is.na(padj)) %>%
      filter(padj < padj_thresh) %>%
      filter(abs(log2FoldChange) > logfc_thresh)
    
    if(regulation == "Up") df <- df %>% filter(log2FoldChange > logfc_thresh)
    if(regulation == "Down") df <- df %>% filter(log2FoldChange < -logfc_thresh)
    
    df
  }
  
  deg1 <- reactive({ get_DEGs(input$file1, input$padj_thresh, input$logfc_thresh, input$regulation) })
  deg2 <- reactive({ get_DEGs(input$file2, input$padj_thresh, input$logfc_thresh, input$regulation) })
  
  venn_sets <- reactive({
    list(
      set1 = deg1()$gene_symbol,
      set2 = deg2()$gene_symbol
    )
  })
  
  output$vennPlot <- renderPlot({
    req(venn_sets())
    sets <- venn_sets()
    
    size1 <- length(sets$set1)
    size2 <- length(sets$set2)
    cross_size <- length(intersect(sets$set1, sets$set2))
    
    venn.plot <- draw.pairwise.venn(
      area1 = size1,
      area2 = size2,
      cross.area = cross_size,
      category = c(input$name1, input$name2),
      fill = c("skyblue", "pink"),
      lty = "blank",
      ind = FALSE,
      scaled = TRUE,
      cex = 1.5,
      cat.pos = c(input$label_pos1, input$label_pos2),
      cat.dist = c(input$label_dist1, input$label_dist2)
    )
    grid.draw(venn.plot)
  })
  
  output$download_venn <- downloadHandler(
    filename = function() { paste0("VennDiagram_", input$regulation, ".pdf") },
    content = function(file){
      pdf(file)
      sets <- venn_sets()
      size1 <- length(sets$set1)
      size2 <- length(sets$set2)
      cross_size <- length(intersect(sets$set1, sets$set2))
      
      venn.plot <- draw.pairwise.venn(
        area1 = size1,
        area2 = size2,
        cross.area = cross_size,
        category = c(input$name1, input$name2),
        fill = c("skyblue", "pink"),
        lty = "blank",
        ind = FALSE,
        scaled = TRUE,
        cex = 1.5,
        cat.pos = c(input$label_pos1, input$label_pos2),
        cat.dist = c(input$label_dist1, input$label_dist2)
      )
      grid.draw(venn.plot)
      dev.off()
    }
  )
  
  output$download_excel <- downloadHandler(
    filename = function() { paste0("DEG_Venn_Sets_", input$regulation, ".xlsx") },
    content = function(file){
      set1 <- deg1()
      set2 <- deg2()
      
      only1 <- set1 %>% filter(!gene_symbol %in% set2$gene_symbol)
      only2 <- set2 %>% filter(!gene_symbol %in% set1$gene_symbol)
      
      # Shared: keep both contrasts’ values, side-by-side
      shared <- inner_join(set1, set2, by = c("gene_symbol", "ensembl", "entrezgene", 
                                              "gene_biotype", "chromosome_name"),
                           suffix = c(paste0("_", input$name1), paste0("_", input$name2)))
      
      wb <- createWorkbook()
      addWorksheet(wb, paste0(input$name1, "_only"))
      addWorksheet(wb, paste0(input$name2, "_only"))
      addWorksheet(wb, "shared")
      
      writeData(wb, sheet=1, only1)
      writeData(wb, sheet=2, only2)
      writeData(wb, sheet=3, shared)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
