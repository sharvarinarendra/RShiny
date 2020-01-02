library(shiny)
library(DESeq2)
library(pacman)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(genefilter)
library(plyr)
library(rowr)
library(gplots)
library(VennDiagram)

#source("helpers.R")

PCA <- function(rld, dds, title) {
  
  plotPCA(rld, intgroup = "group") + geom_label(aes(label = dds$ID)) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  
}

varClust <- function(rld, n, title) {
  
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), n)
  mat <- assay(rld)[topVarGenes,]
  
  mat <- mat - rowMeans(mat)
  
  anno <- as.data.frame(colData(rld)[, c("group")])
  rownames(anno) <- colData(rld)[, c("ID")]
  colnames(anno) <- "group"
  pheatmap(mat, annotation_col = anno, main = title)
}

countPlot <- function(gene, res_ann, dds) {
  x <- gene
  y <- as.character(str_remove_all(x, "\\s-(.*)"))
  top <- res_ann[y, c("pvalue", "padj", "symbol")]
  cnts <- plotCounts(dds, gene = y, intgroup = "group",
                     returnData = TRUE)
  print(ggplot(cnts, aes(x = group, y = count, color = group)) +
          geom_label(aes(label = dds$ID),
                     position = position_jitter(width = 0.1, height = 0)) +
          ggtitle(paste0(y, " (", top$symbol, ")"),
                  subtitle = paste0("p = ", sprintf("%.2e", top$pvalue),
                                    "; fdr = ", sprintf("%.2e", top$padj))) +
          theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)))
}

test_res <- res.gene.df
test_res$ensembl <- rownames(test_res)
test_res$new <- paste(test_res$ensembl, "-", test_res$symbol)

test_res_1 <- res.transcript.df
test_res_1$ensembl <- rownames(test_res_1)
test_res_1$new <- paste(test_res_1$ensembl, "-", test_res_1$symbol)

res.gene.sva.df$ensembl <- row.names(res.gene.sva.df)
res.gene.sva.df.1$ensembl <- row.names(res.gene.sva.df.1)
res.gene.sva.df.2$ensembl <- row.names(res.gene.sva.df.2)
res.transcript.sva.df$ensembl <- row.names(res.transcript.sva.df)
res.transcript.sva.df.1$ensembl <- row.names(res.transcript.sva.df.1)
res.transcript.sva.df.2$ensembl <- row.names(res.transcript.sva.df.2)

choice_list <- list("Acute Vs Water Upregulated Genes",
                 "Acute Vs Water Downregulated Genes",
                 "Chronic Vs Water Upregulated Genes",
                 "Chronic Vs Water Downregulated Genes",
                 "Acute Vs Chronic Upregulated Genes",
                 "Acute Vs Chronic Downregulated Genes")

# Define UI ----
ui <- fluidPage(
  titlePanel("Jung's RNA-Seq Analysis Visualization"), # can also change width of fluid row and plot the plots
  
  mainPanel(img(src = "mclean_logo.png", height = 60, width = 280)),
             
  
  br(),
 
    fluidRow(
      column(6, wellPanel(
          h4('PCA plot'),
          selectInput('in1', 
                      label = "Choose data to analyze",
                      choices = list("Gene", 
                                     "Transcript")
                      ),
          plotOutput('out1')
    )
    ),
    column(5, wellPanel(
           h4('Checking the significant and up/downregulated genes'),
           selectInput('in2', 
                       label = "Choose a dataset : ",
                       choices = list("Acute Vs Water (Gene)", "Acute Vs Water (Transcript)", "Chronic Vs Water (Gene)", 
                                      "Chronic Vs Water (Transcript)", "Acute Vs Chronic (Gene)", 
                                      "Acute Vs Chronic (Transcript)")
           ),
           numericInput('in3',
                       label = "Log2 Fold Change Cut-Off for Upregulation (should be greater than zero)",
                       value = 1),
           numericInput('in4',
                        label = "Log2 Fold Change Cut-Off for Downregulation (should be less than zero)",
                        value = -1),
           helpText("Note : You can also input log2 fold change as 0 for both up and downregulation. The result will simply 
                    show the number of genes/transcripts greater than 0 (upregulation) and less than 0 (downregulation)"),
           numericInput('in5',
                        label = "Choose p-value cut-off",
                        value = 0.05),
           actionButton("submit", "Submit"),
           DT::dataTableOutput("out2")
    ) 
          ),
    column(6, wellPanel(
           h4('Heatmap showing top genes/transcripts with highest variance across the samples (Clustering)'),
           radioButtons('in6',
                       label = "Choose data : ",
                       choices = list("Gene" = 1, "Transcript" = 2)
                       ),
           numericInput('in7',
                        label = "Enter the number of genes/transcripts to display : ",
                        value = 20
                        ),
           plotOutput('out3')
           )
           ),
    column(6, wellPanel(
           h4('Count Plot of Gene'),
           selectizeInput('in8',
                        label = "Choose Gene ID : ",
                        test_res$new, options = list(create = TRUE)
           ),
           plotOutput('out4')
    )
    ),
    column(6, wellPanel(
           h4('Count Plot of Transcript'),
           selectizeInput('in9',
                       label = "Choose Transcript ID : ",
                       test_res_1$new, options = list(create = TRUE)
                       ),
           plotOutput('out5')
           )
    ),
    column(6, wellPanel(
           h4('Venn Diagram'),
           h5('Choose any four comparisons from the following : '),
           selectInput('in10',
                       label = "Comparison 1 : ",
                       choice_list
           ),
           selectInput('in11',
                       label = "Comparison 2 : ",
                       choice_list
           ),
           selectInput('in12',
                       label = "Comparison 3 : ",
                       choice_list
           ),
           selectInput('in13',
                       label = "Comparison 4 : ",
                       choice_list
           ),
           numericInput('in14',
                        label = "Log2 Fold Change Cut-Off for Upregulation (should be greater than zero)",
                        value = 0),
           numericInput('in15',
                        label = "Log2 Fold Change Cut-Off for Downregulation (should be less than zero)",
                        value = 0),
           helpText("Note : You can also input log2 fold change as 0 for both up and downregulation. The result will simply 
                    show the number of genes/transcripts greater than 0 (upregulation) and less than 0 (downregulation)"),
           numericInput('in16',
                        label = "Choose p-value cut-off",
                        value = 0.05),
           actionButton("submit1", "Submit"),
           plotOutput('out6')
    )
    )
  )
)


  
# Define server logic ----
server <- function(input, output, session) {

  output$out1 <- renderPlot({
    if (input$in1 == "Gene") {
      PCA(rld.gene, dds.gene1, "gene level PCA of samples")
    }
    else if (input$in1 == "Transcript"){
      PCA(rld.transcript, dds.transcript1, "transcript level PCA of samples")
    }
    
  })
  
  datasetInput <- eventReactive(input$submit, {
    switch(input$in2,
           "Acute Vs Water (Gene)" = res.gene.ann,
           "Acute Vs Water (Transcript)" = res.transcript.ann,
           "Chronic Vs Water (Gene)" = res.gene.ann.1,
           "Chronic Vs Water (Transcript)" = res.transcript.ann.1, 
           "Acute Vs Chronic (Gene)" = res.gene.ann.2, 
           "Acute Vs Chronic (Transcript)" = res.transcript.ann.2
           )
  })

  sig <- eventReactive(input$submit, {
    gene <- as.data.frame(datasetInput())
    gene1 <- gene %>% filter(log2FoldChange > input$in3) %>% filter(pvalue < input$in5) %>% filter(!(is.na(symbol))) %>% dplyr::select(symbol)
    if (nrow(gene1) == 0) {
      gene1 <- rbind(gene1, "NA")
      gene1$num <- 1
      colnames(gene1) <- c("Upregulated", "num")
    } else if (nrow(gene1) > 0) {
      colnames(gene1) <- "Upregulated"
      gene1$num <- 1:nrow(gene1)
    }
    gene2 <- gene %>% filter(log2FoldChange < input$in4) %>% filter(pvalue < input$in5) %>% filter(!(is.na(symbol))) %>% dplyr::select(symbol)
    if (nrow(gene2) == 0) {
      gene2 <- rbind(gene2, "NA")
      gene2$num <- 1
      colnames(gene2) <- c("Downregulated", "num")
    } else if (nrow(gene2) > 0) {
      colnames(gene2) <- "Downregulated"
      gene2$num <- 1:nrow(gene2)
    }
    df <- full_join(gene1, gene2)
    df <- df[,-2]
    return(print(df))
  })
  
  output$out2 <- DT::renderDataTable({
    sig()
  })
  
  output$out3 <- renderPlot({
    if (input$in6 == 1) {
      varClust(rld.gene, input$in7, "Gene Level Heatmap")
    }
    else {
      varClust(rld.transcript, input$in7, "Transcript level Heatmap")
    }
  })
  
  output$out4 <- renderPlot({
      countPlot(input$in8, res.gene.ann, dds.gene)
  })
  
  output$out5 <- renderPlot({
    countPlot(input$in9, res.transcript.ann, dds.transcript)
  }) 
  
  sw_in <- reactive(
    "Acute Vs Water Upregulated Genes" = res.gene.sva.df %>% filter(log2FoldChange > input$in14) %>% 
      filter(pvalue < input$in16),
    "Acute Vs Water Downregulated Genes" = res.gene.sva.df %>% filter(log2FoldChange < input$in15) %>% 
      filter(pvalue < input$in16),
    "Chronic Vs Water Upregulated Genes" = res.gene.sva.df.1 %>% filter(log2FoldChange > input$in14) %>% 
      filter(pvalue < input$in16),
    "Chronic Vs Water Downregulated Genes" = res.gene.sva.df.1 %>% filter(log2FoldChange < input$in15) %>% 
      filter(pvalue < input$in16),
    "Acute Vs Chronic Upregulated Genes" = res.gene.sva.df.2 %>% filter(log2FoldChange > input$in14) %>% 
      filter(pvalue < input$in16),
    "Acute Vs Chronic Downregulated Genes" = res.gene.sva.df.2 %>% filter(log2FoldChange < input$in15) %>% 
      filter(pvalue < input$in16)  
  )
  
  datasetInput1 <- eventReactive(input$submit1, {
    switch(input$in10,
          sw_in()
    )
  })
  
  datasetInput2 <- eventReactive(input$submit1, {
    switch(input$in11,
           sw_in()
    )
  })
  
  datasetInput3 <- eventReactive(input$submit1, {
    switch(input$in12,
           sw_in()
    )
  })
  
  datasetInput4 <- eventReactive(input$submit1, {
    switch(input$in13,
           sw_in()
    )
  })
  
  output$out6 <- renderPlot({
    venn.plot <-  venn.diagram(
      x = list(datasetInput1()$ensembl, datasetInput2()$ensembl, datasetInput3()$ensembl, datasetInput4()$ensembl),
      NULL,
      category.names = c("Comp 1", "Comp 2", "Comp 3", "Comp 4"), 
      fill=c(alpha("blue",0.3), alpha('light blue',0.3), alpha('red',0.3), alpha("pink",0.3)),
      cex = 1.5, cat.fontface = 4
    )
    
    grid.draw(venn.plot)
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
