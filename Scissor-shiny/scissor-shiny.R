library(shiny)

ui <- tagList(
  fluidPage(
    titlePanel("Scissor"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        # 对应后面的 output$file2
        uiOutput('file2'),
        # 对应后面的 output$file3
        uiOutput('file3'),
        
        actionButton('reset', 'RESET'),
        hr(),
        radioButtons("Type", "Choose: method",
                     c("COX"="c",'Logistic regression'='l')),
        hr(),
        downloadButton("downloadData", "Download"),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
        br(),
        h5('Cition:'),
        h6('Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data')
      ),
      mainPanel(
        h4("Scissor plot"),
        br(),
        br(),
        shinycssloaders::withSpinner(
          plotOutput("plot")
        )
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  # 上传文件后自动读取文件
  dataInput1 <- reactive({
    sessionEnvir <- sys.frame()
    if (!is.null(input$file1)) eval(parse(text = load(input$file1$datapath, sessionEnvir)))
  })
  
  # 上传文件后自动读取文件
  dataInput2 <- reactive({
    sessionEnvir <- sys.frame()
    if (!is.null(input$file2)) eval(parse(text = load(input$file2$datapath, sessionEnvir)))
  })
  
  # 上传文件后自动读取文件
  dataInput3 <- reactive({
    sessionEnvir <- sys.frame()
    if (!is.null(input$file3)) eval(parse(text = load(input$file3$datapath, sessionEnvir)))
  })
  
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose scRNA RData",
                accept=c('.RData, .Rds')
      )
    })
  }, ignoreNULL = F)
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Step 2: Choose bulk RNA RData",
                accept=c('.RData, .Rds')
      )
    })
  }, ignoreNULL = F)
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file3 <- renderUI({
      fileInput("file3", "Step 3: Choose phenotype RData",
                accept=c('.RData, .Rds')
      )
    })
  }, ignoreNULL = F)
  

  
  datasetInput_cox <- reactive({
    library(Scissor)
    
    # 三个文件上传完自动触发运行
    sc_dataset = dataInput1()
    bulk_dataset = dataInput2()
    phenotype = dataInput3()
    
    
    sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)
    colnames(phenotype) = c('V1','V2','V3')
    
    
    all(colnames(bulk_dataset) == phenotype$V1)
    phenotype <- phenotype[,2:3]
    colnames(phenotype) <- c("time", "status")
    
    infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.05, 
                      family = "cox", Save_file = 'Scissor_inputs.RData')
    
    Scissor_select <- rep(0, ncol(sc_dataset))
    names(Scissor_select) <- colnames(sc_dataset)
    Scissor_select[infos1$Scissor_pos] <- 1
    Scissor_select[infos1$Scissor_neg] <- 2
    sc_dataset_1 <<- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
  })
  
  datasetInput_lr <- reactive({
    library(Scissor)
    
    # 三个文件上传完自动触发运行
    sc_dataset = dataInput1()
    bulk_dataset = dataInput2()
    phenotype = dataInput3()
    
    sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)
    phenotype <- phenotype
    tag <- c('wild-type', 'TP53 mutant')
    infos4 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = 0.5, 
                      family = "binomial", Save_file = "Scissor_LUAD_TP53_mutation.RData")
    
    infos5 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = 0.2, 
                      family = "binomial", Load_file = "Scissor_LUAD_TP53_mutation.RData")
    
    Scissor_select <- rep(0, ncol(sc_dataset))
    names(Scissor_select) <- colnames(sc_dataset)
    Scissor_select[infos5$Scissor_pos] <- 1
    Scissor_select[infos5$Scissor_neg] <- 2
    sc_dataset_2 <<- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
  })
  
  output$plot <- renderPlot({
    
    if(input$Type == 'c'){
      datasetInput_cox()
      
      DimPlot(sc_dataset_1, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'),
              pt.size = 1.2, order = c(2,1))
    } else if(input$Type == 'l'){
      datasetInput_lr()
      
      DimPlot(sc_dataset_2, reduction = 'umap', group.by = 'scissor', 
              cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      if (input$Type == 'c') {
        write.csv(sc_dataset_1@meta.data,file,row.names = T,quote = F)
      } else if (input$Type == 'l') {
        write.csv(sc_dataset_2@meta.data,file,row.names = T,quote = F)
      }
    }
  )
 
   
}

# Run the application 
shinyApp(ui = ui, server = server)

