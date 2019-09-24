Sys.umask("002") # this ensures that all files created during the current session can be edited and deleted by all members of the group

# !!! http://shiny.rstudio.com/gallery/update-input-demo.html

library(shiny)
library(shinythemes) # for the tab layout of the website
library(data.table) # for data manipulation
library(utils) # for unzip
library(RJDBC) # to access the databases
library(xtable) # to display tables
require(XLConnect) # to export results in excel
options(shiny.maxRequestSize=50*1024^2) 
options(scipen = 999) ### to turn off scientific notation; to turn on - options(scipen = 0)

path2data <- file.path("","opt","app","data", "mra")
#path2data <- file.path("data")
source("lib_validation.R")

ui <- navbarPage(
  "Assessment of RMS Sample Micro Representativeness"
  , theme = shinytheme('cerulean')
  
  #######################################################################################
  #######################################################################################
  ### file manager tab
  
  , tabPanel(
    "File Manager"
    , sidebarPanel(
      uiOutput("new_folder")
      , actionButton("create_folder_button","Create dataset")
      , uiOutput("filemanager_choose_folder")
      , actionButton("delete_folder_button","Delete dataset")
    )
    , mainPanel(
      checkboxInput("unzip","Unzip uploaded files",value=TRUE)
      , fileInput("upload","Upload")
      , uiOutput("choose_files")
      , actionButton("show_button","Show")
      , downloadButton("download_files")
      , actionButton("delete_button","Delete")
      , textOutput("delete_files")
      , tableOutput("show_files")
    )
  )
  
  #######################################################################################
  #######################################################################################
  ### validation tab
  
  , tabPanel(
    "Input Validation"
    , htmlOutput("validation_cur_folder")
    , actionButton("validation_run_button","Validate selected")
    , h4("Explanation of Valid flag:")
    , h5("1 = Valid.")
    , h5("2 = More than 5% stores in panel and / or 20% stores in universe have missing information for secondary store characteristic.")
    , h5("3 = Secondary store characteristic value present in less than 5% stores in universe.")
    , tableOutput("validation_log")
    , dataTableOutput("validation_codebook")
  )
  
  
  #######################################################################################
  #######################################################################################
  ### variable selection tab
  
  , tabPanel(
    "Variable Selection"
    , sidebarPanel(
      renderUI("varsel_cur_folder")
      , actionButton("varsel_run_button","Run variable selection")
      , width = 2
    )
    , mainPanel(
      textOutput("varsel_anova")
      , tableOutput("varsel_scores")
      , dataTableOutput("varsel_codebook")
    )
  )
  
  #######################################################################################
  #######################################################################################
  ### sample check tab
  
  , tabPanel(
    "Micro Rep Analysis"
    , sidebarPanel(
      p(strong("Active dataset:")),
      textOutput("samplecheck_cur_folder"),
      tags$head(tags$style("#samplecheck_cur_folder{color: red;}")),
      
      hr(),
      actionButton("goButton_1", "Generate MBD Importance Summary"),
      
      
      conditionalPanel(
        condition = "output.status == 'TRUE'",
        
        hr(),
        radioButtons('min_retain_ss', 'Minimum Retain Sample Size',
                     c("Minimum of Current Sample Size and Target Sample Size" = "Current.SS", "Zero" = "0"),
                     selected = "Current.SS"),
        hr(),
        
        p("Select the ",strong("Secondary Characteristic"), " (for sub-cell creation) and ", strong("MBD"), " (for sample rotation). 
          If no MBD is selected, sample rotation will be executed as per recommendation from tool."),
        
        uiOutput("choose_sec_char"),
        uiOutput("choose_sample_rotation_mbd"),
        
        hr(),
        
        strong("P-value:"),
        numericInput("pvalue",
                     label = h5("Adjust sub-cell allocation when p-value of Chi-square test is less than: "), 
                     value = 0.05),
        hr(),
        
        actionButton("goButton", "Run Micro Rep Analysis"),
        
        hr()),
      
      hr(),
      
      conditionalPanel(
        condition = "output.status == 'TRUE'",
        h4("Download Results"),
        downloadButton('downloadData', 'Result'),
        helpText("Download results for one secondary characteristics variable at a time.")
      ))
    
    , mainPanel(
      conditionalPanel(
        condition = "output.status == TRUE",
        tabsetPanel( 
          
          tabPanel("MBD Importance Summary",dataTableOutput("output_mbdreport"), hr(),
                   h4(""), htmlOutput("remark"),
                   tags$style(type="text/css", "tfoot {display:none;}")),
          
          tabPanel("MBD Importance Analysis",
                   h2(""), htmlOutput("note_sec1"), hr(),
                   dataTableOutput("output_kpireport"), hr(),
                   h4(""), htmlOutput("remark2"),
                   tags$style(type="text/css", "tfoot {display:none;}")),
          
          tabPanel("Net Sample Compliance",
                   h3(""), htmlOutput("note_sec2"), hr(),
                   h4("Chi Square Test"),dataTableOutput("output_chisq1"), hr(),
                   h3(""), htmlOutput("note_sec3"), hr(),
                   h4("Subcell Allocation"),dataTableOutput("output_subcell1"),
                   tags$style(type="text/css", "tfoot {display:none;}")),
          
          tabPanel("Micro Rep and Sample Rotation", h3(""), htmlOutput("note"), hr(),
                   h3(""), htmlOutput("note_sec4"), hr(),
                   h4("Chi Square Test"), dataTableOutput("output_chisq2"),hr(),
                   h3(""), htmlOutput("note_sec5"), hr(),
                   h4("Subcell Allocation"), dataTableOutput("output_subcell2"),
                   tags$style(type="text/css", "tfoot {display:none;}"))
          
          
        )
      )
    )
    )
)



#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################



server <- function(input, output, session){
  
  options(shiny.maxRequestSize=120*1024^2) # to increase maximum upload file size
  
  #######################################################################################
  #######################################################################################
  ### file manager tab
  
  output$new_folder <- renderUI({
    textInput("new_folder","Name of a new dataset")
  })
  
  values <- reactiveValues(create_folder = FALSE)
  
  observe({
    if (input$create_folder_button == 0) return(NULL)
    values$create_folder = TRUE
  })
  
  create_folder <- reactive({
    if (!values$create_folder) return(NULL)
    dir.create(file.path(path2data,input$new_folder))
    updateTextInput(session, "new_folder","Enter Column Name",'')
    values$create_folder = FALSE
    return(list(result = "done"))
  })
  
  output$filemanager_choose_folder <- renderUI({
    delete_folder()
    create_folder()
    list_folders <- list.dirs(path2data,full.names=FALSE)[-1]
    radioButtons(
      "filemanager_folder"
      , label = "Choose dataset"
      , choices  = list_folders
    )
  })
  
  delete_folder <- reactive({
    input$delete_folder_button
    if (input$delete_folder_button[1] == 0) return(NULL)
    isolate({
      unlink(file.path(path2data,"32456"),recursive=TRUE)
      selected <- input$filemanager_folder
      for (folder in selected) { unlink(file.path(path2data,folder),recursive=TRUE) }
      return("folder deleted")
    })
  })  
  ### files ############################################
  
  output$choose_files <- renderUI({
    delete_files()
    upload_files()
    create_folder()
    list_files <- list.files(file.path(path2data,input$filemanager_folder))
    checkboxGroupInput(
      "files"
      , "choose files"
      , choices  = list_files
      , selected = NULL
    )
  })
  
  upload_files <- reactive({
    inFile <- input$upload
    if (is.null(inFile)) return(NULL)
    isolate({
      unlink(file.path(path2data,input$filemanager_folder,inFile$name))
      file.copy(inFile$datapath, file.path(path2data,input$filemanager_folder,inFile$name))
      if (input$unzip & substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name)) == ".zip")
      {
        unzip(file.path(path2data,input$filemanager_folder,inFile$name),exdir = file.path(path2data,input$filemanager_folder))
        unlink(file.path(path2data,input$filemanager_folder,inFile$name))
      }
      return(list(result = "done"))
    })
  })
  
  output$show_files <- renderUI({
    input$show_button
    if (input$show_button[1] == 0) return(NULL)
    isolate({
      selected <- input$files # selected <- c("codebook.csv","sales_clean.csv") ; path2data <- file.path("","mnt","data","mra_data")
      tables <- list()
      for (fl in selected)
      {
        table <- fread(file.path(path2data,input$filemanager_folder,fl))[1:10,] 
        #table <- fread(file.path(path2data,"!_test_01",fl))[1:10,]
        tables[[as.character(fl)]] <- 
          print(
            xtable(table, caption=fl)
            , type="html"
            , html.table.attributes='class="data table table-bordered table-condensed"'
            ,caption.placement="top"
          )
      }
      big_table <- lapply(tables,paste)
      #return(div(HTML(big_table),class="shiny-html-output"))
      return(div(lapply(big_table,HTML),class="shiny-html-output"))
    })
  })
  
  delete_files <- reactive({
    input$delete_button
    if (input$delete_button[1] == 0) return(NULL)
    isolate({
      selected <- input$files
      for (fl in selected) { file.remove(file.path(path2data,input$filemanager_folder,fl)) }
      return("files deleted")
    })
  })
  
  output$download_files <- downloadHandler(
    filename = function() { paste(paste("download_",format(Sys.time(), format = "%Y%m%d-%H%M%S"),".zip")) },
    content = function(file)
    {
      selected <- input$files
      unlink("tmp", recursive = TRUE, force = TRUE)
      dir.create("tmp")
      tempdir <- strsplit(tempdir(),"/")[[1]][3]
      dir.create(file.path("tmp",tempdir))
      for (fl in selected) { file.copy(file.path(path2data,input$filemanager_folder,fl),file.path("tmp",tempdir,fl)) }
      zip(zip = "zip", zipfile = file, files = file.path("tmp",tempdir,selected),extras = "-j")
    }
  )
  
  #######################################################################################
  #######################################################################################
  ### validation tab
  
  output$validation_cur_folder <- renderUI({
    
    status <- "Non-validated"
    if (isTRUE(file.exists(file.path(path2data,input$filemanager_folder,"validation_log.csv")))) { status <- "Validated" }
    
    HTML(
      paste0(
        "<h4>Active dataset: "
        , "<font color = red>"
        , input$filemanager_folder
        , "</font>"
        , " ("
        , status
        , ")"
        , "</h4>"
        , "<hr>"
      )
    )
  })
  
  output$validation_log <- renderTable({
    validate()
    if (!isTRUE(file.exists(file.path(path2data,input$filemanager_folder,"validation_log.csv")))) return(NULL)
    return(fread(file.path(path2data,input$filemanager_folder,"validation_log.csv")))
  })
  
  output$validation_codebook <- renderDataTable({
    validate()
    if (!isTRUE(file.exists(file.path(path2data,input$filemanager_folder,"codebook.csv")))) return(NULL)
    codebook.dt <- fread(file.path(path2data,input$filemanager_folder,"codebook.csv"))
    # codebook.dt <- fread(file.path(path2data,"CHINA_CVS","codebook.csv"))
    setnames(
      codebook.dt,
      c("mbd","variable","value","res_freq","res_proportion","panel_freq","panel_proportion","valid")
      , c("MBD","Secondary Characteristic","Secondary Characteristic Value","Count (RES)", "Proportion (RES)", "Count (Panel)", "Proportion (Panel)", "Valid flag")
    )
    return(codebook.dt)
  })
  
  validate <- reactive({
    if (input$validation_run_button[1] == 0) return(NULL)
    isolate({
      path <- file.path(path2data,input$filemanager_folder)
      
      # code book
      create_codebook_f(path)
      
      # log
      validate_f(path)
    })
  })
  
  #######################################################################################
  #######################################################################################
  ### variable selection tab
  
  output$varsel_cur_folder <- renderText(return(paste0("Selected dataset: ",input$filemanager_folder)))
  
  build_scores <- reactive({
    #anova()
    if (!file.exists(file.path(path2data,input$filemanager_folder,"codebook.csv"))) return(NULL)
    codebook.DT <- fread(file.path(path2data,input$filemanager_folder,"codebook.csv")) # codebook.DT <- fread(file.path(path2data,"COUNTRY - INDEX - PERIOD","codebook.csv"))
    setkey(codebook.DT,mbd,variable)
    if (!"score" %in% colnames(codebook.DT)) {return(NULL)}
    scores.DT <- dcast.data.table(unique(codebook.DT[,.(mbd,variable,score)]), variable ~ mbd, value.var = "score")
    return(list(scores.DT = scores.DT))
    #return(list(scores.DT = data.table(a=c(1,2,3),b=c(2,3,4))))
  })
  
  output$varsel_anova <- renderText({
    anova()
  })
  
  output$varsel_scores <- renderTable({
    #anova()
    if (is.null(input$filemanager_folder)) return(NULL)
    return(build_scores()$scores.DT)
  })
  
  output$varsel_codebook <- renderDataTable({
    if (is.null(input$filemanager_folder)) return(NULL)
    if (!file.exists(file.path(path2data,input$filemanager_folder,"codebook.csv"))) return(NULL)
    codebook.dt <- fread(file.path(path2data,input$filemanager_folder,"codebook.csv"))
    setnames(
      codebook.dt,
      c("mbd","variable","value","res_freq","res_proportion","panel_freq","panel_proportion","valid")
      , c("MBD","Secondary Characteristic","Secondary Characteristic Value","Count (RES)", "Proportion (RES)", "Count (Panel)", "Proportion (Panel)", "Valid flag")
    )
    return(codebook.dt)
  })
  
  anova <- eventReactive(input$varsel_run_button,{
    path <- file.path(path2data,input$filemanager_folder)
    #path <- file.path(path2data,"COUNTRY - INDEX - PERIOD")
    files <- list.files(path)
    
    codebook.DT <- fread(file.path(path,"codebook.csv"))
    codebook.DT[,score := 0]
    sales.DT <-  fread(
      file.path(path,files[grepl("_sales.csv",files)][1])
      , colClasses = c(outlet = "char")
    )
    sales.DT[,period := NULL]
    source("proc_prepare_sales.R", local = TRUE)
    
    # save the modified sales into the folder
    
    write.csv(
      sales.DT
      , file.path(path,"sales_clean.csv")
      , quote = FALSE
      , row.names = FALSE
    )
    
    chars.DT <- fread(file.path(path,files[grepl("_panel.csv",files)][1]), colClasses = c(outlet = "char"))
    
    withProgress(message = 'Running anova...', value = 0, {
      source("proc_anova_step.R", local = TRUE) # source("proc_anova_step.R")
    })
    
    write.csv(
      codebook.DT
      , file.path(path,"codebook.csv")
      , quote = FALSE
      , row.names = FALSE
    )
  })
  
  #######################################################################################
  #######################################################################################
  ### sample check tab
  
  output$samplecheck_cur_folder <- renderText(return(input$filemanager_folder))
  
  abutton_1 <- eventReactive(input$goButton_1, {
    if(is.null(input$filemanager_folder)) return(NULL)
    
    return("ok")
  })
  
  
  data <- reactive({
    
    if (is.null(abutton_1())) return(NULL)
    isolate({
      path <- file.path(path2data,input$filemanager_folder) # path <- file.path(path2data, "MYM_BASE_201504")
      files <- list.files(path)
      
      sf <- fread(file.path(path,files[grepl("_universe.csv",files)][1]), colClasses = c(outlet = "char"))
      pf <- fread(file.path(path,files[grepl("_panel.csv",files)][1]), colClasses = c(outlet = "char"))
      sd <- fread(file.path(path,files[grepl("_sample_design.csv",files)][1]))
      codebook <- fread(file.path(path,"codebook.csv"))
      
      return(list(sf = sf, pf = pf, sd = sd, codebook=codebook))
    })
  })
  
  
  process <- reactive({
    
    if(is.null(data())) return(NULL)
    
    mbd_rotation <- NULL
    pvalue <- input$pvalue
    min_retain_ss <- input$min_retain_ss
    sf0 <- data()$sf
    pf0 <- data()$pf
    sd <- data()$sd
    codebook <- data()$codebook
    
    withProgress(message = 'Running MBD Importance...', value = 0, {
      source("SubcellAllocation_source.R", local = TRUE)
    })
    
    # source("SubcellAllocation_source.R", local = TRUE)
    list(subcell_all = subcell_all,
         cell_all = cell_all,
         subcell2_all = subcell2_all,
         cell2_all = cell2_all,
         
         output_mbdreport = output_mbdreport,
         output_kpireport_all = output_kpireport_all,
         
         output_chisq1_all = output_chisq1_all,
         output_subcell1_all = output_subcell1_all,
         output_chisq2_all = output_chisq2_all,
         output_subcell2_all = output_subcell2_all)
  })
  
  
  output$status<- renderText({
    if (!is.null(data())&!is.null(process())) { return("TRUE")}
  })
  outputOptions(output, 'status', suspendWhenHidden = FALSE)
  
  output$choose_sec_char <- renderUI({
    if (is.null(data())&is.null(process())) return(NULL)
    sf <- data()$sf
    all_chars <- as.vector(names(sf)[grepl("sec_var_",names(sf))])
    
    selectInput("sec_char", "Choose Secondary Characteristic", as.list(all_chars), selected = NULL, multiple = F)
  })
  
  sec_char <- reactive({
    
    if (is.null(input$sec_char)) return(NULL)
    sec_char <- input$sec_char
  })
  
  output$choose_sample_rotation_mbd <- renderUI({
    
    if (is.null(data())&is.null(process())) return(NULL)
    sf <- data()$sf
    mbd_unique <- as.vector(unique(sf$mbd))
    
    selectInput("mbd_rotation", "Choose MBD (Multiple MBD can be selected)", as.list(mbd_unique), selected = NULL, multiple = T)
  })
  
  mbd_rotation <- reactive({
    
    if (is.null(input$mbd_rotation))  return(NULL)
    mbd_rotation <- input$mbd_rotation
  })
  
  abutton <- eventReactive(input$goButton, {
    sec_char <- sec_char()
    mbd_rotation <- mbd_rotation()
    list(sec_char = sec_char, mbd_rotation = mbd_rotation)
  })
  
  
  
  process2 <- reactive({
    if (is.null(data())&is.null(process())&is.null(abutton()$mbd_rotation)) return(NULL)
    
    min_retain_ss <- input$min_retain_ss
    pvalue <- input$pvalue
    sec_char <- abutton()$sec_char
    mbd_rotation <- abutton()$mbd_rotation
    
    subcell_all = process()$subcell_all
    cell_all = process()$cell_all
    subcell2_all = process()$subcell2_all
    cell2_all = process()$cell2_all
    
    # output_mbdreport = process()$output_mbdreport
    output_kpireport_all = process()$output_kpireport_all
    
    output_chisq1_all = process()$output_chisq1_all
    output_subcell1_all = process()$output_subcell1_all
    output_chisq2_all = process()$output_chisq2_all
    output_subcell2_all = process()$output_subcell2_all
    
    withProgress(message = 'Running Micro Rep...', value = 0, {
      source("SubcellAllocation_source2.R", local = TRUE)
    })
    
    list(output_chisq1 = output_chisq1
         , output_subcell1 = output_subcell1
         , output_chisq2 = output_chisq2
         , output_subcell2 = output_subcell2
         , output_kpireport = output_kpireport
         , output_log = output_log
    )
  })
  
  output$output_mbdreport <- renderDataTable(
    if(is.null(process()$output_mbdreport)){
      return()} else {
        process()$output_mbdreport[, list(`Secondary Characteristic`,`MBD`,`Variable Selection Score`, 
                                          `Max. MBD Importance Score`,`Min. Chi-square p-value`, `Sample Rotation Recommendation`)]
      }
    , options = list(pageLength = 10, autoWidth = TRUE)
  )
  
  output$remark <- renderUI({
    pvalue <- input$pvalue
    if(is.null(process()$output_mbdreport)){
      return(NULL)}else{
        HTML(paste(strong('Sample Rotation Recommendation: '),
                   paste0('Variable Selection Score > 1.5, Maximum MBD Importance Score > 0 and Minimum p-value of Chi-square test within MBD < ',pvalue),
                   sep = '<br/>'))
      }
  })
  
  output$output_chisq1 <- renderDataTable(
    if(is.null(process2()$output_chisq1)){
      return()}else{
        pvalue <- input$pvalue
        source("SubcellAllocation_source_functions.R")
        output_chisq1 <- process2()$output_chisq1[, list(`MBD`, `Cell`, `Number of Sub-cells`, 
                                                         `Count SF`, `Target SS`, `Current SS`, `Min. Retain SS`, `Final Subcell Allocation`
                                                         , `Additional SS`, `Drop SS`, `P-value (current)`, `P-value (after NSC)`)]
        output_chisq1 <- output_chisq_func(output_chisq1,pvalue)
        output_chisq1
      }
    , options = list(pageLength = 10)
  )
  
  output$output_subcell1 <- renderDataTable(
    process2()$output_subcell1[, list(`MBD`, `Cell`, `Secondary Characteristic Value`
                                      , ` Count SF`, `Current SS`, `Final Sub-cell Allocation`
                                      , `Additional SS`,`Drop SS`, `Sub Cell Priority Rank`)]
    , options = list(pageLength = 10)
  )
  
  output$output_kpireport <- renderDataTable(
    process2()$output_kpireport[, list(`MBD`, `Secondary Characteristic Value`, `Panel Size`,`Panel Estimate`,
                                       `Universe`, `Z Score`, `Sample Completeness`,`Universe Completeness`,
                                       `Ratio 1`,`Ratio 2`,`Concern`,`MBD Importance Score`)]
    , options = list(pageLength = 10)
  )
  
  output$note_sec1 <- renderUI({
    sec_char_select <- sec_char()
    if(is.null(process2()$output_kpireport)){
      return(NULL)}else{
        HTML(paste(strong("Secondary Characteristic Selected:  ", span(sec_char_select, style = "color:red"))))
      }
  })
  
  output$note_sec2 <- renderUI({
    sec_char_select <- sec_char()
    if(is.null(process2()$output_kpireport)){
      return(NULL)}else{
        HTML(paste(strong("Secondary Characteristic Selected:  ", span(sec_char_select, style = "color:red"))))
      }
  })
  
  output$note_sec3 <- renderUI({
    sec_char_select <- sec_char()
    if(is.null(process2()$output_kpireport)){
      return(NULL)}else{
        HTML(paste(strong("Secondary Characteristic Selected:  ", span(sec_char_select, style = "color:red"))))
      }
  })
  
  output$note_sec4 <- renderUI({
    sec_char_select <- sec_char()
    if(is.null(process2()$output_kpireport)){
      return(NULL)}else{
        HTML(paste(strong("Secondary Characteristic Selected:  ", span(sec_char_select, style = "color:red"))))
      }
  })
  
  output$note_sec5 <- renderUI({
    sec_char_select <- sec_char()
    if(is.null(process2()$output_kpireport)){
      return(NULL)}else{
        HTML(paste(strong("Secondary Characteristic Selected:  ", span(sec_char_select, style = "color:red"))))
      }
  })
  
  output$remark2 <- renderUI({
    pvalue <- input$pvalue
    if(is.null(process2()$output_kpireport)){
      return(NULL)}else{
        HTML(paste( strong('Concern: '),
                    p(em('"A"'), 'will refer to our uncertainty, whether the incomplete RES &panel observations for the specific mbd* secchar combination could be the source of the high structural observed difference;'),
                    p(em('"B"'), 'indicates that panel RES was incomplete, both in number of stores and numerical importance of the missing observations within the mbd;'),
                    p(em('"C"'), 'is less serious condition, as the missing panel observations have less than 5% numerical weight in the respective universe for the secondary characreristic;'),
                    p(em('"D"'), 'is an indication that universe instruction process is broken; we shall have less confidence in the outcome;'),
                    p(em('"E"'), 'indicates an even bigger discrepancy between the universe figures, trigger for input correction;'),
                    p(em('"F"'), 'indicates that the calculation outcome has higher than 10% sampling error, i.e. the conclusions are based on relatively low sample sizes;'),
                    p(em('"G"'), 'indicates that the calculation outcome has higher than 15% sampling error, i.e. the action should get low priority;'),
                    p(em('"H"'), 'indicates that the given characteristic has relatively low importance within the reported universe (less than 5%),therefore the action should get low priority;'),
                    p(em('"I"'), 'indicates that the importance of the group is still relatively low, the "punishment" will be slightly less, but still indicating the lower importance;'),
                    p(em('"J"'), 'indicates that the relative difference betweeen the "true" and represented mbd*secondary characreristic value universe is beyond 20%, 
                      therefore the action should get higher priority, if all other scores are the same.'),
                    p(strong('Completeness: ')),
                    p(em('Sample Completness'), 'shows the % of available observations in the sample'), 
                    p(em('Universe Completeness'), 'shows what part of the universe is represented by these observations'),
                    p(strong('Ratio: ')),
                    p(em('Ratio 1'), 'is the proportion of all stores in the universe vs universe of those stores where the observation was available'),
                    p(em('Ratio 2'), 'is the ratio of the two estimates: panel based divided by RES based estimation.')))
        # sep = '<br/>'))
      }
  })
  
  output$note <- renderUI({ 
    
    source("SubcellAllocation_source_functions.R")
    
    pvalue <- input$pvalue
    mbdreport <- process()$output_mbdreport
    mbd_rotation <- abutton()$mbd_rotation
    output_chisq1 <- process2()$output_chisq1
    output_chisq2 <- process2()$output_chisq2
    text2 <- NULL
    
    note_text <- sample_rotation_note(pvalue, mbdreport, mbd_rotation, output_chisq1, output_chisq2, text2)
    text1 <- note_text[[1]]
    text2 <- note_text[[2]]
    
    HTML(paste(strong(span(text1, style = "color:blue")), strong(span(text2, style = "color:blue")), sep = '<br/>'))
  })
  
  
  output$output_chisq2 <- renderDataTable({
    
    source("SubcellAllocation_source_functions.R")
    
    pvalue <- input$pvalue
    if(is.null(process2()$output_chisq2)){
      return()} else {
        output_chisq2 <- as.data.table(process2()$output_chisq2)[, list(`MBD`,`Cell`,`Number of Sub-cells`,`Count SF`
                                                                        ,`Target SS`,`Min. Retain SS`,`Current SS`,`Final Subcell Allocation`,
                                                                        `Additional SS`, `Drop SS`,`P-value (current)`,`P-value (after NSC)`, `P-value (after rotation)`)]
        output_chisq2 <- output_chisq_func(output_chisq2,pvalue)
        return(output_chisq2)
      }
  }
  , options = list(pageLength = 10)
  )
  
  output$output_subcell2 <- renderDataTable(
    
    if(is.null(process2()$output_subcell2)){
      return()} else{
        return(as.data.table(process2()$output_subcell2)[, list(`MBD`,`Cell`,`Secondary Characteristic Value`,`Count SF`,
                                                                `Current SS`,`Final Subcell Allocation`
                                                                ,`Additional SS`,`Drop SS`,`Sub Cell Priority Rank`)])
      }
    ,  options = list(pageLength = 10)
  )
  
  output$downloadData <- downloadHandler(
    
    filename = function() { paste0("MRA Report - ",input$sec_char," - ",format(Sys.time(), format = "%Y%m%d-%H%M%S"),".xlsx") },
    
    content = function(file){
      source("SubcellAllocation_source_functions.R")
      fname <- paste(file,"xlsx",sep=".")
      wb <- loadWorkbook(fname, create = TRUE)
      
      pvalue_style_list <- pvalue_style(wb)
      pvalue1 <- pvalue_style_list[[1]]
      pvalue2 <- pvalue_style_list[[2]]
      pvalue3 <- pvalue_style_list[[3]]
      
      ### Log sheet
      createSheet(wb, name = "Log")
      writeWorksheet(wb, process2()$output_log, sheet = "Log")
      
      ### MBD Importance Report
      createSheet(wb, name = "MBD Importance Summary")
      writeWorksheet(wb, process()$output_mbdreport, sheet = "MBD Importance Summary")
      
      ### KPI Report
      createSheet(wb, name = "MBD Importance Analysis")
      writeWorksheet(wb, process2()$output_kpireport, sheet = "MBD Importance Analysis")
      
      ### Output without Sample Rotation
      createSheet(wb, name = "NSC - Chi-Sqr Test")
      writeWorksheet(wb, process2()$output_chisq1, sheet = "NSC - Chi-Sqr Test")
      createSheet(wb, name = "NSC - Subcell Allocation")
      writeWorksheet(wb, process2()$output_subcell1, sheet = "NSC - Subcell Allocation") 
      setColumnWidth(wb, sheet = "NSC - Subcell Allocation", column = 1, width = 6000)
      
      colIndex <- which(names(process2()$output_chisq1) %in%c("P-value (current)", "P-value (after NSC)"))
      highlight_func(wb,process2()$output_chisq1,colIndex[1],0.1,pvalue3,"NSC - Chi-Sqr Test")
      highlight_func(wb,process2()$output_chisq1,colIndex[1],0.05,pvalue2,"NSC - Chi-Sqr Test")
      highlight_func(wb,process2()$output_chisq1,colIndex[1],0.01,pvalue1,"NSC - Chi-Sqr Test")
      
      highlight_func(wb,process2()$output_chisq1,colIndex[2],0.1,pvalue3,"NSC - Chi-Sqr Test")
      highlight_func(wb,process2()$output_chisq1,colIndex[2],0.05,pvalue2,"NSC - Chi-Sqr Test")
      highlight_func(wb,process2()$output_chisq1,colIndex[2],0.01,pvalue1,"NSC - Chi-Sqr Test")
      
      
      
      if(!is.null(process2()$output_chisq2)){
        
        createSheet(wb, name = "Micro Rep - Chi Sqr Test")
        writeWorksheet(wb, process2()$output_chisq2, sheet = "Micro Rep - Chi Sqr Test")
        createSheet(wb, name = "Micro Rep - Subcell Allocation")
        writeWorksheet(wb, process2()$output_subcell2, sheet = "Micro Rep - Subcell Allocation") 
        setColumnWidth(wb, sheet = "Micro Rep - Subcell Allocation", column = 1, width = 6000)
        
        colIndex <- which(names(process2()$output_chisq2) %in%c("P-value (current)", "P-value (after NSC)", "P-value (after rotation)"))
        highlight_func(wb,process2()$output_chisq2,colIndex[1],0.1,pvalue3,"Micro Rep - Chi Sqr Test")
        highlight_func(wb,process2()$output_chisq2,colIndex[1],0.05,pvalue2,"Micro Rep - Chi Sqr Test")
        highlight_func(wb,process2()$output_chisq2,colIndex[1],0.01,pvalue1,"Micro Rep - Chi Sqr Test")
        
        highlight_func(wb,process2()$output_chisq2,colIndex[2],0.1,pvalue3,"Micro Rep - Chi Sqr Test")
        highlight_func(wb,process2()$output_chisq2,colIndex[2],0.05,pvalue2,"Micro Rep - Chi Sqr Test")
        highlight_func(wb,process2()$output_chisq2,colIndex[2],0.01,pvalue1,"Micro Rep - Chi Sqr Test")
        
        highlight_func(wb,process2()$output_chisq2,colIndex[3],0.1,pvalue3,"Micro Rep - Chi Sqr Test")
        highlight_func(wb,process2()$output_chisq2,colIndex[3],0.05,pvalue2,"Micro Rep - Chi Sqr Test")
        highlight_func(wb,process2()$output_chisq2,colIndex[3],0.01,pvalue1,"Micro Rep - Chi Sqr Test")
      }
      saveWorkbook(wb)
      file.rename(fname,file)
    }
  )
  
  
}


shinyApp(ui = ui, server = server)
