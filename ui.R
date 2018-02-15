library(shiny)
library(shinyjs)
library(prospectr)
library(shinyFiles)
library(DT)
library(httr)
#library(caret)
#library(jsonlite)
library(xlsx)
library(soil.spec)
library(reshape2)
library(ggplot2)
library(dplyr)
library(plotly)

library(xlsx)

library(shinyjs)

appCSS <- "#loading-content {  position: absolute;  background: #000000;  opacity: 0.9;  z-index: 100;  left: 0;  right: 0;  height: 100%; text-align: center; color: #FFFFFF;}"



ui<-fluidPage(	
  useShinyjs(),
    inlineCSS(appCSS),


      # Loading message
div(
  id = "loading-content",
    h2("Loading...")
      ),
        theme = shinytheme("cerulean"),titlePanel(fluidRow(column(8, col = "blue","SpecWeb: App for Spectroscopy Data Processing and Analysis",offset = 1),
           column(1, img(height = 80, width  = 200, src = "logo.png")))),



navbarPage(title = "Infrared Spectral Processing Methods",inverse = FALSE,
  tabPanel("Read OPUS files",
    sidebarLayout(
      sidebarPanel(
        shinyDirButton('diropus', 'Select OPUS folder', 'Please select a folder'),
        br(),
        br(),
        br(),
        br(),
          shinySaveButton('saveop', 'Save file', 'Save file as ...,', filetype = list (text="csv","txt")),    
            tags$p(),
              p("OPUS (OPtical User Software) is the Bruker data collection and analysis program for Alpha, Multi-Purpose   Analyzer (MPA) and Tensor 27 FT-IR (HTS-xt) spectrometers. Spectral data files stored in this format are characterized by a numeric extension, usually a zero unless there are duplicates which are given increamented by one. There are lots of spectral details  stored in form of data-blocks  inside these files ranging from date and time of measurement, type of instrument used, Absorbance values, IR regions and much, please see OPUS_ver.xx BasePackage manual."),width = 4),
                mainPanel(
                  br(),
                    p("This tab contain a module for reading infrared spectral files of the format OPUS files recorded using three Bruker spectrometers into R to be used in chemometric methods e.g. developing calibration models, for prediction, clustering, etc. CO2 bands (2380-2351 cm-1) within the MIR data should be removed during calibration model fitting from the original OPUS format to obtain unbiased comparable models from different parts of the world with varying CO2 concentration in the atmosphere."),
                      h4("A cross-section of converted OPUS files to be shown here..."),
                        dataTableOutput('fileopus')
                          )
                            )
                              ),



navbarMenu("Exploratory Analysis",  
  tabPanel("View converted spectra",
    sidebarLayout(
      sidebarPanel(
        fileInput('file0', 'Choose CSV File',accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))
        ),
          mainPanel(
            h3("Visualize a few raw spectra"),
       	      tableOutput('view.raw'),
                h3("Raw spectral signatures"),
       	          plotOutput('raw.plot')
		                )
                      )
                        ),

      
  tabPanel("Preprocess raw spectra",
    sidebarLayout(
      sidebarPanel(
        fileInput('preprospec', 'Choose CSV File',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
          uiOutput("preprospecsel")
            ),
              mainPanel(
                div(DT::dataTableOutput("preprospectable"), style = "font-size: 100%; width: 92%")
                    )
                      )
                        ),



    tabPanel("Principal Component Analysis",
      sidebarLayout(
        sidebarPanel(
          fileInput('file1', 'Choose CSV File', accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
           tags$hr(),
            tags$hr(),
              textInput("fileNamer", "Raw File Name"),
                shinySaveButton('saveraw', 'Save file', 'Save file as ...,', filetype = c("png","pdf")),   
                  tags$hr(),
                    textInput("fileNamep", "Preprocessed File Name"),
                      shinySaveButton('saveprocessed', 'Save file', 'Save file as ...,', filetype = c("png","pdf"))  
                        ),
                          mainPanel(
                            h3("Visualize a few raw spectra"),
                              tableOutput('pca.raw'),
                                h3("Raw spectral signatures"),
                                  plotOutput('pca.plot1'),
                                    h3("Preprocessed spectral signatures"),
                                      plotOutput('pca.plot2'),
                                        h3("PCA scores"),
                                          plotOutput('pca.plot3', click = "plot_click"),
                                            verbatimTextOutput("clicked"),
                                              h3("Visualize list of selected outliers"),
                                                tableOutput('outliers')
                                                  )
                                                    )
                                                      ),
  tabPanel("Outliers")

    ),
      

navbarMenu("Data Quality Control",  
  tabPanel("Internal lab stds",
    sidebarLayout(
      sidebarPanel(
        selectInput('alphadt',"Select standard",choices=c("Mua","Whitesand")),
          br(),
            shinyFilesButton('file', "Select directory with std files", "Please select directory with std files",multiple=TRUE),
              textOutput('txt'),
                br(),
                  h5("Uploaded new file is:"),
                    verbatimTextOutput('table1')),
                      mainPanel(
                        plotOutput('contents')
                          )
                            )
                              ),
      
  tabPanel("View replicate spectra")
    ),



navbarMenu("Tools",  
  tabPanel("QR-codes",
    sidebarPanel(
      numericInput('n', 'Number of Qr codes required:',min = 1, max = 100000, n),
        textInput('study', "Study/Project Prefix:", study),
          shinySaveButton("saveqr", "Save QR codes", "Save file as ...", filetype=list(csv="csv",xlsx="xlsx"))  
            ),
              mainPanel(
                dataTableOutput("qrcodes")
                  )
                    ),
      
  tabPanel("Kennard_Stone",
    sidebarLayout(
      sidebarPanel(
        fileInput('filek', 'Choose CSV file with infrared data',accept=c('text/csv','text/comma-separated-values,text/plain','.csv'), width="900px"),
          checkboxInput('header','Header',TRUE),radioButtons('sep','Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
            radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"'),
              helpText("Note: File read should have spectra starting from second column"),
                sliderInput(inputId="perc", label=h3("Proportion for the reference samples",align= "center"), min=0.01, max=0.99, value=0.1, step=0.01),
                  helpText("Note: value should be between more than 0 and less than 1."),
                    tags$hr(),
                      radioButtons("filetype", "file type:", choices = c("csv","tsv")),
                        shinySaveButton("saveks", "Save", "Save As ...", filetype=list(csv="csv",xlsx="xlsx",xls="xls"))
                          ),
                            mainPanel(
                              h3("Visualize a few raw spectra"),
                                tableOutput('contents2'),
                                  h3("Raw spectral signatures"),
                                    plotOutput('plot1'),
                                      h3("Sample selection by Kennard_Stone"),
                                        plotOutput('plot2'),
                                          h3("Visualize list of selected samples"),
                                            tableOutput('selected')
                                              )
                                                )
                                                  ),
  tabPanel("Fertilizer_Screening")
  

    ),


tabPanel("Calibration",sidebarLayout(
  sidebarPanel(
    fileInput('ircalfile', 'Choose IR CSV File',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
      tags$hr(),
        tags$hr(),
          fileInput('refcalfile', 'Choose Ref CSV File',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
            tags$hr(),
              tags$hr(),
                shinyDirButton("dir", "Chose working directory", "working directory"),
                  tags$hr(),
                    radioButtons("model","Choose regression method for calibration",choices = c("PLS","RF"), width = 400,selected = "PLS"),
                      tags$hr(),
                        actionButton("runaction", "Calibrate")
                          ),
                            mainPanel(
                              dataTableOutput('ircalfiletable'),
                                dataTableOutput('refcalfiletable')
                                  )
                                    )
                                      ),

             
tabPanel("Infrared_Predictions",sidebarLayout(
  sidebarPanel(
    fileInput('predfile', 'Choose New CSV file to predict',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
      tags$hr(),
        shinyDirButton("dirmodels", "Choose folder with models", "Models folder"),
          actionButton("comparedata", "Compare"),
            tags$hr(), 
              actionButton("predaction", "Predict"),
	             tags$hr(),
                uiOutput("make_box"),            
                  tags$hr(),
                    dataTableOutput("mytable")
                      , width = 3), 
                        mainPanel(
                          verbatimTextOutput("dirpaths"),
                            plotOutput('predfile'),
                              div(DT::dataTableOutput("predfiletable"), style = "font-size: 100%; width: 92%")
                                )
                                  )
                                    )
                                      )
                                        )