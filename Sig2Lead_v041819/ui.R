#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(visNetwork)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
  
  # Application title
  titlePanel("Sig2Lead"),
      fluidPage(
        tabsetPanel(
          tabPanel("Search", value ="search", icon = NULL,
                   # fluidRow(
############# Add this ##############################################     
                   fluidRow(
                     column(3, 
                            wellPanel(
                              
#####################################################################                   
                     br(),
                            selectInput("Signature", NULL, list("Input a Gene", "Upload a Signature"), selected = "Input a Gene"),
                     br(),
                     
                     conditionalPanel(
                       condition = "input.Signature == 'Input a Gene'",
                            textInput("gene_knockdown", "Input a Gene", value="bcl2a1")
                     ),
                     conditionalPanel(
                        condition = "input.Signature == 'Upload a Signature'",
                            fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".csv"), width = NULL, 
                                         buttonLabel = "Browse...", placeholder = "No file selected")
                     ),
                            br(),
                     fileInput("AddedCompounds", "Add compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
                               buttonLabel = "Browse...", placeholder = "No file selected"),
                   uiOutput("mytest"),
                    # conditionalPanel(
                    #      condition = "input.AddedCompounds.name == null",
                    #         textInput("AddedLabel", "Label for added compounds", value="Added"),
                    #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
                    # ),
                            br(),
                   selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
                   br(),
                   #conditionalPanel(
                  #     condition = "if(!is.null(input.AddedCompounds))",
                  #     checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
                  #),     
                   #br(),
                            actionButton(label="Go!", "Go")
############# Add this ######################################################################################                  
                            )),
#fluidRow(
column(9, 
       htmlOutput("notify1"),
       visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
       htmlOutput("notify2")
       
))


# )

#############################################################################################################

           # )
        ),
        
#        mainPanel(
#        textOutput("ItWorked"),
          tabPanel("LINCS Compounds", value = "LINCS Compounds", icon =NULL,
          DT::dataTableOutput("SMILES"),
          downloadButton("SMILESDownload", label = "Download SMILES"),
          downloadButton("SDFDownload", label = "Download SDF"),
          downloadButton("ConTableDownload", label="Download Table")
        ),
      tabPanel("Heatmap", value="heatmap", icon = NULL, 
#               actionButton(label="GenerateHeatmap", "GenerateHeatmap"),
               textInput("CutHeight", "Tanimoto Similarity", value="0.75"),
               textInput("ClusterSize", "Minimum Cluster Size", value="3"),
               actionButton(label="Get MDS Plot", "GetRepresentatives"),
               plotOutput("distPlot", width = "850px", height = "800px")
               
              
        ),
      tabPanel("Combined Score", value="Combined", icon=NULL,
               DT::dataTableOutput("Combined_Score"),
               downloadButton("Max_Scores", label = "Download Scores")),
      #        #downloadButton("NCIDownload", label = "Download Similar NCI"))


      tabPanel("MDS Plot", value="MDS", icon=NULL,
               plotOutput("MDSPlot", width = "900px", height = "800px"),
               DT::dataTableOutput("Representatives"),
               downloadButton("RepDownload", label="Download Representatives"),
               downloadButton("ClusterDownload", label = "Download Clusters"),
               actionButton(label="Get Related NCI Compounds", "GetNCI")
               ),


      tabPanel("Similar Compounds", value="NCI", icon=NULL,
              DT::dataTableOutput("Similar_NCI"),
              downloadButton("NCIDownload", label = "Download Similar NCI")),
#### add this script to UI #### 

tabPanel("Global STITCH Network", value = "mSTITCH", icon = NULL,
         fluidRow(
           column(3, wellPanel(
             actionButton(label="View Global STITCH Network", "GlobalSTITCH"),
             checkboxInput("all_clusters", "Shows all clusters", value = FALSE, width = '400px')
           )),
           column(9, 
                  br(),      
                  textOutput("numberClust")    
           )
         ),
         fluidRow(
           column(12, 
                  visNetworkOutput("mSTITCHPlot", width = "1000px", height = "750px")
           )
         )),

tabPanel("Selected Cluster Network", value = "STITCH", icon = NULL,
         fluidRow(
           column(2, wellPanel(
             textInput("ClusterNumber", "Cluster Number", value = 1),
             br(),
             selectInput("confidence", "Confidence Score", choices = c(150, 400, 700, 950), selected = 400),
             br(),
             sliderInput("Connections", "Numbers of Connections", min = 2, max = 25, value = 10, ticks = FALSE),
             br(),
             # textInput("Gene", "Gene of Interest", value = "bcl2a1"),
             checkboxInput("fixedNet", "Fixed node position", value = FALSE, width = '400px'),
             checkboxInput("show_undetectNodes", "Shows all", value = FALSE, width = '400px'),
             actionButton(label="View with STITCH", "GetSTITCH")
           )),
           column(9, 
                  br(),
                  visNetworkOutput("STITCHPlot", width = "900px", height = "600px")
           )
         ))


############################
        )
      )
  )
)


