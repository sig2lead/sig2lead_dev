#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
# #########################################################
# # Install Bioconductor Packages
# ###############################
# 
# list.of.packages_bc <- c("ChemmineOB", "ChemmineR")
# new.packages_bc <- list.of.packages_bc[!(list.of.packages_bc %in% installed.packages()[,"Package"])]
# 
# if(length(new.packages_bc)){ 
#   
#   if (!requireNamespace("BiocManager", quietly=TRUE)){
#     BiocManager::install(version="3.11")
#     BiocManager::install(new.packages_bc, ask=FALSE)}
# }

#######################################################################
# Install all required packages
###############################

list.of.packages <- c(
  #"shiny",
  "shinyjs",
  "ggplot2",
  "plotly",
  "httr",
  "jsonlite",
  "DT",
  "visNetwork",
  #"ChemmineOB",
  #"ChemmineR",
  "plyr",
  "dendextend",
  "colorspace",
  "ggforce",
  "rlist",
  "scatterpie",
  "ggrepel",
  "bazar",
  "XML",
  "RCurl",
  "bitops",
  "scrapeR",
  "igraph",
  "circlize",
  "enrichR",
  "readr",
  "dplyr",
  "gplots"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE, INSTALL_opts = '--no-lock')

############## Bioconductor ################################################################
# 

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version="3.12")

list.of.packages_bm <- "BiocManager"
new.packages_bm <- list.of.packages_bm[!(list.of.packages_bm %in% installed.packages()[,"Package"])]
if(length(new.packages_bm)) install.packages(new.packages_bm, dependencies=TRUE, INSTALL_opts = '--no-lock', version="3.12")


list.of.packages_bc <- c("ChemmineOB", "ChemmineR")
new.packages_bc <- list.of.packages_bc[!(list.of.packages_bc %in% installed.packages()[,"Package"])]

if(length(new.packages_bc)){
  
  BiocManager::install(new.packages_bc, ask=FALSE)
}

##########################################################
# Load Libraries
#################
library(shiny)
library(shinyjs)
library(plotly)
library(visNetwork)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
  
  # Application title
  titlePanel("Sig2Lead"),
      fluidPage(useShinyjs(),
        tabsetPanel(id = "tabs",
          tabPanel("Signature Connectivity Analysis", value ="search", icon = NULL,
                   # fluidRow(
############# Add this ##############################################     
                   
                     #fluidRow(
                     #column(3, 
                    #        wellPanel(
                     sidebarPanel(position="left",        
#####################################################################                   
                     #br(),
                            selectInput("Signature", NULL, list("Define Target Gene", "Upload a Signature", "Find Analogs in LINCS"), selected = "Input a Gene"),
                     br(),
                     
                     conditionalPanel(
                       condition = "input.Signature == 'Define Target Gene'",
                            textInput("gene_knockdown", "Input a Target Gene", value="bcl2a1")
                     ),
                     conditionalPanel(
                        condition = "input.Signature == 'Upload a Signature'",
                            fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".txt", ".csv"), width = NULL, 
                                         buttonLabel = "Browse...", placeholder = "No file selected")
                     ),
                            br(),
                     fileInput("AddedCompounds", "Add candidate compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
                               buttonLabel = "Browse...", placeholder = "No file selected"),
                   uiOutput("mytest"),
                    # conditionalPanel(
                    #      condition = "input.AddedCompounds.name == null",
                    #         textInput("AddedLabel", "Label for added compounds", value="Added"),
                    #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
                    # ),
                            br(),
                   #selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
                   br(),
                  actionButton(label="Show Advanced Options", "Advanced_button"),
                  br(),
                  br(),

                  
                   #checkboxInput("Advanced", "Advanced Options", selected),
                   conditionalPanel(condition="input.Advanced_button ==  1",
                        selectInput("ConOrDiscon","Activation or Inhibition", list("Inhibit", "Activate"), selected = "Inhibit"),            
                        conditionalPanel(condition="input.ConOrDiscon == 'Inhibit'",
                          selectInput("Algorithm","Chemical Similarity Search", list("minSim", "fpSim"), selected = "minSim"),
                          numericInput("Concordance", "Concordance Threshold", value=0.2, min=0.2, max=0.5, step=0.1)),
                        conditionalPanel(condition="input.ConOrDiscon == 'Activate'",
                                         selectInput("Algorithm","Chemical Similarity Search", list("minSim", "fpSim"), selected = "minSim"),
                                         numericInput("Concordance", "Discordance Threshold", value=-0.2, min=-0.5, max=-0.2, step=-0.1))
                        
                   ),  
                  br(),
                  # conditionalPanel(
                  #     condition = "if(!is.null(input.AddedCompounds))",
                  #     checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
                  #),     
                   #br(),
                            actionButton(label="Go!", "Go"),

                  br(),
                  br(),
                  br(),
                  # output$downloadCandidates <- renderUI({
                  #   req("input.Go==1")
                  #   downloadButton("CandidatesDownload", label = "My Candidates Ranked")
                  #   }),
                  uiOutput("downloadCandidates"),
                  #br(),
                  uiOutput("downloadLINCS")
                  #downloadButton("LINCSDownload", label="LINCS Candidates Ranked"),
                  
                  #downloadButton("SMILESDownload", label = "Download SMILES"),
                  #downloadButton("SDFDownload", label = "Download SDF"),
                  #downloadButton("ConTableDownload", label="Download Table"),
                  #downloadButton("NCIDownload", label = "Download Similar NCI")
############# Add this ######################################################################################                  
                ),
#fluidRow(
#column(9, 
#       htmlOutput("notify1"),
#       visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
#       htmlOutput("notify2")
       
#)),


# )

#############################################################################################################
      mainPanel(position="right",
                
          uiOutput("check_kd"),  
          br(),
          br(),
          DT::dataTableOutput("CANDIDATES"),
          #downloadButton("CandidatesDownload", label = "Download Candidates"),
          #downloadButton("SMILESDownload", label = "Download SMILES"),
          #downloadButton("SDFDownload", label = "Download SDF"),
          br(),
          br(),
          DT::dataTableOutput("SMILES"),
          br(),
          br(),
          #downloadButton("ConTableDownload", label="Download Table"),
          br(),
          br(),
          DT::dataTableOutput("Similar_NCI")
          #downloadButton("NCIDownload", label = "Download Similar NCI")
          
          
         
      )     # )
    ), #up to here
        
#        mainPanel(
#        textOutput("ItWorked"),
      #tabPanel("Candidate Compounds", value = "Candidate Compounds", icon =NULL,
               
       #  DT::dataTableOutput("CANDIDATES"),
         #downloadButton("SMILESDownload", label = "Download SMILES"),
         #downloadButton("SDFDownload", label = "Download SDF"),
         #downloadButton("ConTableDownload", label="Download Table")
      #),
      #tabPanel("LINCS Compounds", value = "LINCS Compounds", icon =NULL,
      #    DT::dataTableOutput("SMILES")
      #    downloadButton("SMILESDownload", label = "Download SMILES"),
      #    downloadButton("SDFDownload", label = "Download SDF"),
      #    downloadButton("ConTableDownload", label="Download Table")
      #  ),
#tabPanel("Consensus Ranking", value="Combined", icon=NULL,
         #          fluidRow(
         #            column(3, 
         #                   wellPanel(
         #                     
         #                     #####################################################################                   
         #                     br(),
         #                     selectInput("Signature", NULL, list("Input a Gene", "Upload a Signature", "Similarity Search"), selected = "Input a Gene"),
         #                     br(),
         #                     
         #                     conditionalPanel(
         #                       condition = "input.Signature == 'Input a Gene'",
         #                       textInput("gene_knockdown", "Input a Gene", value="bcl2a1")
         #                     ),
         #                     conditionalPanel(
         #                       condition = "input.Signature == 'Upload a Signature'",
         #                       fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".csv"), width = NULL, 
         #                                 buttonLabel = "Browse...", placeholder = "No file selected")
         #                     ),
         #                     br(),
         #                     fileInput("AddedCompounds", "Add compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
         #                               buttonLabel = "Browse...", placeholder = "No file selected"),
         #                     uiOutput("mytest"),
         #                     # conditionalPanel(
         #                     #      condition = "input.AddedCompounds.name == null",
         #                     #         textInput("AddedLabel", "Label for added compounds", value="Added"),
         #                     #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
         #                     # ),
         #                     br(),
         #                     #selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
         #                     br(),
         #                     selectInput("Advanced", "See Advanced Options", list("Yes", "No"), selected = "No"),
         #                     conditionalPanel(condition="input.Advanced == 'Yes'",
         #                                      selectInput("Algorithm","Algorithm", list("minSim", "fpSim"), selected = "minSim"),
         #                                      selectInput("ConOrDiscon","Activation or Inhibition", list("Concordant", "Discordant"), selected = "Concordant"),
         #                                      numericInput("Concordance", "Concordance Threshold", value=0.2)
         #                     ),        
         #                     #conditionalPanel(
         #                     #     condition = "if(!is.null(input.AddedCompounds))",
         #                     #     checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
         #                     #),     
         #                     #br(),
         #                     actionButton(label="Go!", "Go")
         #                     ############# Add this ######################################################################################                  
         #                   )),
         #            #fluidRow(
         #            column(9, 
         #                   htmlOutput("notify1"),
         #                   visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
         #                   htmlOutput("notify2")
         #                   
         #            ))
         #          
         #          
         #          # )
         #          
         #          #############################################################################################################
         #          
         #          # )
#), #up to here
#DT::dataTableOutput("Combined_Score"),
#downloadButton("Max_Scores", label = "Download Scores")),
#        #downloadButton("NCIDownload", label = "Download Similar NCI"))

      tabPanel("Chemical Similarity Analysis", value="heatmap", icon = NULL,
               #br(),
               #br(),
          sidebarPanel(position="left", width=5,
               
               #plotOutput("distPlot", width="850px", height="800px")
      #          fluidRow(
      #            column(3,
      #                   wellPanel(
      # 
      #                     #####################################################################
      #                     br(),
      #                     selectInput("Signature", NULL, list("Input a Gene", "Upload a Signature", "Similarity Search"), selected = "Input a Gene"),
      #                     br(),
      # 
      #                     conditionalPanel(
      #                       condition = "input.Signature == 'Input a Gene'",
      #                       textInput("gene_knockdown", "Input a Gene", value="bcl2a1")
      #                     ),
      #                     conditionalPanel(
      #                       condition = "input.Signature == 'Upload a Signature'",
      #                       fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".csv"), width = NULL,
      #                                 buttonLabel = "Browse...", placeholder = "No file selected")
      #                     ),
      #                     br(),
      #                     fileInput("AddedCompounds", "Add compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
      #                               buttonLabel = "Browse...", placeholder = "No file selected"),
      #                     uiOutput("mytest"),
      #                     # conditionalPanel(
      #                     #      condition = "input.AddedCompounds.name == null",
      #                     #         textInput("AddedLabel", "Label for added compounds", value="Added"),
      #                     #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
      #                     # ),
      #                     br(),
      #                     #selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
      #                     br(),
      #                     selectInput("Advanced", "See Advanced Options", list("Yes", "No"), selected = "No"),
      #                     conditionalPanel(condition="input.Advanced == 'Yes'",
      #                                      selectInput("Algorithm","Algorithm", list("fpSim","minSim"), selected = "fpSim"),
      #                                      selectInput("ConOrDiscon","Activation or Inhibition", list("Concordant", "Discordant"), selected = "Concordant"),
      #                                      numericInput("Concordance", "Concordance Threshold", value=0.2)
      #                     ),
      # 
      #                     actionButton(label="Go!", "Go")
      #                     ############# Add this ######################################################################################
      #                   )),
      #            #fluidRow(
      #            column(9,
      #                   htmlOutput("notify1"),
      #                   visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
      #                   htmlOutput("notify2")
      # 
      #            ))
      # 
      # 
      #          # )
      # 
      #          #############################################################################################################
      # 
      #          # )
      # ), #up to here
                actionButton("Cluster", label="Run SAR"),
                
                #br(),
                #br(),
#               actionButton(label="GenerateHeatmap", "GenerateHeatmap"),
               #textInput("CutHeight", "Tanimoto Similarity", value="0.75"),
               #textInput("ClusterSize", "Minimum Cluster Size", value="3"),
               #actionButton(label="MultiDimensional Scaling", "GetRepresentatives"),
               #br(),
               #br(),
               #br(),
                  actionButton(label="Show Advanced Options", "Advanced_button_mds"),
                  br(),
                  uiOutput("downloadReps"),
                  br(),
                 #br(),
                  uiOutput("downloadClusters"),
                  br(),
                  conditionalPanel(condition="input.Advanced_button_mds ==  1",
                      textInput("CutHeight", "Tanimoto Similarity", value="0.75"),
                      textInput("ClusterSize", "Minimum Cluster Size", value="3"),
                      checkboxInput("cluster_all_cmpds", value=FALSE, label="Cluster all Compounds"),
                  br(),
                  br()
                  
                 
               #downloadButton("RepDownload", label="Download Representatives"),
               #downloadButton("ClusterDownload", label = "Download Clusters"),

          )),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),

          column(12,
               #plotOutput("distPlot", width = "850px", height = "800px"),
               #plotOutput("MDSPlot", width = "900px", height = "800px"),
               #uiOutput("cluster_check"),
               uiOutput("check_analogs"),
               uiOutput("check_clustersize"),
               uiOutput("check_clusternumber"),
               splitLayout(cellWidths = c("50%", "50%"),
               plotOutput("distPlot", width = "500px", height = "500px"),
               plotOutput("MDSPlot", width = "500px", height = "500px")),
               DT::dataTableOutput("Representatives")
               #downloadButton("RepDownload", label="Download Representatives"),
               #downloadButton("ClusterDownload", label = "Download Clusters"),
               #actionButton(label="Get Related NCI Compounds", "GetNCI")
               
              
          )
     
      ),
      # tabPanel("MDS Plot", value="MDS", icon=NULL,
      #          plotOutput("MDSPlot", width = "900px", height = "800px"),
      #          DT::dataTableOutput("Representatives"),
      #          downloadButton("RepDownload", label="Download Representatives"),
      #          downloadButton("ClusterDownload", label = "Download Clusters"),
      #          actionButton(label="Get Related NCI Compounds", "GetNCI")
      #          ),


#      tabPanel("Similar Compounds", value="NCI", icon=NULL,
       #       DT::dataTableOutput("Similar_NCI"),
         #     downloadButton("NCIDownload", label = "Download Similar NCI")),
#### add this script to UI #### 

tabPanel("Network Analysis", value = "mSTITCH", icon = NULL,
         #fluidRow(
        #   column(3, wellPanel(
          sidebarPanel( 
            selectInput("STITCH", NULL, list("Global STITCH", "Cluster Network STITCH"), multiple=FALSE),
            
            conditionalPanel(
              condition = "input.STITCH == 'Global STITCH'",
              checkboxInput("all_clusters", "Shows all clusters", value = FALSE, width = '400px'),
            actionButton(label="View Global STITCH Network", "GlobalSTITCH")
             ),
            br(),
            br(),
            br(),
            conditionalPanel(
              condition = "input.STITCH == 'Cluster Network STITCH'",
            textInput("ClusterNumber", "Cluster Number", value = 1),
            br(),
            selectInput("confidence", "Confidence Score", choices = c(150, 400, 700, 950), selected = 400),
            br(),
            sliderInput("Connections", "Numbers of Connections", min = 2, max = 25, value = 10, ticks = FALSE),
            br(),
            # textInput("Gene", "Gene of Interest", value = "bcl2a1"),
            checkboxInput("fixedNet", "Fixed node position", value = FALSE, width = '400px'),
            checkboxInput("show_undetectNodes", "Shows all", value = FALSE, width = '400px'),
            actionButton(label="View Selected Cluster Netowrk", "GetSTITCH")
            ) 
            
           #)),
           #column(9, 
            #        br(),      
            #      textOutput("numberClust"),    
           #)
         #),
         #fluidRow(
          #   column(12, 
           #       visNetworkOutput("mSTITCHPlot", width = "1000px", height = "750px"),
          # )
         ),
         mainPanel( 
           uiOutput("check_clustering"),  
           textOutput("numberClust"),
           visNetworkOutput("mSTITCHPlot", width = "1000px", height = "750px"),
           
           br(),
           br(),
           
           visNetworkOutput("STITCHPlot", width = "900px", height = "600px")
           
           #br(),
           #br()
           
          )
          )

# tabPanel("Selected Cluster Network", value = "STITCH", icon = NULL,
#          fluidRow(
#            column(2, wellPanel(
#              textInput("ClusterNumber", "Cluster Number", value = 1),
#              br(),
#              selectInput("confidence", "Confidence Score", choices = c(150, 400, 700, 950), selected = 400),
#              br(),
#              sliderInput("Connections", "Numbers of Connections", min = 2, max = 25, value = 10, ticks = FALSE),
#              br(),
#              # textInput("Gene", "Gene of Interest", value = "bcl2a1"),
#              checkboxInput("fixedNet", "Fixed node position", value = FALSE, width = '400px'),
#              checkboxInput("show_undetectNodes", "Shows all", value = FALSE, width = '400px'),
#              actionButton(label="View with STITCH", "GetSTITCH")
#            )),
#            column(9, 
#                   br(),
#                   visNetworkOutput("STITCHPlot", width = "900px", height = "600px")
#            )
#          ))


############################
        )
      )
  )
)


