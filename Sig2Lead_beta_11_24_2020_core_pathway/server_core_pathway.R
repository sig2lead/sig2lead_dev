#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

#########Load all required libraries here, they cannot be sourced
library(shiny)
library(httr)
library(jsonlite)
library(DT)
#library(heatmaply)
#library(shinyHeatmaply)
library(visNetwork)
library(ChemmineOB)
library(ChemmineR)
library(plyr)
library("dendextend")
library("colorspace")
library(ggforce)
library(rlist)
library(scatterpie)
library(ggrepel)
library(bazar)
library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(igraph)

#source("http://bioconductor.org/biocLite.R")

#biocLite("ChemmineOB")
#biocLite("ChemmineR")

load("./minSim_apfp_RObjects.RData")
source("lib/chemmineR_2_option.R")
source("./lib/cluster_fpsim.R")
######Discordant or Concordant
######Benchmarking -- Behrouz

begin_time <<- Sys.time()
options(shiny.maxRequestSize=200*1024^2)

shinyServer(function(input, output) {
  
  #output$mytest <- renderUI(
  #  if(!is.null(input$AddedCompounds)) {
  #    div(textInput("AddedLabel", "Label for Added Compounds", value="My Candidates"),
  #        checkboxInput("Bypass", "Bypass Clustering", value=TRUE)
  #    )
  #  }
  #)

 observeEvent(input$Go,
#  output$Go <- eventReactive(input$Go,  
    {
      Bypass = TRUE 
      #withProgress(message = "Running Sig2Lead", value = 0, {
        source("lib/GeneKD.R", local = TRUE)
        source("lib/ConDisconSignatures.R", local = TRUE)
        #lsm_rows <- 1:41572
      #1. Query iLINCS for user input gene
        if (input$Signature == "Input a Gene"){
          define_knockdown(input$gene_knockdown)
          print("Found your knockdowns") 
        
          if(input$ConOrDiscon=="Concordant"){
          incProgress(1/7, message = "Finding Concordant Signatures")
          ConDisconSignatures(sigid)}
          
          ConDisconSignatures <<- ConDisconSignatures(sigid)
          
          ##############################################################
          ##############################################################
          # Redo with LINCS fingerprints loaded
          load("./lincs_fps.RData")
          
          
          source("lib/GetSMILES.R", local = TRUE)
          #3. Generate SMILES of compounds
             GetSMILES()
             #print("SMILES acquired")
             #incProgress(1/7, message = "Converting SMILES to SDF...Slow")
          # #4. Generate Data Table
             display <<- final_compounds
             colnames(display) <- c("LSM_ID", "SMILES")
          #   
            #output$SMILES <<- DT::renderDataTable(display)
            #output$SMILES <<- DT::renderDataTable({DT::datatable(display, selection = list(selection = "single", target = "cell"))}, escape = FALSE)
          #   
          #   makeLink <- function(val) {
          #     #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          #     paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          #   }
          #   
          #   makeLink2 <- function(val) {
          #     #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          #     paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
          #   }
          #   #lincs <<-  data.frame(matrix(unlist(con_comma[[1]]), nrow= 457, byrow=T), stringsAsFactors = FALSE)
          #   #lincs <<- con[order(con[,2], decreasing = TRUE),c(11,6,2)]
          #   #colnames(lincs) <- c("LSM-ID", "Compound", "Concordance")
          #   
          #   makeLink3 <- function(val) {
          #     #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          #     paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          #   }
          #   
          #   makeLink4 <- function(val) {
          #     #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          #     paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
          #   }
          ##################################################################################
          ##################################################################################
          # Need to tweak with LINCS fingerprints pre-loaded
           df <- list()
           for (d in 1:length(con_comma)){
           
             df[[d]] <- as.data.frame(matrix(unlist(con_comma[[d]]), nrow=length(con_comma[[d]][[1]])))
           
           }
           df_concordances <<- ldply(df, data.frame)
           df_concordances_2 <- df_concordances[,c(11,6,2,3)]
           colnames(df_concordances_2) <- c("LSM-ID", "Compound", "Concordance", "Significance")
           df_concordances_2$Concordance <- as.numeric(as.character(df_concordances_2$Concordance))
           df_concordances_o <<- df_concordances_2[order(df_concordances_2$Concordance, decreasing = TRUE ),]
          # 
          # 
           df_concordances_no_dups <<- df_concordances_o[!duplicated(df_concordances_o["LSM-ID"]),]
           display_no_dups <- display[!duplicated(display["LSM_ID"]),]
          # 
           df_concordances_o_smiles <- merge(df_concordances_no_dups, display_no_dups, by.x="LSM-ID", by.y="LSM_ID", all.x=FALSE, all.y=FALSE)
           df_concordances_o_smiles<<- df_concordances_o_smiles[order(df_concordances_o_smiles$Concordance, decreasing = TRUE ), c(1,5,2,3,4)]
          # 
          # #############Output to LINCS Compounds Table;
          # 
          # #################
          # # Max Concordance
          # #################
           max_cons <- df_concordances_o_smiles[order(df_concordances_o$`LSM-ID`, -abs(df_concordances_o$Concordance)), ]
           max_cons_2 <- max_cons[!duplicated(max_cons$`LSM-ID`),]
           max_cons_3 <<- max_cons_2[order(max_cons_2$Concordance, decreasing = TRUE), -5]
           max_cons_4 <<- max_cons_3[-which(is.na(max_cons_3$`LSM-ID`)),]
           con_df <- as.data.frame((matrix(unlist(con_comma),byrow=T)))
           max_cons_4$Concordance <- round(max_cons_4$Concordance,3)
           output$SMILES <<- DT::renderDataTable(max_cons_4[,-4], caption="LINCS Candidates", rownames=FALSE)
           load("./lincs_fps.RData")
          # 
           lsm_rows <<- which(rownames(lincs_fps_2) %in% max_cons_4$`LSM-ID`)
           #output$SMILES <- DT::renderDataTable({
          # 
          #   #display$link <- makeLink(display$LSM_ID)
          #   #display$LSM_ID <- makeLink(display$LSM_ID)
          #   #df_concordances_o_smiles$`LSM-ID`<- makeLink(df_concordances_o_smiles$`LSM-ID`)
          #   #max_cons_3$`LSM-ID`<- makeLink(max_cons_3$`LSM-ID`)
          #   #return(display)
          #   #return(display)
          #   #return(df_concordances_o_smiles)
          #   return(max_cons_3)
          # },
          # escape = FALSE, rownames = FALSE)
          # 
           #http://lincsportal.ccs.miami.edu/SmallMolecules/view/LSM-5467
           display <<- display
          # 
           source("lib/ChemmineOB.R", local = TRUE)
          # 
          # ###This is the slow step, see if we can go straight from smiles to fingerprints and seperately convert to SDF
           #print("Converting SMILES to SDF")
           Get_SDF()
          
          ###
          #   
          #   print("SDF conversion complete")
          #   incProgress(3/7, message = "Clustering Related Compounds")
          #5. Cluster compounds identified through LINCS, generate heatmaps
          # If Bypass Clu
          #source("lib/chemmineR.R", loca
        } 
        
        #####################################################################
        ##### Input List of Genes ###########################################
        #####################################################################
        if (input$Signature == "Input List of Genes"){
          #max_cons_4 <- list()
          #display <- list()
          for (k in 1:3){
            if (k == 1){
              define_knockdown(input$gene_knockdown_1)}
            else if (k == 2){
              define_knockdown(input$gene_knockdown_2)}
            else if (k == 3){
              define_knockdown(input$gene_knockdown_3)}
            
          print("Found your knockdowns") 
          
          if(input$ConOrDiscon=="Concordant"){
            incProgress(1/7, message = "Finding Concordant Signatures")
            ConDisconSignatures(sigid)}
          
          ConDisconSignatures <<- ConDisconSignatures(sigid)
          
          ##############################################################
          ##############################################################
          # Redo with LINCS fingerprints loaded
          load("./lincs_fps.RData")
          
          
          source("lib/GetSMILES.R", local = TRUE)
          #3. Generate SMILES of compounds
          GetSMILES()
          
          # #4. Generate Data Table
          display <- final_compounds
          colnames(display) <- c("LSM_ID", "SMILES")
          #   
          #output$SMILES <<- DT::renderDataTable(display)
          #output$SMILES <<- DT::renderDataTable({DT::datatable(display, selection = list(selection = "single", target = "cell"))}, escape = FALSE)
          #   
          #   makeLink <- function(val) {
          #     #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          #     paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          #   }
          #   
          #   makeLink2 <- function(val) {
          #     #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          #     paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
          #   }
          #   #lincs <<-  data.frame(matrix(unlist(con_comma[[1]]), nrow= 457, byrow=T), stringsAsFactors = FALSE)
          #   #lincs <<- con[order(con[,2], decreasing = TRUE),c(11,6,2)]
          #   #colnames(lincs) <- c("LSM-ID", "Compound", "Concordance")
          #   
          #   makeLink3 <- function(val) {
          #     #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          #     paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          #   }
          #   
          #   makeLink4 <- function(val) {
          #     #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          #     paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
          #   }
          ##################################################################################
          ##################################################################################
          # Need to tweak with LINCS fingerprints pre-loaded
          df <- list()
          for (d in 1:length(con_comma)){
            
            df[[d]] <- as.data.frame(matrix(unlist(con_comma[[d]]), nrow=length(con_comma[[d]][[1]])))
            
          }
          df_concordances <<- ldply(df, data.frame)
          df_concordances_2 <- df_concordances[,c(11,6,2,3)]
          colnames(df_concordances_2) <- c("LSM-ID", "Compound", "Concordance", "Significance")
          df_concordances_2$Concordance <- as.numeric(as.character(df_concordances_2$Concordance))
          df_concordances_o <<- df_concordances_2[order(df_concordances_2$Concordance, decreasing = TRUE ),]
          # 
          # 
          df_concordances_no_dups <<- df_concordances_o[!duplicated(df_concordances_o["LSM-ID"]),]
          display_no_dups <- display[!duplicated(display["LSM_ID"]),]
          # 
          df_concordances_o_smiles <- merge(df_concordances_no_dups, display_no_dups, by.x="LSM-ID", by.y="LSM_ID", all.x=FALSE, all.y=FALSE)
          df_concordances_o_smiles<<- df_concordances_o_smiles[order(df_concordances_o_smiles$Concordance, decreasing = TRUE ), c(1,5,2,3,4)]
          # 
          # #############Output to LINCS Compounds Table;
          # 
          # #################
          # # Max Concordance
          # #################
          max_cons <- df_concordances_o_smiles[order(df_concordances_o$`LSM-ID`, -abs(df_concordances_o$Concordance)), ]
          max_cons_2 <- max_cons[!duplicated(max_cons$`LSM-ID`),]
          max_cons_3 <<- max_cons_2[order(max_cons_2$Concordance, decreasing = TRUE), -5]
          max_cons_4 <- max_cons_3[-which(is.na(max_cons_3$`LSM-ID`)),]
          con_df <- as.data.frame((matrix(unlist(con_comma),byrow=T)))
          max_cons_4$Concordance <- round(max_cons_4$Concordance,3)
          #output$SMILES <<- DT::renderDataTable(max_cons_4[[k]][,-4], caption="LINCS Candidates", rownames=FALSE)
          load("./lincs_fps.RData")
          # 
          lsm_rows <<- which(rownames(lincs_fps_2) %in% max_cons_4$`LSM-ID`)
        
          source("lib/ChemmineOB.R", local = TRUE)
          # 
          # ###This is the slow step, see if we can go straight from smiles to fingerprints and seperately convert to SDF
          #print("Converting SMILES to SDF")
          Get_SDF()
          
          
          }
          stop()
        }  
        else if (input$Signature == "Upload a Signature"){
          source("lib/UploadSignature.R")
          uploadsig <- input$UploadSignature
          sigUpload(paste(getwd(), uploadsig$name, sep="/"))
          print("Found your signature")
          incProgress(1/7, message = "Finding Concordant Signatures")
          UploadConDiscon(rjson)
          
        }
        
        else if(input$Signature == "Similarity Search")
        {
          lsm_rows <- 1:41572
          #output$SMILES <<- NULL
          
          adds<<-input$AddedCompounds
          if(!(grepl(".sdf", adds$name))){
            adds_SMI <<- read.SMIset(adds$datapath)
            #adds_SMI<-read.SMIset(paste(getwd(), adds$name, sep="/"))
            #adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
            #adds_csv <- adds_csv[c(2,1)]
            #colnames(adds_csv) <- c("LSM_ID", "SMILES")
            #adds_csv <<- adds_csv
            #test2 <<- rbind(adds_csv, display)
            #colnames(test2) <<- c("Compound_ID", "SMILES")
            sdfset_add <<- smiles2sdf(adds_SMI)
            #sdf_smiles <<- c(sdf_smiles, sdfset_add)
            bypass_clustering(sdf_smiles, sdfset_add)}
          else if(grepl(".sdf", adds$name)){
            print("Added compounds in sdf format")
            max_cons_4 <<- NULL
            }
        }
          
            
          
        
        
        #ConDisconSignatures <<- ConDisconSignatures(sigid)
        #}
        
        
        #incProgress(1/7, detail = "Finding Concordant Signatures")
        
      #2. Get a list of all concordant compounds
        
        
       #print("Found concordant signatures")
      #   incProgress(1/7, message = "Obtaining SMILES")
      #   source("lib/GetSMILES.R", local = TRUE)
      # #3. Generate SMILES of compounds
      #   GetSMILES()
      #   print("SMILES acquired")
      #   incProgress(1/7, message = "Converting SMILES to SDF...Slow")
      # #4. Generate Data Table
      #   display <<- final_compounds
      #   colnames(display) <- c("LSM_ID", "SMILES")
      #   
      #   #output$SMILES <<- DT::renderDataTable(display)
      #   #output$SMILES <<- DT::renderDataTable({DT::datatable(display, selection = list(selection = "single", target = "cell"))}, escape = FALSE)
      #   
      #   makeLink <- function(val) {
      #     #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
      #     paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
      #   }
      #   
      #   makeLink2 <- function(val) {
      #     #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
      #     paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
      #   }
      #   #lincs <<-  data.frame(matrix(unlist(con_comma[[1]]), nrow= 457, byrow=T), stringsAsFactors = FALSE)
      #   #lincs <<- con[order(con[,2], decreasing = TRUE),c(11,6,2)]
      #   #colnames(lincs) <- c("LSM-ID", "Compound", "Concordance")
      #   
      #   makeLink3 <- function(val) {
      #     #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
      #     paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
      #   }
      #   
      #   makeLink4 <- function(val) {
      #     #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
      #     paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
      #   }

        # df <- list()
        # for (d in 1:length(con_comma)){
        # 
        #   df[[d]] <- as.data.frame(matrix(unlist(con_comma[[d]]), nrow=length(con_comma[[d]][[1]])))
        # 
        # }
        # df_concordances <<- ldply(df, data.frame)
        # df_concordances_2 <- df_concordances[,c(11,6,2,3)]
        # colnames(df_concordances_2) <- c("LSM-ID", "Compound", "Concordance", "Significance")
        # df_concordances_2$Concordance <- as.numeric(as.character(df_concordances_2$Concordance))
        # df_concordances_o <<- df_concordances_2[order(df_concordances_2$Concordance, decreasing = TRUE ),]
        # 
        # 
        # df_concordances_no_dups <<- df_concordances_o[!duplicated(df_concordances_o["LSM-ID"]),]
        # display_no_dups <- display[!duplicated(display["LSM_ID"]),]
        # 
        # df_concordances_o_smiles <- merge(df_concordances_no_dups, display_no_dups, by.x="LSM-ID", by.y="LSM_ID", all.x=FALSE, all.y=FALSE)
        # df_concordances_o_smiles<<- df_concordances_o_smiles[order(df_concordances_o_smiles$Concordance, decreasing = TRUE ), c(1,5,2,3,4)]
        # 
        # #############Output to LINCS Compounds Table;
        # 
        # #################
        # # Max Concordance
        # #################
        # max_cons <- df_concordances_o_smiles[order(df_concordances_o$`LSM-ID`, -abs(df_concordances_o$Concordance)), ]
        # max_cons_2 <- max_cons[!duplicated(max_cons$`LSM-ID`),]
        # max_cons_3 <<- max_cons_2[order(max_cons_2$Concordance, decreasing = TRUE), -5]
        # max_cons_4 <<- max_cons_3[-which(is.na(max_cons_3$`LSM-ID`)),]
        # #con_df <- as.data.frame((matrix(unlist(con_comma),byrow=T)))
        # load("./lincs_fps.RData")
        # 
        # lsm_rows <- which(rownames(lincs_fps_2) %in% max_cons_4$`LSM-ID`)
        # output$SMILES <- DT::renderDataTable({
        # 
        #   #display$link <- makeLink(display$LSM_ID)
        #   #display$LSM_ID <- makeLink(display$LSM_ID)
        #   #df_concordances_o_smiles$`LSM-ID`<- makeLink(df_concordances_o_smiles$`LSM-ID`)
        #   #max_cons_3$`LSM-ID`<- makeLink(max_cons_3$`LSM-ID`)
        #   #return(display)
        #   #return(display)
        #   #return(df_concordances_o_smiles)
        #   return(max_cons_3)
        # },
        # escape = FALSE, rownames = FALSE)
        # 
        # #http://lincsportal.ccs.miami.edu/SmallMolecules/view/LSM-5467
        # display <<- display
        # 
        # source("lib/ChemmineOB.R", local = TRUE)
        # 
        # ###This is the slow step, see if we can go straight from smiles to fingerprints and seperately convert to SDF
        # print("Converting SMILES to SDF")
        # Get_SDF()

        ###
      #   
      #   print("SDF conversion complete")
      #   incProgress(3/7, message = "Clustering Related Compounds")
      #5. Cluster compounds identified through LINCS, generate heatmaps
        # If Bypass Clu
        #source("lib/chemmineR.R", local = TRUE)
        
#    if(input$Bypass == FALSE){
        # Add Compounds and Cluster
    if(!is.null(input$AddedCompounds)){
      adds<<-input$AddedCompounds
      #5.1. Include Added Compounds if applicable
         #if (!is.null(input$AddedCompounds))      {
         if (Bypass == FALSE)      {                     
                      #adds<<-input$AddedCompounds
                      
                      if(!(grepl(".sdf", adds$name))){
                        adds_SMI <<- read.SMIset(adds$datapath)
                        adds_csv <- read.csv(adds$datapath, header = FALSE, sep = "\t")
                       # adds_SMI <<- read.SMIset(paste(getwd(), adds$name, sep="/"))
                        #adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
                        adds_csv <- adds_csv[,c(2,1)]
                        colnames(adds_csv) <- c("LSM_ID", "SMILES")
                        adds_csv <<- adds_csv
                        test2 <<- rbind(adds_csv, display)
                        colnames(test2) <<- c("Compound_ID", "SMILES")
                        
                        sdfset_add <<- smiles2sdf(adds_SMI)
                        sdf_smiles <<- c(sdf_smiles, sdfset_add)
                        #sdf_smiles <<- c(sdf_smiles, adds_SDF)
                        cluster_compounds()
                        print("Your compounds have been clustered with LINCS compounds")
                        source("lib/ColorMap.R", local = TRUE)
                        color_dend()
                      }
                      else if(grepl(".sdf", adds$name)){
                        print("Added compounds in sdf format")
                        #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
                        #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
                        lsm_smiles_2 <- lsm_smiles
                        sdfset_add <<- read.SDFset(adds$datapath)
                        added_smiles <<- sdf2smiles(sdfset_add)
                        added_labels <- cid(added_smiles)
                        added_smiles_2 <- as.character(added_smiles[1:length(added_smiles)])
                        added_smiles_3 <- cbind(added_labels, added_smiles_2)
                        
                        #added_smiles_4 <- unname(added_smiles_3)
                        colnames(added_smiles_3) <- c("Compound_ID", "SMILES")
                        
                        lsm_smiles_2_labels <- cid(lsm_smiles_2)
                        lsm_smiles_3 <- as.character(lsm_smiles_2[1:length(lsm_smiles_2)])
                        lsm_smiles_3b <- unname(lsm_smiles_3)
                        lsm_smiles_4 <- cbind(lsm_smiles_2_labels, lsm_smiles_3b)
                        colnames(lsm_smiles_4) <- c("Compound_ID", "SMILES")
                        colnames(display) <-  c("Compound_ID", "SMILES")
                        
                        test2 <<- rbind(added_smiles_3, display)
                        
                        colnames(test2) <<- c("Compound_ID", "SMILES")
                        
                        
                        sdf_smiles <<- c(sdf_smiles, sdfset_add)
                        cluster_compounds()
                        print("Your compounds have been clustered with LINCS compounds")
                        source("lib/ColorMap.R", local = TRUE)
                        adds_SMI <<- added_smiles
                        color_dend()
                        }
                      #color_dend() currently not rendering
                      #output$distPlot<<- renderPlot(heatmap.2(dist_mat, Rowv=dend, Colv=dend, colRow = heatmap_colors, colCol = heatmap_colors, col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
                }#,
      #5.2. Don't add compounds
      #5.2. Add compounds, but bypass clustering
        else
                {
                  adds<<-input$AddedCompounds
                  if(!(grepl(".sdf", adds$name))){
                    adds_SMI <<- read.SMIset(adds$datapath)
                    #adds_SMI<-read.SMIset(paste(getwd(), adds$name, sep="/"))
                    #adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
                    #adds_csv <- adds_csv[c(2,1)]
                    #colnames(adds_csv) <- c("LSM_ID", "SMILES")
                    #adds_csv <<- adds_csv
                    #test2 <<- rbind(adds_csv, display)
                    #colnames(test2) <<- c("Compound_ID", "SMILES")
                    sdfset_add <<- smiles2sdf(adds_SMI)
                    #sdf_smiles <<- c(sdf_smiles, sdfset_add)
                    bypass_clustering(sdf_smiles, sdfset_add)}
                  else if(grepl(".sdf", adds$name)){
                    print("Added compounds in sdf format")
                    # inFile <- input$AddedCompounds
                    # old_name <- inFile$datapath
                    # dirstr <- dirname(inFile$datapath)
                    # new_name <- paste(dirstr, inFile$name, sep="/")
                    # file.rename(old_name,new_name)
                    #sdfset_add <<- read.SDFset(adds$name)
                    #print(dim(sdfset_add))
                    #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                    #added_labels <- sdfid(sdfset_add)
                    
                    ####################################
                    # ChemmineR Atom Pair Fingerprints
                    ##################################
                    #apset_added <-sdf2ap(sdfset_add)
                    #fpset_added <- as.matrix(desc2fp(apset_added))
                    
                    #sdfset_add <<- as.matrix(fpset_added)
                    #print(dim(sdfset_add))
                    #sdf_smiles <<- c(sdf_smiles, sdfset_add[1:10])
                    
                    
                    # If fpsim <- UI widget == fpsim
                    if ((input$Algorithm == "fpSim") & (input$Signature == "Input a Gene")){
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      test_bypass_fpsim <-  bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      test_bypass_fpsim <- cbind(added_labels, test_bypass_fpsim)
                      colnames(test_bypass_fpsim) <- c("Compound", "LSM_ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      for (i in 1:nrow(test_bypass_fpsim)){
                        Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      }
                      test_bypass_fpsim$Concordance <- Concordance
                      test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      #test_bypass_fpsim <<- test_bypass_fpsim
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
                      #stopApp()
                    }
                    
                    #else minsim
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Input a Gene"))
                    { #load("./minSim_apfp_RObjects.RData")
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      test_bypass_minsim <- bypass_clustering_minsim(sdfset_add, lsm_rows,max_cons_4)
                      test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      colnames(test_bypass_minsim) <- c("Compound", "LSM_ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      for (i in 1:nrow(test_bypass_minsim)){
                        Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      }
                      test_bypass_minsim$Concordance <- Concordance
                      test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,4], test_bypass_minsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_minsim <<- test_bypass_minsim
                      #stopApp()
                    }
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Similarity Search"))
                    { #load("./minSim_apfp_RObjects.RData")
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      test_bypass_minsim <- bypass_clustering_minsim(sdfset_add, lsm_rows)
                      test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      colnames(test_bypass_minsim) <- c("Compound", "LSM_ID", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      #Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance 
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      #stopApp()
                    }
                    
                    else if ((input$Algorithm == "fpSim") & (input$Signature == "Similarity Search"))
                    { #load("./minSim_apfp_RObjects.RData")
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      test_bypass_fpsim <- bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      test_bypass_fpsim <- data.frame(added_labels, test_bypass_fpsim)
                      colnames(test_bypass_fpsim) <- c("Compound", "LSM_ID", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      #Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance 
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                    }
                    # cluster_compounds()
                  # print("LINCS clustering complete!")
                  # output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
                }
    }#)
    }
#    else{print("Clustering bypassed")
    else{print("No added compounds, must cluster.")
      #if (!is.null(input$AddedCompounds))      {
        
        cluster_compounds()
        print("LINCS clustering complete!")
        output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
        # adds<<-input$AddedCompounds
        # if(!(grepl(".sdf", adds$name))){
        #   adds_SMI<-read.SMIset(paste(getwd(), adds$name, sep="/"))
        #   adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
        #   adds_csv <- adds_csv[c(2,1)]
        #   colnames(adds_csv) <- c("LSM_ID", "SMILES")
        #   adds_csv <<- adds_csv
        #   test2 <<- rbind(adds_csv, display)
        #   colnames(test2) <<- c("Compound_ID", "SMILES")
        #   adds_SDF<-smiles2sdf(adds_SMI)
        #   sdf_smiles <<- c(sdf_smiles, adds_SDF)
        #   bypass_clustering()}
        # else if(grepl(".sdf", adds$name)){
        #   print("Added compounds in sdf format")
        #   sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
        #   sdf_smiles <<- c(sdf_smiles, sdfset_add[1:10])
        #   bypass_clustering()
          
          
        #}  
        #print("Your compounds have been clustered with LINCS compounds")
        #source("lib/ColorMap.R", local = TRUE)
        #color_dend()
        #color_dend() currently not rendering
        #output$distPlot<<- renderPlot(heatmap.2(dist_mat, Rowv=dend, Colv=dend, colRow = heatmap_colors, colCol = heatmap_colors, col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
      #}#,
      #5.2. Don't add compounds
      # else
      # {
      #   bypass_clustering()
      #   #print("LINCS clustering complete!")
      #   #output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
      # }
      
    }
      
      
      
      
      
          
        
      #5.2 Get combined score
        # output Combined Scores
        #source("./lib/Combined_Score_rev_no_auto.R")
        # source("./lib/Combined_Score_auto_2.R")
        # source("./lib/Combined_Score_auto_3.R")
        # if (!is.null(input$AddedCompounds)){
        #   if(input$Bypass == TRUE){   
        #       combined_score_output <- combined_score_2()
        #       output$Combined_Score <- DT::renderDataTable({
        #       combined_score_output$`LSM-ID`<- makeLink(combined_score_output$`LSM-ID`)
        #       combined_score_output$`My Compounds`<- makeLink2(combined_score_output$`My Compounds`)
        #       return(combined_score_output)
        #       #df_concordances_o$`LSM-ID`<- makeLink(df_concordances_o$`LSM-ID`)
        #       #return(display)
        #       #return(df_concordances_o)
        #       
        # },
        # escape = FALSE, rownames = FALSE)  
        #   }else{
        #     print("Combined Score with Clustering")
        #     combined_score_output <- combined_score_3()
        #     output$Combined_Score <- DT::renderDataTable({
        #     
        #       combined_score_output$`LSM-ID`<- makeLink3(combined_score_output$`LSM-ID`)
        #       combined_score_output$`My Compounds`<- makeLink4(combined_score_output$`My Compounds`)
        #       return(combined_score_output)}, escape=FALSE, rownames=FALSE)
        #   }} 
        # else{
        #   print("No added compounds so can't calculate combined score")
        # }
        # 
#        cluster_compounds()
        
#        output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
       incProgress(7/7, message = "Complete!")
    #}) Progress bar to lin 62
    }) # End input$Go
  
  observeEvent(input$Cluster,{
    withProgress(message = "Clustering", value = 0, {
    if(!(grepl(".sdf", adds$name))){
      adds_SMI <<- read.SMIset(adds$datapath)
      adds_csv <- read.csv(adds$datapath, header = FALSE, sep = "\t")
      # adds_SMI <<- read.SMIset(paste(getwd(), adds$name, sep="/"))
      #adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
      adds_csv <- adds_csv[,c(2,1)]
      colnames(adds_csv) <- c("LSM_ID", "SMILES")
      adds_csv <<- adds_csv
      test2 <<- rbind(adds_csv, display)
      colnames(test2) <<- c("Compound_ID", "SMILES")
      
      sdfset_add <<- smiles2sdf(adds_SMI)
      sdf_smiles <<- c(sdf_smiles, sdfset_add)
      #sdf_smiles <<- c(sdf_smiles, adds_SDF)
      fpsim_cluster(fps)
      print("Your compounds have been clustered with LINCS compounds")
      source("lib/ColorMap.R", local = TRUE)
      color_dend()
    }
    else if(grepl(".sdf", adds$name)){
      print("Added compounds in sdf format")
      #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
      #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
      #lsm_smiles_2 <- lsm_smiles
      #sdfset_add <<- read.SDFset(adds$datapath)
      #added_smiles <<- sdf2smiles(sdfset_add)
      #added_labels <- cid(added_smiles)
      #added_smiles_2 <- as.character(added_smiles[1:length(added_smiles)])
      #added_smiles_3 <- cbind(added_labels, added_smiles_2)
      
      #added_smiles_4 <- unname(added_smiles_3)
      #colnames(added_smiles_3) <- c("Compound_ID", "SMILES")
      
      #lsm_smiles_2_labels <- cid(lsm_smiles_2)
      #lsm_smiles_3 <- as.character(lsm_smiles_2[1:length(lsm_smiles_2)])
      #lsm_smiles_3b <- unname(lsm_smiles_3)
      #lsm_smiles_4 <- cbind(lsm_smiles_2_labels, lsm_smiles_3b)
      #colnames(lsm_smiles_4) <- c("Compound_ID", "SMILES")
      #colnames(display) <-  c("Compound_ID", "SMILES")
      
      #test2 <<- rbind(added_smiles_3, display)
      
      #colnames(test2) <<- c("Compound_ID", "SMILES")
      
      
      #sdf_smiles <<- c(sdf_smiles, sdfset_add)
      #cluster_compounds()
      #print("Your compounds have been clustered with LINCS compounds")
      
      ############# Hierarchical CLustering ###########################
      source("lib/ColorMap.R", local = TRUE)
      #adds_SMI <<- added_smiles
      #color_dend()
      
      adds<<-input$AddedCompounds
      sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
      added_labels <- sdfid(sdfset_add)
      
      ####################################
      # ChemmineR Atom Pair Fingerprints
      ##################################
      apset_added <-sdf2ap(sdfset_add)
      fpset_added <- as.matrix(desc2fp(apset_added))
      
      #sdfset_add <- as.matrix(fpset_added)
      #print(dim(sdfset_add))
      #lsm_rows <- 1:100
      fpsim_cluster(fpset_added, lsm_rows)
      #incProgress(1/2, detail = "Clustering")
      #test_cluster_fpsim <- cbind(added_labels, test_bypass_fpsim)
      #colnames(test_bypass_fpsim) <- c("Compound", "LSM_ID", "Similarity")
      #test_bypass_fpsim <<- test_bypass_fpsim
      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
      output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), 
                                               Colv=as.dendrogram(hc), 
                                               col=colorpanel(40, "white","yellow","red"), 
                                               density.info="none", 
                                               trace="none", 
                                               #labCol=rownames(fpset), 
                                               #labRow=rownames(fpset), 
                                               cexRow=0.5, 
                                               cexCol=0.5))
    
    
      incProgress(1/1, message = "Hierarchical Clustering Complete")
    }
      ########## MDS Plot #########################################################################
      added_smiles <<- sdf2smiles(sdfset_add)
      added_labels <- cid(added_smiles)
      added_smiles_2 <- as.character(added_smiles[1:length(added_smiles)])
      added_smiles_3 <- cbind(added_labels, added_smiles_2)
      
      added_smiles_4 <- unname(added_smiles_3)
      colnames(added_smiles_3) <- c("Compound_ID", "SMILES")
      
      #lsm_smiles_2_labels <- cid(lsm_smiles_2)
      #lsm_smiles_3 <- as.character(lsm_smiles_2[1:length(lsm_smiles_2)])
      #lsm_smiles_3b <- unname(lsm_smiles_3)
      #lsm_smiles_4 <- cbind(lsm_smiles_2_labels, lsm_smiles_3b)
      #colnames(lsm_smiles_4) <- c("Compound_ID", "SMILES")
      colnames(display) <-  c("Compound_ID", "SMILES")
      
      test2 <<- rbind(added_smiles_3, display)
      
      colnames(test2) <<- c("Compound_ID", "SMILES")
      
      Added_Label <<- input$AddedLabel
      withProgress(message = "Plotting Representatives", value=0, {
        #source("lib/Centroid.R", local=TRUE)
        source("lib/Centroid2.R", local=TRUE)
        #Add a threshold option for cut_tree and Cluster size
        #6. Cut tree at user specified similarity score, identify representative based on minimum distance from all other members
        cut_tree((1-as.numeric(input$CutHeight)))
        #6.1. Include Added compounds
        if (!is.null(input$AddedCompounds)){
          find_centroid_adds(ClusterMembers, (as.numeric(input$ClusterSize)-1))
        }
        #6.2. No compounds were added
        else{
          find_centroid(ClusterMembers, (as.numeric(input$ClusterSize)-1))
        }
        df_centroid <- as.data.frame(unlist(centroid), stringsAsFactors = FALSE)
        colnames(df_centroid)="Representative"
        df_centroid<-df_centroid
        centroid_clusters(centroid, ClusterMembers)
        df_centroid <- cbind(centroid = df_centroid, Cluster = unlist(ClusterCentroid))
        
        
        makeLink <- function(val) {
          #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
        }
        
        makeLink2 <- function(val) {
          #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
        }
        
        
        
        
        
        if (!is.null(input$AddedCompounds)){
          test5 <- data.frame("Compound_ID"=0, "SMILES"=0)
          for (i in 1:length(df_centroid$Representative)){    
            if (length(which(test2$Compound_ID %in% df_centroid$Representative[[i]]))>0){
              test3 <- which(test2$Compound_ID %in% df_centroid$Representative[[i]])
              for (j in 1:length(test3)){
                test4 <- test2[test3[[j]], c("Compound_ID", "SMILES")]
                if (test5$Compound_ID[[1]]!=0){
                  test5 <- rbind(test5, test4)
                }
                else{
                  test5 <- test4
                }
              }
            }
          }
        }
        else {
          test5 <- data.frame("LSM_ID"=0, "SMILES"=0)
          for (i in 1:length(df_centroid$Representative)){    
            if (length(which(display$LSM_ID %in% df_centroid$Representative[[i]]))>0){
              test3 <- which(display$LSM_ID %in% df_centroid$Representative[[i]])
              for (j in 1:length(test3)){
                test4 <- display[test3[[j]], c("LSM_ID", "SMILES")]
                if (test5$LSM_ID[[1]]!=0){
                  test5 <- rbind(test5, test4)
                }
                else{
                  test5 <- test4
                }
              }
            }
          }
        }
        colnames(test5) <- c("Representative", "SMILES")
        test5 <<- test5
        
        
        df_centroid <- merge(test5, df_centroid, stringsAsFactors = FALSE)
        df_centroid <- data.frame(lapply(df_centroid, as.character), stringsAsFactors = FALSE)
        df_centroid$Cluster <- as.integer(df_centroid$Cluster)
        df_centroid <- df_centroid[order(df_centroid$Cluster),]
        #df_centroid<<-df_centroid
        df_centroid <- df_centroid[,-2]
        for(i in 1:nrow(df_centroid)){
            df_centroid$Compound[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == df_centroid$Representative[i]), "Compound"]
        } 
                     
                     
        df_centroid <- df_centroid[,c(1,3,2)]
        df_centroid <<- df_centroid
                 
        
        output$Representatives <<- DT::renderDataTable({
          for (i in 1:length(df_centroid$Representative)){
            if (!is.na((grep("LSM", df_centroid$Representative[[i]])) && (grep("LSM", df_centroid$Representative[[i]])==1))){
              df_centroid$Representative[[i]] <- makeLink(df_centroid$Representative[[i]])
            }
            else if(!is.null(grep("LSM", df_centroid$Representative[[i]]))){
              df_centroid$Representative[[i]] <- gsub("ZINC0", "ZINC", df_centroid$Representative[[i]])              
              df_centroid$Representative[[i]] <- makeLink2(df_centroid$Representative[[i]])
            }
          }
          return(df_centroid)
          #df_centroid <- df_centroid[,-2]
          df_centroid <<- df_centroid
        }, escape = FALSE, rownames = FALSE)
        #print(df_centroid$Representative)
        
        
        print(paste("Got the centroids as ", centroid[[1]], sep=""))
        incProgress(1/2, detail = "Plotting MDS")
        
        source("lib/MDS.R")
        #7. Generate MDS plot from representatives distance to one another, pie chart radius corresponds to cluster size
        MDS_plot(centroid)
        #            centroid_clusters(centroid, ClusterMembers)
        #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
        #                             + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
        output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
                                     + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "My Candidates")))
        #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid), hjust=1, vjust=2)) 
        #                           + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
        incProgress(2/2, detail = "Complete!")                                                         
      })
      
      })
    
  }
  )
  
  
    output$CandidatesDownload <- downloadHandler(
      filename = function() {
        'candidate_compounds.csv'
      },
      content = function(filename){
      if(exists("test_bypass_minsim")){
        write.csv(test_bypass_minsim[order(test_bypass_minsim[,4], test_bypass_minsim[,3], decreasing=TRUE),], filename, row.names=FALSE)
      }
      else{
        write.csv(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], filename, row.names=FALSE)
      }
      }
      )
  
  output$LINCSDownload <- downloadHandler(
    filename = function() {
      paste(input$gene_knockdown, '_lincs_candidates.csv', sep='')
    },
    content = function(filename){
      write.csv(max_cons_4[,-4], filename)
    }
  )
  
  # output$SMILESDownload <- downloadHandler(
  #   filename = function() {
  #     paste(input$gene_knockdown, '_compounds.csv', sep='')
  #   },
  #   content = function(filename){
  #     write.csv(final_compounds, filename)
  #   }
  #   )
  # 
  # output$SDFDownload <- downloadHandler(
  #   filename = function() {
  #     paste(input$gene_knockdown, '.sdf', sep='')
  #   },
  #   content = function(filename){
  #     write.SDF(sdf_smiles, filename)
  #   }
  # )
  # 
  # output$ConTableDownload <- downloadHandler(
  #   
  #   filename = "concordance_table.csv",
  #   content = function(filename){
  #     write.csv(df_concordances_o_smiles, filename)
  #   }
  # )
#  output$Representatives <- eventReactive(input$GetRepresentatives,
#   observeEvent(input$GetRepresentatives,
#         {
#           Added_Label <<- input$AddedLabel
#           withProgress(message = "Plotting Representatives", value=0, {
#             #source("lib/Centroid.R", local=TRUE)
#             source("lib/Centroid2.R", local=TRUE)
#             #Add a threshold option for cut_tree and Cluster size
#           #6. Cut tree at user specified similarity score, identify representative based on minimum distance from all other members
#             cut_tree((1-as.numeric(input$CutHeight)))
#           #6.1. Include Added compounds
#             if (!is.null(input$AddedCompounds)){
#               find_centroid_adds(ClusterMembers, (as.numeric(input$ClusterSize)-1))
#             }
#           #6.2. No compounds were added
#             else{
#             find_centroid(ClusterMembers, (as.numeric(input$ClusterSize)-1))
#             }
#             df_centroid <- as.data.frame(unlist(centroid), stringsAsFactors = FALSE)
#             colnames(df_centroid)="Representative"
#             df_centroid<-df_centroid
#             centroid_clusters(centroid, ClusterMembers)
#             df_centroid <- cbind(centroid = df_centroid, Cluster = unlist(ClusterCentroid))
#             
#             
#             makeLink <- function(val) {
#               #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
#               paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
#             }
#             
#             makeLink2 <- function(val) {
#               #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
#               paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
#             }
#             
#             
#             
#       
# 
#             if (!is.null(input$AddedCompounds)){
#               test5 <- data.frame("Compound_ID"=0, "SMILES"=0)
#             for (i in 1:length(df_centroid$Representative)){    
#                           if (length(which(test2$Compound_ID %in% df_centroid$Representative[[i]]))>0){
#                             test3 <- which(test2$Compound_ID %in% df_centroid$Representative[[i]])
#                                  for (j in 1:length(test3)){
#                                    test4 <- test2[test3[[j]], c("Compound_ID", "SMILES")]
#                                    if (test5$Compound_ID[[1]]!=0){
#                                    test5 <- rbind(test5, test4)
#                                    }
#                                    else{
#                                      test5 <- test4
#                                    }
#                                    }
#                           }
#             }
#             }
#             else {
#               test5 <- data.frame("LSM_ID"=0, "SMILES"=0)
#               for (i in 1:length(df_centroid$Representative)){    
#                 if (length(which(display$LSM_ID %in% df_centroid$Representative[[i]]))>0){
#                   test3 <- which(display$LSM_ID %in% df_centroid$Representative[[i]])
#                   for (j in 1:length(test3)){
#                     test4 <- display[test3[[j]], c("LSM_ID", "SMILES")]
#                     if (test5$LSM_ID[[1]]!=0){
#                     test5 <- rbind(test5, test4)
#                     }
#                     else{
#                       test5 <- test4
#                     }
#                   }
#                 }
#               }
#             }
#             colnames(test5) <- c("Representative", "SMILES")
#             test5 <<- test5
#               
#             
#             df_centroid <- merge(test5, df_centroid, stringsAsFactors = FALSE)
#             df_centroid <- data.frame(lapply(df_centroid, as.character), stringsAsFactors = FALSE)
#             df_centroid$Cluster <- as.integer(df_centroid$Cluster)
#             df_centroid <- df_centroid[order(df_centroid$Cluster),]
#             #df_centroid <- df_centroid[,-2]
#             for(i in 1:nrow(df_centroid)){
#               df_centroid$Compound[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == df_centroid$Representative[i]), "Compound"]
#             } 
#             
#             
#             df_centroid <- df_centroid[,-2]
#         
#             output$Representatives <- DT::renderDataTable(df_centroid)
#             df_centroid <<- df_centroid
#             # for (i in 1:length(df_centroid$Representative)){
#             #   if (!is.na((grep("LSM", df_centroid$Representative[[i]])) && (grep("LSM", df_centroid$Representative[[i]])==1))){
#             #     df_centroid$Representative[[i]] <- makeLink(df_centroid$Representative[[i]])
#             #   }
#             #   else if(!is.null(grep("LSM", df_centroid$Representative[[i]]))){
#             #     df_centroid$Representative[[i]] <- gsub("ZINC0", "ZINC", df_centroid$Representative[[i]])              
#             #     df_centroid$Representative[[i]] <- makeLink2(df_centroid$Representative[[i]])
#             #   }
#             # }
#             #  
#             #   return(df_centroid)
#             #   #df_centroid <- df_centroid[,-2]
#             #   #df_centroid <<- df_centroid
#             # }, escape = FALSE, rownames = FALSE)
#             #print(df_centroid$Representative)
#             
#             
#             print(paste("Got the centroids as ", centroid[[1]], sep=""))
#             incProgress(1/2, detail = "Plotting MDS")
#             
#             source("lib/MDS.R")
#           #7. Generate MDS plot from representatives distance to one another, pie chart radius corresponds to cluster size
#             MDS_plot(centroid)
# #            centroid_clusters(centroid, ClusterMembers)
#             #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
#             #                             + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
#             output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
#                                          + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "My Candidates")))
#             #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid), hjust=1, vjust=2)) 
#               #                           + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
#             incProgress(2/2, detail = "Complete!")                                                         
#             })
#         })
#             #p <- (ggplot() + geom_text_repel(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid))) 
#             #    + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added")))
#             
#             #move_layers(p, geom_text_repel(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid)), position = "top"))
#             
#             #output$MDSPlot<<- renderPlot(p)
#                                                                                                       
#             #                            })
#                                          
            
  #####Set it up so you can download Clusters the Representatives are part of
  
  output$RepDownload <- downloadHandler(
    filename = "representatives.csv",
    content = function(filename){
      write.csv(df_centroid, filename)
    }
  
)
  output$ClusterDownload <- downloadHandler(
    filename = "clusters.csv",
    content = function(filename){
      write.csv(ClusterMembers, filename)
    }
  )
  
  output$Max_Scores<- downloadHandler(
    filename = "max_scores.csv",
    content = function(filename){
      write.csv(combined_score_output, filename)
    }
    
  )  
# output Combined Scores
#  source("./lib/Combined_Score.R")
#  combined_score_output <<- combined_score(simMA, max_cons_2)
#  output$Combined_Score <- DT::renderDataTable({
#    
#    return(combined_score_output)
    #df_concordances_o$`LSM-ID`<- makeLink(df_concordances_o$`LSM-ID`)
    #return(display)
    #return(df_concordances_o)
#  },
# escape = FALSE, rownames = FALSE)
#  output$Similar_NCI <- eventReactive(input$GetNCI, {
  observeEvent(input$GetNCI, {
    print("Finding similar NCI compounds")
    source("./lib/SimilaritySearch.R")
    nci_list <- list()
    for(i in 1:length(centroid)){
      similarity_search(centroid[i])
      nci_list[[i]] <- output_nci}
    #nci_df <- as.data.frame(cbind(unlist(nci_list, recursive = FALSE), unlist(ClusterCentroid)))
    nci_df <- mapply(c,ClusterCentroid, nci_list)
    nci_df2 <- list()
    for(i in 1:length(nci_df)){
      if(is.character(nci_df[[i]])){ 
        nci_df2 <- list.append(nci_df2, nci_df[[i]])}
    }
    
    nci_df3 <- data.frame()
    for (j in 1:length(nci_df2)){
      for (k in 2:length(nci_df2[[j]])){
        cluster <- nci_df2[[j]][1]
        nsc <- nci_df2[[j]][k]
        cmpd_row <- cbind(nsc, cluster)
        nci_df3 <- rbind(nci_df3, cmpd_row)
      }
    } 
    #cluster_number <- unlist(ClusterCentroid)
    #nci_unlist <- unlist(nci_list)
    
    makeLink2 <- function(val) {
      #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
      paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
    }
    output$Similar_NCI <- DT::renderDataTable({
      
      #nci_df3$link <- makeLink2(nci_df3$nsc)
      nci_df3$nsc <- makeLink2(nci_df3$nsc)
      return(nci_df3)
      
    }, escape = FALSE)
    
    
    
#    output$Similar_NCI <- renderDataTable(nci_df3)
    #print(nci_list)})
  })

######Not yet working, df2 is not a global variable yet.
    observeEvent(input$GetSTITCH, {
    print("Getting Script")
    #source("./lib/StitchNetFn3.R")
    source("./lib/StitchNetFn7.R")
    print("sourced")
    #stitchNet("sirolimus", "mtor")
    SelectedCluster <- input$ClusterNumber
    h <- 1
    CompoundLSM <- list()
    for (i in 1:length(ClusterMembers)){
#      if (ClusterMembers[[i]] %in% LINCSCompounds)
      if (ClusterMembers[[i]] == SelectedCluster) {
        CompoundLSM[[h]] <- labels(ClusterMembers)[[1]][[i]]
        h <- h+1
      }
    }
    CompoundLSM <<- CompoundLSM
    print("IDs obtained")
    h <- 1

    CompoundName <- list()
    for (j in 1:length(CompoundLSM)){
      if (grep("LSM-", CompoundLSM[[j]])==1){
        #for (f in 1:length(df2$lincspertid)){
          if (CompoundLSM[[j]] %in% df2$lincspertid){
            index <- which(df2$lincspertid==CompoundLSM[[j]])
            if (length(index)>1){
              index <- index[[1]]
            }
            #CompoundName[[h]] <- df2$stitchID[[index]]
            CompoundName[[h]] <- df2$compound[[index]]
            h <- h+1
        #  }
        }
      }
      else if (!("LSM" %in% CompoundLSM[[j]])){
        CompoundName[[h]] <- CompoundLSM[[j]]
        h <- h+1
      }
      
    }
    CompoundName <<- CompoundName
    
    chemInput <- unlist(CompoundName)
    stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections)
    print("function ran")
    output$STITCHPlot <- renderVisNetwork(visNetwork(molnodes, interScore) %>% 
                                            visNodes(physics = FALSE) %>%
                                            visEdges(physics= FALSE, smooth = FALSE, color = list(color = "gray", highlight = "blue")) %>%
                                            visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
                                                      color = list(background = "ffa7a7", border = "black", 
                                                                   highlight = list(background = "ffa7a7", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Gene KD", shape = "circle",
                                                      color = list(background = "fdff00", border = "black",
                                                                   highlight = list(background = "feff8f", 
                                                                                    border = "red")) ) %>%
                                            visGroups(groupname = "STITCH connectors", shape = "dot", size = 10,
                                                      color = list(background = "ddd", border = "black", 
                                                                   highlight = list(background = "caff44", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                                      color = list(background = "4e91fc", border = "253085", 
                                                                   highlight = list(background = "4e91fc", 
                                                                                    border = "red"))) %>%
                                            visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = FALSE) %>%
                                            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                            visInteraction(navigationButtons = TRUE) %>%
                                            visLegend(width = 0.25) %>%
                                            visExport(type="png", name="export-network", float="left")
                                          )
    print("finished")

  })
########
    
  })
