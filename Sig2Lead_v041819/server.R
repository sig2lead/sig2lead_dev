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

## add this #########
library(readr)
load("./lib/KDgenes") 
######################

#source("http://bioconductor.org/biocLite.R")

#biocLite("ChemmineOB")
#biocLite("ChemmineR")


######Discordant or Concordant
######Benchmarking -- Behrouz

begin_time <<- Sys.time()
options(shiny.maxRequestSize=200*1024^2)

shinyServer(function(input, output) {
  
  output$mytest <- renderUI(
    if(!is.null(input$AddedCompounds)) {
      div(textInput("AddedLabel", "Label for added compounds", value="Added"),
          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
      )
    }
  )  
 observeEvent(input$Go,
#  output$Go <- eventReactive(input$Go,  
    {
      withProgress(message = "Running Sig2Lead", value = 0, {
        source("lib/GeneKD.R", local = TRUE)
        source("lib/ConDisconSignatures.R", local = TRUE)
      #1. Query iLINCS for user input gene

#### add this ####################################################################################################                   
        
        if ( (input$Signature == "Input a Gene") & !(toupper(input$gene_knockdown) %in% toupper(KDgenes)) ){
          incProgress(3/7, message = "Not found the queried gene knockdown")    
          source("./lib/KDnet.R")
          gene_knockdown <<- toupper(input$gene_knockdown)
          incProgress(5/7, message = "Retrieve STITCH interactors") 
          #KDnet(gene_knockdown)
          output$notify1 <- renderUI({
            HTML(paste(c("<h3> <b> <font color =", '"blue"', ">", toupper(gene_knockdown), " </font> </b> knockdown <b> NOT </b> present in LINCS !!! </h3>",
                         "<h4> The known interactors of", toupper(gene_knockdown), "are retrieved from the STITCH database. </h4>")
            ))
          })               
          
          output$KDnet_Plot <- renderVisNetwork(KDnet(gene_knockdown))

          output$notify2 <- renderUI({
            HTML(paste(c("<h4> Please try again using <b> <font color =", '"red"', "> Genes </font> </b> in the <b> <font color =", '"green"', "> Green Nodes </font> </b> (as knockdowns presented in LINCS). </h4>")
            ))
          })          
          incProgress(7/7, message = "Finish") 
        } else if( (input$Signature == "Input a Gene" & toupper(input$gene_knockdown) %in% toupper(KDgenes)) | (input$Signature == "Upload a Signature") ){
          
#######################################################################################################                   
          
        if (input$Signature == "Input a Gene"){
          define_knockdown(input$gene_knockdown)
        print("Found your knockdowns") 
        incProgress(1/7, message = "Finding Concordant Signatures")
        ConDisconSignatures(sigid)
        #ConDisconSignatures <<- ConDisconSignatures(sigid)
        }
        else if (input$Signature == "Upload a Signature"){
          source("lib/UploadSignature.R")
          uploadsig <- input$UploadSignature
          sigUpload(paste(getwd(), uploadsig$name, sep="/"))
          print("Found your signature")
          incProgress(1/7, message = "Finding Concordant Signatures")
          UploadConDiscon(rjson)
          
        }
        
        #incProgress(1/7, detail = "Finding Concordant Signatures")
        
      #2. Get a list of all concordant compounds
        
        
        print("Found concordant signatures")
        incProgress(1/7, message = "Obtaining SMILES")
        source("lib/GetSMILES.R", local = TRUE)
      #3. Generate SMILES of compounds
        GetSMILES()
        print("SMILES acquired")
        incProgress(1/7, message = "Converting SMILES to SDF...Slow")
      #4. Generate Data Table
        display <<- final_compounds
        colnames(display) <- c("LSM_ID", "SMILES")
        
        #output$SMILES <<- DT::renderDataTable(display)
        #output$SMILES <<- DT::renderDataTable({DT::datatable(display, selection = list(selection = "single", target = "cell"))}, escape = FALSE)
        
        makeLink <- function(val) {
          #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
        }
        
        makeLink2 <- function(val) {
          #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
        }
        #lincs <<-  data.frame(matrix(unlist(con_comma[[1]]), nrow= 457, byrow=T), stringsAsFactors = FALSE)
        #lincs <<- con[order(con[,2], decreasing = TRUE),c(11,6,2)]
        #colnames(lincs) <- c("LSM-ID", "Compound", "Concordance")
        
        makeLink3 <- function(val) {
          #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
        }
        
        makeLink4 <- function(val) {
          #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
          paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
        }
        
        df <- list()
        for (d in 1:length(con_comma)){
          
          df[[d]] <- as.data.frame(matrix(unlist(con_comma[[d]]), nrow=length(con_comma[[d]][[1]])))
          
        }
        df_concordances <<- ldply(df, data.frame)
        # Change 4-16-2019
        df_concordances_2 <- df_concordances[,c(11,6,2,3)]
        colnames(df_concordances_2) <- c("LSM-ID", "Compound", "Concordance", "Significance")
        #df_concordances_2 <- df_concordances[,c(11,6,2)]
        #colnames(df_concordances_2) <- c("LSM-ID", "Compound", "Concordance")
        #df_concordances_2$Concordance <- round(as.numeric(as.character(df_concordances_2$Concordance)), 3)
  #################      
        
        
        df_concordances_2$Concordance <- round(as.numeric(as.character(df_concordances_2$Concordance)),3)
        df_concordances_o <<- df_concordances_2[order(df_concordances_2$Concordance, decreasing = TRUE ),]
        
        
        df_concordances_no_dups <<- df_concordances_o[!duplicated(df_concordances_o["LSM-ID"]),]
        display_no_dups <- display[!duplicated(display["LSM_ID"]),]
        
        df_concordances_o_smiles <- merge(df_concordances_no_dups, display_no_dups, by.x="LSM-ID", by.y="LSM_ID", all.x=FALSE, all.y=FALSE)
        df_concordances_o_smiles<<- df_concordances_o_smiles[order(df_concordances_o_smiles$Concordance, decreasing = TRUE ), c(1,5,2,3,4)]
        
        #############Output to LINCS Compounds Table;
        
        #################
        # Max Concordance
        #################
        max_cons <- df_concordances_o_smiles[order(df_concordances_o$`LSM-ID`, -abs(df_concordances_o$Concordance)), ]
        max_cons_2 <- max_cons[!duplicated(max_cons$`LSM-ID`),]
        max_cons_3 <<- max_cons_2[order(max_cons_2$Concordance, decreasing = TRUE), -5]
        max_cons_4 <<- max_cons_3[-which(is.na(max_cons_3$`LSM-ID`)),]
        #con_df <- as.data.frame((matrix(unlist(con_comma),byrow=T)))
        
        
        output$SMILES <- DT::renderDataTable({
          
          #display$link <- makeLink(display$LSM_ID)
          #display$LSM_ID <- makeLink(display$LSM_ID)
          #df_concordances_o_smiles$`LSM-ID`<- makeLink(df_concordances_o_smiles$`LSM-ID`)
          #max_cons_3$`LSM-ID`<- makeLink(max_cons_3$`LSM-ID`)
          #return(display)
          #return(display)
          #return(df_concordances_o_smiles)
          return(max_cons_3)
        },
        escape = FALSE, rownames = FALSE)
        
        #http://lincsportal.ccs.miami.edu/SmallMolecules/view/LSM-5467
        display <<- display
        
        source("lib/ChemmineOB.R", local = TRUE)
        
        ###This is the slow step, see if we can go straight from smiles to fingerprints and seperately convert to SDF
        print("Converting SMILES to SDF")
        Get_SDF()
        
        ###
        
        print("SDF conversion complete")
        incProgress(3/7, message = "Clustering Related Compounds")
      #5. Cluster compounds identified through LINCS, generate heatmaps
        # If Bypass Clu
        source("lib/chemmineR.R", local = TRUE)
        
#    if(input$Bypass == FALSE){
    if(!is.null(input$AddedCompounds)){
      adds<<-input$AddedCompounds
      #5.1. Include Added Compounds if applicable
         #if (!is.null(input$AddedCompounds))      {
         if (input$Bypass == FALSE)      {                     
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
                    sdfset_add <<- read.SDFset(adds$datapath)
                    #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
                    #sdf_smiles <<- c(sdf_smiles, sdfset_add[1:10])
                    bypass_clustering(sdf_smiles,sdfset_add)
                  # cluster_compounds()
                  # print("LINCS clustering complete!")
                  # output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
                }
    }#)
    }
#    else{print("Clustering bypassed")
    else{print("No added compounds, must cluster.")
      #if (!is.null(input$AddedCompounds))      {
        adds <<- NULL
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
        source("./lib/Combined_Score_auto_2.R")
        source("./lib/Combined_Score_auto_3.R")
        if (!is.null(input$AddedCompounds)){
          if(input$Bypass == TRUE){   
              combined_score_output <- combined_score_2()
              output$Combined_Score <- DT::renderDataTable({
                
                combined_score_output$Concordance <- round(as.numeric(as.character(combined_score_output$Concordance)),3)
                combined_score_output$Similarity <- round(as.numeric(as.character(combined_score_output$Similarity)),3)
                combined_score_output$'Combined Score' <- round(as.numeric(as.character(combined_score_output$'Combined Score')),3)
                
                
              
                
                combined_score_output$`LSM-ID`<- makeLink(combined_score_output$`LSM-ID`)
              combined_score_output$`My Compounds`<- makeLink2(combined_score_output$`My Compounds`)
          
              return(combined_score_output)
              #df_concordances_o$`LSM-ID`<- makeLink(df_concordances_o$`LSM-ID`)
              #return(display)
              #return(df_concordances_o)
              
        },
        escape = FALSE, rownames = FALSE)  
          }else{
            print("Combined Score with Clustering")
            combined_score_output <- combined_score_3()
            output$Combined_Score <- DT::renderDataTable({
            
              combined_score_output$Concordance <- round(as.numeric(as.character(combined_score_output$Concordance)),3)
              combined_score_output$Similarity <- round(as.numeric(as.character(combined_score_output$Similarity)),3)
              combined_score_output$'Combined Score' <- round(as.numeric(as.character(combined_score_output$'Combined Score')),3)
              
              combined_score_output$`LSM-ID`<- makeLink3(combined_score_output$`LSM-ID`)
              combined_score_output$`My Compounds`<- makeLink4(combined_score_output$`My Compounds`)
              #combined_score_output$)
              return(combined_score_output)}, escape=FALSE, rownames=FALSE)
          }} 
        else{
          print("No added compounds so can't calculate combined score")
        }
        
#        cluster_compounds()
        
#        output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
       incProgress(7/7, message = "Complete!")

############ add this #####################################                  
        }
######################   
        
    })
    }) # End input$Go
 
  output$SMILESDownload <- downloadHandler(
    filename = function() {
      paste(input$gene_knockdown, '_compounds.csv', sep='')
    },
    content = function(filename){
      write.csv(final_compounds, filename)
    }
    )
  
  output$SDFDownload <- downloadHandler(
    filename = function() {
      paste(input$gene_knockdown, '.sdf', sep='')
    },
    content = function(filename){
      write.SDF(sdf_smiles, filename)
    }
  )
  
  output$ConTableDownload <- downloadHandler(
    
    filename = "concordance_table.csv",
    content = function(filename){
      write.csv(df_concordances_o_smiles, filename)
    }
  )
#  output$Representatives <- eventReactive(input$GetRepresentatives,
  observeEvent(input$GetRepresentatives,
        {
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
            df_centroid<<-df_centroid
        
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
              df_centroid <<- df_centroid
            }, escape = FALSE, rownames = FALSE)
            #print(df_centroid$Representative)
            
            
            print(paste("Got the centroids as ", centroid[[1]], sep=""))
            incProgress(1/2, detail = "Plotting MDS")
            
            source("lib/MDS.R")
          #7. Generate MDS plot from representatives distance to one another, pie chart radius corresponds to cluster size
            MDS_plot(centroid)
#            centroid_clusters(centroid, ClusterMembers)
            output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
                                         + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS")))
            if(!is.null(input$AddedCompounds)){
              output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2))
                                         + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", Added_Label)))
              incProgress(2/2, detail = "Complete!")
            }
            else{
              output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2))
                                           + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added")))
              incProgress(2/2, detail = "Complete!")
            }
            #}
            #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid), hjust=1, vjust=2)) 
              #                           + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
            #incProgress(2/2, detail = "Complete!")                                                         
            })
        })
            #p <- (ggplot() + geom_text_repel(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid))) 
            #    + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added")))
            
            #move_layers(p, geom_text_repel(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid)), position = "top"))
            
            #output$MDSPlot<<- renderPlot(p)
                                                                                                      
            #                            })
                                         
            
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
    source("./lib/StitchNetFn9.R")
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
     # if (grep("LSM-", CompoundLSM[[j]])==1){

## add this ###
      if (length(grep("LSM-", CompoundLSM[[j]]))!=0 && grep("LSM-", CompoundLSM[[j]])==1){
###############
        
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

######## add this script to server ######
    if(length(CompoundName) == 0){
      CompoundName <- CompoundLSM
    }
#########################################
    
    CompoundName <<- CompoundName
    
    chemInput <- unlist(CompoundName)
    #stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections)

## add this ########################        
    chemInput <<- chemInput
    stitchNet(chem = chemInput, gene = input$gene_knockdown, limit = input$Connections, confidence = input$confidence)
#################################
       
    print("function ran")

##### add/replace this script to server #######################################    
    observeEvent(input$show_undetectNodes, {
      if(!input$show_undetectNodes){
        select_molnodes <- molnodes %>% filter(group != "Undetected by STITCH")
        assign("select_molnodes", select_molnodes, envir = .GlobalEnv)
      } else if(input$show_undetectNodes){
        select_molnodes <- molnodes 
        assign("select_molnodes", select_molnodes, envir = .GlobalEnv)
      }
    })
    
    # output$STITCHPlot <- renderVisNetwork(stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections))
    
    output$STITCHPlot <- renderVisNetwork(visNetwork(select_molnodes, interScore) %>%
                                            visNodes(physics = !(input$fixedNet)) %>%
                                            visEdges(physics = !(input$fixedNet), smooth = TRUE, color = list(color = "#008080", highlight = "blue")) %>% #"grey
                                            visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
                                                      font = list(color = "660000"), #"133f08"
                                                      color = list(background = "ffa7a7", border = "black", #"aefc6a" 
                                                                   highlight = list(background = "ffa7a7", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Gene KD", shape = "circle",
                                                      font =  list(color = "blue"),
                                                      color = list(background = "fdff00", border = "black",
                                                                   highlight = list(background = "feff8f", 
                                                                                    border = "red")) ) %>%
                                            visGroups(groupname = "STITCH gene/compound", shape = "dot", size = 10,
                                                      font = list(color = "magenta", size = 10),
                                                      color = list(background = "ddd" , border = "black", #background ="#bb84ff" 
                                                                   highlight = list(background = "caff44", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                                      color = list(background = "82a9ff", border = "253085",  #4e91fc
                                                                   highlight = list(background = "4e91fc", 
                                                                                    border = "red"))) %>%
                                            visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = !(input$fixedNet), smooth = TRUE) %>% #"layout_nicely"
                                            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                            visInteraction(navigationButtons = TRUE) %>%
                                            visLegend(width = 0.25) #%>%
                                          #visExport only works under Shiny environment, so enable a code below when running in Shiny environment
                                          #visExport(type="pdf", name="export-network", label = "Save network", float="left", background = "#fff")
    )
    ###########################################################################    
    
    
    
#    output$STITCHPlot <- renderVisNetwork(visNetwork(molnodes, interScore) %>% 
 #                                           visNodes(physics = FALSE) %>%
  #                                          visEdges(physics= FALSE, smooth = FALSE, color = list(color = "gray", highlight = "blue")) %>%
   #                                         visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
    #                                                  color = list(background = "ffa7a7", border = "black", 
     #                                                              highlight = list(background = "ffa7a7", 
      #                                                                              border = "red"))) %>%
       #                                     visGroups(groupname = "Gene KD", shape = "circle",
        #                                              color = list(background = "fdff00", border = "black",
         #                                                          highlight = list(background = "feff8f", 
          #                                                                          border = "red")) ) %>%
           #                                 visGroups(groupname = "STITCH connectors", shape = "dot", size = 10,
            #                                          color = list(background = "ddd", border = "black", 
             #                                                      highlight = list(background = "caff44", 
              #                                                                      border = "red"))) %>%
               #                             visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                #                                      color = list(background = "4e91fc", border = "253085", 
                 #                                                  highlight = list(background = "4e91fc", 
                  #                                                                  border = "red"))) %>%
                   #                         visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = FALSE) %>%
                    #                        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                     #                       visInteraction(navigationButtons = TRUE) %>%
                      #                      visLegend(width = 0.25) %>%
                        #                    visExport(type="png", name="export-network", float="left")
                       #                   )
  #  print("finished")

  })
########

### add this ##########    
    
    observeEvent(input$GlobalSTITCH, {
      print("Starting mSTITCH")
      source("./lib/mStitchFn.R")
      withProgress(message = "Starting mSTITCH", value = 0, {
        Sys.sleep(0.25)
        
        S2L <- data.frame('lincspertid' = rownames(ClusterMembers), 
                          'clusterID' = ClusterMembers[ ,1], stringsAsFactors = FALSE)
        key <- df2
        lincsID <- key %>% select('compound', 'lincspertid') %>% unique()
        tmp <- left_join(S2L, lincsID, by = "lincspertid")
        ind <- which(is.na(tmp$compound), arr.ind = TRUE)
        tmp$compound[ind] <- tmp$lincspertid[ind]
        
        clusterList <- list()
        for(i in 1:length(unique(tmp$clusterID))){
          clusterList[i] <- tmp %>% filter(clusterID == i) %>% select(compound)
        }
        
        numbers <- as.numeric()
        for (i in 1:length(clusterList)){
          numbers[i] <- length(clusterList[[i]])
        }
        
        clus_size <- as.numeric()
        if(!(input$all_clusters)){
          clus_size <- 3
        } else {
          clus_size <- 1
        }
        
        ind <- as.data.frame(numbers) %>% dplyr::filter(numbers >= clus_size) %>% nrow()
        clusterList_tom <- clusterList[1:ind]
        
        chemX <- unlist(clusterList_tom)
        genes <- input$gene_knockdown
        
        incProgress(0.6, message = "Gethering global gene-chemical interactions...Slow")
        
        stitch_all <- lapply(clusterList_tom, mStitch, gene = genes)
        doc <- bind_rows(stitch_all)  
        
        tmp <- strsplit(doc$X15, "|", fixed = TRUE) %>% sapply("[", 1) %>% 
          strsplit(":", fixed = TRUE) %>% sapply("[", 2) %>% 
          as.numeric()
        interScore <- data.frame("from" = doc$X3, "to" = doc$X4, "score" = tmp, "width" = (tmp*2)^2, 
                                 "title" = paste0("<p>", "Score:", "<br>", tmp, "</p>"), stringsAsFactors = FALSE)
        #interScore <<- interScore
        
        mol <- unique(c(interScore$from, interScore$to)) 
        group <- toupper(mol) %>% replace(. %in% toupper(chemX), "Input chemicals present in STITCH") %>%
          replace(. %in% toupper(genes), "Gene KD") %>%
          replace(!. %in% c("Input chemicals present in STITCH", "Gene KD"), "STITCH gene/compound")
        
        link_modnodes <- character()
        for (i in seq_along(mol)) for (j in seq_along(chemX)) {
          if(grepl(toupper(genes), toupper(mol[i])) ){
            link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=', mol[i], '>', mol[i], '</a>')
          } else if(grepl("CHEMBL", mol[i]) ){
            link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.ebi.ac.uk/chembl/beta/compound_report_card/', mol[i], '>', mol[i], '</a>')
          } else if(grepl("ZINC", mol[i]) ){
            link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'http://zinc15.docking.org/substances/',  mol[i], '>', mol[i], '</a>')
          } else if(grepl("LSM", mol[i]) ){
            link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', mol[i], '>', mol[i], '</a>')
          } else if(toupper(mol[i]) %in% toupper(chemX[j]) ){
            link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://pubchem.ncbi.nlm.nih.gov/search/#query=', mol[i], '>', mol[i], '</a>')
          } else if( is.na(link_modnodes[i]) ){
            link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://en.wikipedia.org/wiki/', mol[i], '>', mol[i], '</a>')
          }
        }
        
        molnodes <- data.frame(label = mol, id = mol, group = group, title = link_modnodes)
        
        # get stitch undetected compounds (caused by whatever reasons i.e., name/id confusing between chemInput and stitch) and their landing pages 
        ind <- toupper(chemX) %in% toupper(mol)
        stitch_undetect <- chemX[!ind]
        stitch_undetect <- if(as.numeric(length(stitch_undetect))==0){"NULL"} else{stitch_undetect} 
        
        # create landing pages for undetect_Node
        link_undetect_Node <- character()
        for (i in seq_along(stitch_undetect)){
          if(grepl("CHEMBL", stitch_undetect[i]) == 1){
            link_undetect_Node[i] <- paste0('<a target="_blank"href=https://www.ebi.ac.uk/chembl/beta/compound_report_card/', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
          } else if(grepl("ZINC", stitch_undetect[i]) == 1){
            link_undetect_Node[i] <- paste0('<a target="_blank"href=http://zinc15.docking.org/substances/', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
          } else if(grepl("LSM", stitch_undetect[i]) == 1){
            link_undetect_Node[i] <- paste0('<a target="_blank"href=http://lincsportal.ccs.miami.edu/SmallMolecules/view/', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
          } else {
            link_undetect_Node[i] <- paste0('<a target="_blank"href=https://pubchem.ncbi.nlm.nih.gov/search/#query=', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
          }
        }
        
        undetect_Node <- data.frame(label = stitch_undetect, id = stitch_undetect, group = "Undetected by STITCH", title = link_undetect_Node) 
        
        # when need unmatched inputs posted in the network, adding undetected compounds into molnodes
        molnodes <- rbind(molnodes, undetect_Node)
        molnodes <- unique(molnodes)
        # molnodes <<- molnodes
        
        incProgress(0.2, message = "Generating global network")
        Sys.sleep(0.5)
        
        output$numberClust <- renderText(paste0("Global STITCH network plotted from ", length(clusterList_tom), 
                                                " clusters (>=", clus_size, " compounds/cluster) (",  ceiling((100*length(clusterList_tom)/length(clusterList))), 
                                                "%)"))
        
        
        output$mSTITCHPlot <- renderVisNetwork(visNetwork(molnodes, interScore) %>%
                                                 visNodes(physics = TRUE) %>%
                                                 visEdges(physics= TRUE, smooth = TRUE, color = list(color = "#008080", highlight = "blue")) %>% 
                                                 visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
                                                           font = list(color = "660000"), 
                                                           color = list(background = "ffa7a7", border = "black", 
                                                                        highlight = list(background = "ffa7a7", 
                                                                                         border = "red"))) %>%
                                                 visGroups(groupname = "Gene KD", shape = "circle",
                                                           font =  list(color = "blue"),
                                                           color = list(background = "fdff00", border = "black",
                                                                        highlight = list(background = "feff8f", 
                                                                                         border = "red")) ) %>%
                                                 visGroups(groupname = "STITCH gene/compound", shape = "dot", size = 10,
                                                           font = list(color = "magenta", size = 10),
                                                           color = list(background = "ddd" , border = "black", 
                                                                        highlight = list(background = "caff44", 
                                                                                         border = "red"))) %>%
                                                 visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                                           color = list(background = "82a9ff", border = "253085",  
                                                                        highlight = list(background = "4e91fc", 
                                                                                         border = "red"))) %>%
                                                 visIgraphLayout(randomSeed = 9, layout = "layout_components", physics = FALSE, smooth = FALSE) %>% 
                                                 visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                                 visInteraction(navigationButtons = TRUE) %>%
                                                 visLegend(width = 0.25) %>%
                                                 visExport(type="png", name="export-network", float="left") # note: visExport only works under Shiny environment
        )
        incProgress(0.1, message = "Finish!")
        Sys.sleep(0.5)
        
      })
      
      
    })    
    
 #############################   
        
  })
