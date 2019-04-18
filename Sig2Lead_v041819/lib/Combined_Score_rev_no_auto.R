combined_score <- function(){
  
  #######Extract Vectors of Similarites for Added Compounds from Similarity Matrix
  
      #adds_csv
  
      if(!(grepl(".sdf", adds$name))){
        added_sims <- simMA[(rownames(simMA) %in% adds_csv[,2]),]
      }
      else if(grepl(".sdf", adds$name)){
        added_ids <<- unlist(lapply(1:length(sdfset_add),function(x) sdfset_add@SDF[[x]]@header[1]))
        added_sims <<- simMA[(rownames(simMA) %in% added_ids),]
      }
      
  
  #######Extract Concordances for LSM-IDs in Similarity Matrix
      
      
  
  ########Calculate Concordance Matrix###############
  # Added Compounds X LINCS Compounds
  ##################################################
      
      # added_sims for LSM x Con for LSM
      
      #combined_score <-  max_cons_2[which(max_cons_2$`LSM-ID`%in% colnames(added_sims)),]
      
      combined_score <-  max_cons_4[which(max_cons_4$`LSM-ID`%in% colnames(added_sims)),]
    
      added_sims_2 <- added_sims[,(colnames(added_sims) %in% combined_score[,1])]
      
      added_sims_o <<- added_sims_2[,order(colnames(added_sims_2))]
      print(paste("Dimensions of added_sims_o:", dim(added_sims_o), sep= " "))
            
      
  ####################################
  # Order Concordances
  ####################################
      
      combined_score_o <- combined_score[order(as.character(combined_score[,1])),]
      concordance <- combined_score_o
      
      
  ################################
  # Multiply concordance x Similarity
  ################################    
      
  score_df <-  data.frame(matrix(0, nrow=nrow(added_sims_o), ncol=ncol(added_sims_o)))
  rownames(score_df) <- rownames(added_sims_o)
  colnames(score_df) <- colnames(added_sims_o)
  print("Made it this far...")
  for (i in 1:nrow(added_sims_o)){    
      for (j in 1:length(concordance[,1])){
      score_df[i,j] <- added_sims_o[i,j] * concordance[j,3]
      }
    print(paste("Compound Score Interation: ", i, sep=""))
  }
  print(" Made it this far too..")
  #####################
  combinations  <<- as.vector(as.matrix(score_df)) 
  added_sims_vec <<- as.vector(as.matrix(added_sims_o))
  concordance_vec <<- rep(concordance[,3], each=(length(combinations)/nrow(concordance)))
  #added_sims_vec <- as.vector(as.matrix(added_sims_o))
  
  #times_row <- nrow
  rows_rep <- rep(rownames(score_df), times=(length(combinations)/nrow(score_df)))
  cols_rep <- rep(colnames(score_df), each = length(combinations)/ncol(score_df))

   
  ranked_scores <- score_df
   #ranked_scores <- data.frame(rank(score_df)) 
   #return(score_df)
  scores_combined <- as.data.frame(cbind(rows_rep, cols_rep, concordance_vec, added_sims_vec, combinations)) 
  #scores_combined <- as.data.frame(cbind(rows_rep, cols_rep, combinations))
  scores_combined$combinations <- as.numeric(as.character(scores_combined$combinations))
  ranked_scores <- scores_combined[order(scores_combined[,5], decreasing = TRUE),]
  colnames(ranked_scores) <- c("My Compounds", "LINCS Compounds", "Concordance", "Similarity", "Combined Score")
  ranked_scores_2 <<- ranked_scores[!duplicated(ranked_scores$'My Compounds'),]
  end_time <<- Sys.time()
return(ranked_scores_2) 
  
}

