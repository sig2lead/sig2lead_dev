###########################
# minLSH
########

minlsh <- function(fpset_added, b, block, lsm_rows, lincs_fps_2){
  
  load("./precomputed_objects/lsh_lincs_precomputed_64.Rdata")
  #load("./lincs_fps.RData")
  
  m=1024
  
  #lincs_fps <- lincs_fps_2
  #combined_fps <- rbind(fpset_antipsych, fps)
  #combined_fps <- rbind(fps)
  
  #####################################################################################
  ######## Bin query for each block ################################
  
  
  start_time <- Sys.time()
  query <- fpset_added
  query_blocks <- list()
  query_test <- list()
  for (i in 1:(m/b)){
    
    query_blocks[[i]] <- query[,((i*b)-(b-1)):(i*b)]
    
  }
  
  for (j in 1:(m/b)){
    query_blocks[[j]] <- as.data.frame(apply(query_blocks[[j]],1,paste0,collapse=""))
  }
  
  query_test <- vector("integer", length=nrow(query))
  start_time <- Sys.time()
  
  for (k in 1:nrow(query)){
    query_test[k] <- reps[[block]][which(reps[[block]]$hash == query_blocks[[block]][k,1]), "bin"]
  }
  
  ######################################################################################################
  # Compute Tanimoto with Query and Candidates in Same Bin
  ########################################################
  
  
  ################################################################################################################
  ######## Compute Tanimoto of compounds in query bin ###############
  test_fpsim_block <- list()
  #l <- 1
  
  candidates_df <- as.data.frame(candidates[[block]])
  for (l in 1:nrow(fpset_added)){
    lincs_candidates <- candidates_df[which(candidates_df$bin == query_test[[l]]),"id"]
    lincs_concordant <- rownames(lincs_fps_2)[lsm_rows]
    lincs_candidates_concordant <- lincs_candidates[which(lincs_candidates %in% lincs_concordant)]
    test_fpsim_block[[l]] <- (fpSim(query[l,], lincs_fps_2[which(rownames(lincs_fps_2) %in% lincs_candidates_concordant), ]))
    test_fpsim_block[[l]] <- test_fpsim_block[[l]][which.max(test_fpsim_block[[l]])]  
    
  }
  
  
  end_time <- Sys.time()
  return(test_fpsim_block)
  query_time <- end_time - start_time
}