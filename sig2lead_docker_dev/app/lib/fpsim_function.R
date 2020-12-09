fpsim <- function(fpset_added, lsm_rows){
  load("./lincs_fps.RData")
  #fpSim_start_time <- Sys.time()
  simmA_fpSim <-sapply(rownames(fpset_added), function(x) max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  lincs_fpsim <-sapply(rownames(fpset_added), function(x) which.max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #lincs_fpsim_2 <- lincs_fps_2[lincs_fps_2]
  #fpSim_end_time <- Sys.time()
  #total_fpSim_time[[j]] <- fpSim_end_time - fpSim_start_time
  #print(paste("Total time to run fpSim: ", total_fpSim_time[[j]], sep=""))
  
  lincs_fp_ind <- lsm_rows[lincs_fpsim]
  lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  simmA_fpSim_df <- data.frame(lincs_vec,simmA_fpSim)
  return(simmA_fpSim_df)
  
}