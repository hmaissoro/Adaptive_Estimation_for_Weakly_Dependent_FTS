################################################################################
# File:             03_estimate_mean.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     23.01.2024
# 
# Contains R-script for mean estimation for numerical study: FTS Model 2.
# 
# For a detailed description see:
#   Maissoro, H., Patilea, V. and Vimond, M. (2024). Adaptive estimation for 
#     weakly dependent functional time series.
################################################################################

# Load packages and functions
library(data.table)
library(matrixStats)

all_func <- list.files("../R/")
lapply(all_func, function(func){
  if (func != "zzz.R") source(file = file.path("../R/", func), echo = FALSE)
})

# Get the script path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# Mean function estimation function ----
t0 <- c(0.2, 0.4, 0.7, 0.8)
## Estimate mean functions risk
estim_mean_risk_fun <- function(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", t0 = t0, n_cores = 30){
  
  # Load local regularity estimates
  locreg_file_name <- paste0(script_path, "/FAR/locreg_estimates/dt_locreg_FAR_",
                             white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(locreg_file_name)
  
  # Define bandwidth grid
  K <- 40
  b0 <- 4 * (N * lambda) ** (- 0.9)
  bK <- 5 / 10
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()
  
  if (white_noise == "mfBm") {
    # Estimate mean risk by mc
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda))
    index_mc <- gsub(pattern = paste0("dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                     replacement = "", x = data_file_list)
    index_mc <- as.numeric(index_mc)
    
    res <- parallel::mclapply(index_mc, function(mc_i, dt_locreg, Ni, lambdai, bw_grid, t0){
      # Load data
      data_file_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", Ni, "lambda", lambdai,
                               "/dt_mc_FAR_mfBm_N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt_random_mc <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
      
      # Estimate the risk of the mean function
      bw <- unique(dt_random_mc[, .(id_curve, presmooth_bw)])[, presmooth_bw]
      Ht <- dt_locreg[id_mc == mc_i & order(t), Ht]
      Lt <- dt_locreg[id_mc == mc_i & order(t), Lt]
      
      dt_risk_muhat <- estimate_mean_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = 'tobs', ycol = "X",
        t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
      )
      
      # Save
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_risk_muhat)
      rm(dt_risk_muhat) ; gc()
      
      file_save <- paste0(script_path, "/FAR/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda,
                          "/dt_mean_risk_FAR_mfBm_N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      saveRDS(object = dt_res, file = file_save)
      rm(dt_res) ; gc()
      return(paste0("Mean risk : done for MC = ", mc_i, ", N = ", N, " and lambda = ", lambda))
    }, mc.cores = n_cores, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0)
    
  } else if (white_noise == "fBm") {
    # Not need for the paper ...
  }
  
  rm(locreg_file_name, dt_locreg) ; gc()
  return(paste0("Mean risk : done for FTS Model ", gsub("d", "", design), ", N = ", N, " and lambda = ", lambda, " at ", Sys.time()))
}

## Estimate mean function only for "plus_mean"
estim_mean_fun <- function(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", t0, n_cores = 30){
  
  if (white_noise == "mfBm") {
    # Get the MC index in the mean risk file name 
    data_risk_file_list <- list.files(paste0(script_path, "/FAR/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda))
    id_mc_data_risk <- gsub(pattern = paste0("dt_mean_risk_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                            replacement = "", x = data_risk_file_list)
    id_mc_data_risk <- as.numeric(id_mc_data_risk)
    
    # Get the MC index in the data file name 
    data_file_list <- list.files(paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)
    
    # Get the already done MC repetition
    index_mc <- id_mc_data[which(id_mc_data_risk %in% id_mc_data)]
    
    # Estimate mean by mc
    dt_mean_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, Ni, lambdai, t0){
      # Get optimal bandwidth smoothing parameter
      data_risk_file <- paste0(script_path, "/FAR/mean_estimates/", design, "_mean_risk/N", Ni, "lambda", lambdai, 
                               "/dt_mean_risk_FAR_mfBm_N=", Ni, "_lambda=", lambdai, "_id_mc=", mc_i, "_", design, ".RDS")
      dt_mean_risk <- readRDS(data_risk_file)
      dt_optbw <- dt_mean_risk[, .("optbw" = h[which.min(mean_risk)]), by = c("id_mc", "t")]
      
      # Load data
      data_file_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", Ni, "lambda", lambdai,
                               "/dt_mc_FAR_mfBm_N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      
      # Extract data
      dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
      optbw <- dt_optbw[id_mc == mc_i][order(t), optbw]
      t0 <- dt_optbw[id_mc == mc_i][order(t), t]
      
      # Estimate the mean function
      dt_mean <- estimate_mean(
        data = dt_random_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, optbw = optbw)
      
      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_mean[, .(t, optbw, PN, muhat)])
      dt_res <- data.table::merge.data.table(
        x = dt_res,
        y = unique(dt[ttag == "tcommon", .("t" = tobs, "mutrue" = process_mean)]),
        by = "t")
      rm(optbw, dt_mean, dt_mean_risk) ; gc()
      return(dt_res)
    }, mc.cores = n_cores, Ni = N, lambdai = lambda))
    
  } else if (white_noise == "fBm") {
    # Not need for the paper
  }
  
  ## Save
  file_name <- paste0(script_path, "/FAR/mean_estimates/dt_mean_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(file_name, dt_mean_mc) ; gc()
  
  return(paste0("Done : dt_mean_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate mean : FTS Model 2 ----
## Estimate mean risk function
estim_mean_risk_fun(N = 150, lambda = 40, white_noise = "mfBm", design = "d2", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 1000, lambda = 40, white_noise = "mfBm", design = "d2", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 1000, lambda = 1000, white_noise = "mfBm", design = "d2", t0 = t0, n_cores = 30)

## Estimate mean function
estim_mean_fun(N = 150, lambda = 40, white_noise = "mfBm", design = "d2", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 40, white_noise = "mfBm", design = "d2", n_cores = 30)
estim_mean_fun(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 1000, white_noise = "mfBm", design = "d2", n_cores = 30)

