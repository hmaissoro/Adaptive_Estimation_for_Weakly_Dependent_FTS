################################################################################
# File:             06_estimate_autocov.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     22.01.2024
# 
# Contains R-script for autocovariance estimation for numerical study: 
#                 FTS Model 2 with zero-mean.
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

# Autocovariance estimation function ----
## Parameters
t0 <- c(0.2, 0.4, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]

rm(dt_st) ; gc()

## Autocovariance estimation function 
estim_autocov_risk_fun <- function(
  N = 150, lambda = 40, white_noise = "mfBm", design = "d4",
  s0 = s0, t0 = t0, lag = 1, n_cores = 30){
  
  # Load local regularity estimates
  locreg_file_name <- paste0(script_path, "/FAR/locreg_estimates/dt_locreg_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(locreg_file_name)
  
  # Define bandwidth grid
  K <- 40
  b0 <- 4 * max((N * lambda) ** (- 0.9), (N * (lambda ** 2)) ** (- 0.9))
  bK <- 5 / 10
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()
  
  if (white_noise == "mfBm") {
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda))
    index_mc <- gsub(pattern = paste0("dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                     replacement = "", x = data_file_list)
    index_mc <- as.numeric(index_mc)
    index_mc <- sort(index_mc)
    
    # Estimate the aucovariance risk by mc
    res <- parallel::mclapply(index_mc, function(mc_i, dt_locreg, N, lambda, white_noise, design){
      # Load data
      data_file_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda,
                               "/dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt_random_mc <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
      
      # Estimate the risk of the gamma function ----
      ## Extract local regularity parameters
      dt_locreg_mc <- dt_locreg[id_mc == mc_i]
      dt_locreg_mc <- data.table::merge.data.table(
        x = data.table::data.table("s" = s0, "t" = t0),
        y = dt_locreg[id_mc == mc_i, .("t" = t, Ht, Lt)],
        by = "t"
      )
      dt_locreg_mc <- data.table::merge.data.table(
        x = dt_locreg_mc,
        y = dt_locreg[id_mc == mc_i, .("s" = t, "Hs" = Ht, "Ls" = Lt)],
        by = "s"
      )
      
      bw <- unique(dt_random_mc[id_mc == mc_i, .(id_curve, presmooth_bw)])[, presmooth_bw]
      Ht <- dt_locreg_mc[order(s, t), Ht]
      Lt <- dt_locreg_mc[order(s, t), Lt]
      Hs <- dt_locreg_mc[order(s, t), Hs]
      Ls <- dt_locreg_mc[order(s, t), Ls]
      
      ## Estimate the risk function
      dt_risk_gammahat <- estimate_autocov_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = s0, t = t0, lag = 1, bw_grid = bw_grid,
        smooth_ker = epanechnikov, Hs = Hs, Ls = Ls,
        Ht = Ht, Lt = Lt, Delta = NULL, h = bw)
      
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lambda" = lambda, dt_risk_gammahat)
      rm(dt_risk_gammahat) ; gc()
      
      file_save <- paste0(script_path, "/FAR/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda,
                          "/dt_autocov_risk_FAR_mfBm_N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      saveRDS(object = dt_res, file = file_save)
      rm(dt_res) ; gc()
      return(paste0("Autocov risk : done for MC = ", mc_i, ", N = ", N, " and lambda = ", lambda))
    }, mc.cores = n_cores, dt_locreg = dt_locreg, N = N, lambda = lambda, white_noise = white_noise, design = design)
  } else if (white_noise == "fBm") {
    # Not need for the paper.
  }
  rm(locreg_file_name, dt_locreg, res) ; gc()
  
  return(paste0("Done : dt_autocov_risk_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

## Estimate mean function
estim_autocov_fun <- function(N = 150, lambda = 40, design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30){
  # Get the MC index in the autocov risk file name 
  data_risk_file_list <- list.files(paste0(script_path, "/FAR/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda))
  id_mc_data_risk <- gsub(pattern = paste0("dt_autocov_risk_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                          replacement = "", x = data_risk_file_list)
  id_mc_data_risk <- as.numeric(id_mc_data_risk)
  
  # Get the MC index in the data file name 
  data_file_list <- list.files(paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda))
  id_mc_data <- gsub(pattern = paste0("dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                     replacement = "", x = data_file_list)
  id_mc_data <- as.numeric(id_mc_data)
  
  # Get the already done MC repetition
  index_mc <- id_mc_data[which(id_mc_data_risk %in% id_mc_data)]
  
  # Load \widetilde{\gamma}_N(t) obtained from the script: 04_estimate_true_autocov.R 
  autocovtilde_file_name <- paste0(script_path, "/FAR/autocov_estimates/dt_autocovtilde_FAR_mfBm_N=5000_", design,".RDS")
  dt_autocovtilde <- readRDS(autocovtilde_file_name)
  
  # Estimate autocovariance by mc
  dt_autocov_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, Ni, lambdai, lag, s0, t0){
    ## Load the risk data
    id_mc_data_risk <- paste0(script_path, "/FAR/autocov_estimates/", design, "_autocov_risk/N", Ni, "lambda", lambdai, 
                              "/dt_autocov_risk_FAR_mfBm_N=", Ni, "_lambda=", lambdai, "_id_mc=", mc_i, "_", design, ".RDS")
    dt_autocov_risk <- readRDS(id_mc_data_risk)
    dt_optbw <- unique(dt_autocov_risk[!is.nan(autocov_risk), .("optbw" = h[which.min(autocov_risk)]), by = c("id_mc", "s", "t")])
    dt_optbw <- data.table::merge.data.table(x = dt_optbw, y = dt_autocovtilde, by = c("s", "t"))
    
    # Add optimal mean function
    ## Add \mu(t)
    dt_optbw[, c("muhat_opt_s", "mutrue_opt_s", "muhat_opt_t", "mutrue_opt_t") := .(0, 0, 0, 0)]
    
    ## Load raw data
    data_file_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda,
                             "/dt_mc_FAR_mfBm_N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
    dt <- readRDS(data_file_name)
    
    # Extract data
    dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
    optbw <- dt_optbw[id_mc == mc_i][order(s, t), optbw]
    muhat_opt_s <- dt_optbw[id_mc == mc_i][order(s, t), muhat_opt_s]
    muhat_opt_t <- dt_optbw[id_mc == mc_i][order(s, t), muhat_opt_t]
    
    # Estimate the mean function
    dt_autocov <- estimate_autocov(
      data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
      s = s0, t = t0, lag = lag, optbw = optbw, bw_grid = NULL,
      Hs = NULL, Ls = NULL, Ht = NULL, Lt = NULL,
      Delta = NULL, h = NULL, center = TRUE,
      mean_estimates_s = muhat_opt_s, mean_estimates_t = muhat_opt_t,
      smooth_ker = epanechnikov)
    
    # Add estimate of the true gamma
    dt_autocov <- data.table::merge.data.table(
      x = dt_autocov,
      y = dt_optbw[id_mc == mc_i, .(s, t, mutilde_s, mutilde_t, gammatilde, autocovtilde)],
      by = c("s", "t")
    )
    # Return and clean
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "time_taken" = time_taken, dt_autocov)
    return(dt_res)
  }, mc.cores = n_cores, Ni = N, lambdai = lambda, lag = lag, s0 = s0, t0 = t0))
  
  ## Save
  file_name <- paste0(script_path, "/FAR/autocov_estimates/dt_autocov_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_autocov_mc, file = file_name)
  rm(dt_autocovtilde, dt_autocov_mc) ; gc()
  
  return(paste0("Done : dt_autocov_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate autocovariance: FTS Model 3 with zero-mean ----

## Estimate autocovariance risk function
estim_autocov_risk_fun(N = 150, lambda = 40, white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
estim_autocov_risk_fun(N = 1000, lambda = 40, white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1, n_cores = 30)

## Estimate autocovariance function
estim_autocov_fun(N = 150, lambda = 40, design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
estim_autocov_fun(N = 1000, lambda = 40, design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30)

