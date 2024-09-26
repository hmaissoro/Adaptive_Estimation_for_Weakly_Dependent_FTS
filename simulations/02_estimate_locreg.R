################################################################################
# File:             02_estimate_locreg.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     23.01.2024
# 
# Contains R-script for local regularity estimation for numerical study: FTS Model 2.
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

# local regularity estimation function ----
t0 <- c(0.2, 0.4, 0.7, 0.8)
estim_locreg_fun <- function(N = 150, lambda = 40, white_noise = "mfBm", design = "d2", t0, nc_cores = 30){

  if (white_noise == "mfBm") {
    data_file_list <- list.files(paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)

    # Nmc <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda)))
    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(id_mc_data, function(mc_i, Ni, lambdai, t0){
      data_file_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", Ni, "lambda", lambdai,
                               "/dt_mc_FAR_mfBm_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]

      # Extract and sort data
      dt <- dt[order(id_curve, tobs)]
      bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

      # Estimate the local regularity
      ## Delta
      lambdahat <- mean(dt[, .N, by = id_curve][, N])
      delta <- min(exp(- log(lambdahat) ** ( 1 / 3)), 0.2)

      ## Centered process
      dt_locreg <- estimate_locreg(
        data = dt, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, Delta = delta, h = bw,
        smooth_ker = epanechnikov, center = TRUE)

      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg)
      rm(dt_locreg) ; gc()
      return(dt_res)
    }, mc.cores = nc_cores, Ni = N, lambdai = lambda, t0 = t0))

  } else if (white_noise == "fBm") {
    # Not Need for FTS Model 2
  }
  ## Save
  file_name <- paste0(script_path, "/FAR/locreg_estimates/dt_locreg_FAR_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_reg_mc, file = file_name)
  rm(file_name, dt_reg_mc) ; gc() ; gc()

  return(paste0("Done : dt_locreg_FAR_mfBm_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate local regularity : FTS Model 2 ----
estim_locreg_fun(N = 150, lambda = 40, white_noise = "mfBm", design = "d2", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 40, white_noise = "mfBm", design = "d2", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 1000, white_noise = "mfBm", design = "d2", t0 = t0, nc_cores = 30)


# Estimate local regularity : FTS Model 2 with zero ----
estim_locreg_fun(N = 150, lambda = 40, white_noise = "mfBm", design = "d4", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 40, white_noise = "mfBm", design = "d4", t0 = t0, nc_cores = 30)

