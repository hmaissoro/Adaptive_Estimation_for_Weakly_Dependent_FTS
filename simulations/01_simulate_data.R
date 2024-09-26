################################################################################
# File:             01_simulate_data.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     22.01.2024
# 
# Contains R-script for data generation for numerical study: FTS Model 2.
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

# Simulation global parameters----
sig <- 0.25
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

# Simulation function ----
simulate_data_fun <- function(mc_i, Ni, lambdai, t0, sig = 0.25,
                              process_ker = get_real_data_far_kenel,
                              process_mean = get_real_data_mean,
                              white_noise = "mfBm", Lt = 4,
                              hurst = Hlogistic, Hconst = 0.4){
  # Generate data according the choixe of the white noise
  if (white_noise == "mfBm") {
    ## Generate FAR(1)
    dt_gen <- simulate_far(N = Ni, lambda = lambdai,
                           tdesign = "random",
                           Mdistribution = bounded_uniform,
                           tdistribution = runif,
                           tcommon = t0,
                           hurst_fun = hurst,
                           L = Lt,
                           far_kernel = process_ker,
                           far_mean = process_mean,
                           int_grid = 100L,
                           burnin = 100L,
                           remove_burnin = TRUE)
  } else if (white_noise == "fBm") {
    hurst <-  function(t) Hconst + 0 * t
    ## Generate FAR(1)
    dt_gen <- simulate_far(N = Ni, lambda = lambdai,
                           tdesign = "random",
                           Mdistribution = bounded_uniform,
                           tdistribution = runif,
                           tcommon = t0,
                           hurst_fun = hurst,
                           L = Lt,
                           far_kernel = process_ker,
                           far_mean = process_mean,
                           int_grid = 100L,
                           burnin = 100L,
                           remove_burnin = TRUE)
  }
  ## Change the name of the mean function
  data.table::setnames(dt_gen, "far_mean", "process_mean")
  
  ## Add noise
  dt_gen[ttag == "trandom", X := X + rnorm(n = .N, mean = 0, sd = sig), by = "id_curve"]
  
  ## Get pre-smoothing bandwidth
  ### Define and exponential bandwidth grid
  lambdahat <- mean(dt_gen[ttag == "trandom", .N, by = id_curve][, N])
  K <- 30
  b0 <- 1 / lambdahat
  bK <- lambdahat ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(b0, bK, a, K) ; gc()
  
  ### Get optimal bw for each curve
  dt_gen[ttag == "trandom", X := X + rnorm(n = .N, mean = 0, sd = sig), by = "id_curve"]
  index <- dt_gen[, unique(id_curve)]
  index_last30 <- tail(sort(index), 30)
  
  presmooth_bw <- median(unlist(lapply(index_last30, function(id, dtt, bw_grid){
    # Filter data
    d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]
    
    # Get optimal bandwidth
    bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
    return(bw)
  }, dtt = dt_gen, bw_grid = bw_grid)))
  
  dt_gen[, presmooth_bw := presmooth_bw]
  
  # Add MC index
  if (white_noise == "mfBm") {
    dt_gen[, c("id_mc", "N", "lambda") := .(mc_i, Ni, lambdai)]
    data.table::setcolorder(
      x = dt_gen,
      neworder = c("id_mc", "N", "lambda", "id_curve", "tobs", "ttag", "process_mean", "X", "presmooth_bw"))
    
    # Save RDS
    file_save_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", Ni, "lambda", lambdai,
                             "/dt_mc_FAR_mfBm_N=", Ni, "_lambda=", lambdai, "_id_mc=", mc_i, "_", design,".RDS")
    
  } else if (white_noise == "fBm") {
    dt_gen[, c("id_mc", "N", "lambda", "Htrue") := .(mc_i, Ni, lambdai, Hconst)]
    data.table::setcolorder(
      x = dt_gen,
      neworder = c("id_mc", "N", "lambda", "Htrue", "id_curve", "tobs", "ttag", "process_mean", "X", "presmooth_bw")
    )
    # Save RDS
    file_save_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", Ni, "lambda", lambdai, "/H", Hconst,
                             "/dt_mc_FAR_fBm_H=", Hconst, "_N=", Ni, "_lambda=", lambdai, "_id_mc=", mc_i, "_", design,".RDS")
    
  }
  saveRDS(object = dt_gen, file = file_save_name)
  
  # return the result
  return(paste0("Done :", file_save_name))
}

# Simulate all MC process

simulate_data <- function(Nmc = 1:400, Ni = 400, lambdai = 300, t0, sig = 0.25,
                          process_ker = get_real_data_far_kenel,
                          process_mean = get_real_data_mean,
                          white_noise = "mfBm", Lt = 4,
                          hurst = Hlogistic, Hconst = Hconst, design = "d2"){
  parallel::mclapply(Nmc, function(mc_i, Ni, lambdai, t0, sig,
                                   process_ker, process_mean, white_noise, hurst, Hconst){
    dt_ <- simulate_data_fun(mc_i = mc_i, Ni = Ni, lambdai = lambdai, t0 = t0, sig = sig,
                             process_ker = process_ker, process_mean = process_mean, 
                             white_noise = white_noise, Lt = Lt,
                             hurst = hurst, Hconst = Hconst)
    return(dt_)
  }, Ni = Ni, lambdai = lambdai, t0 = t0, sig = sig,
  process_ker = process_ker, process_mean = process_mean,
  white_noise = white_noise, hurst = hurst, Hconst = Hconst, mc.cores = 30)
  
  print(paste0("Done : dt_mc_FAR_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", design,".RDS at ", Sys.time()))
}

# Data generation : FTS Model 2 ----

## Mean function
mean_d2 <- function(t) 4 * sin(3 * pi * t / 2)

## Autoregressive kernel
ker_d2 <- function(s,t, operator_norm = 0.5){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}


simulate_data(Nmc = 1:400, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process_ker = ker_d2, process_mean = mean_d2,
              white_noise = "mfBm", hurst = Hlogistic,
              Hconst = Hconst, design = "d2")

simulate_data(Nmc = 1:400, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process_ker = ker_d2, process_mean = mean_d2,
              white_noise = "mfBm", hurst = Hlogistic,
              Hconst = Hconst, design = "d2")

simulate_data(Nmc = 1:400, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process_ker = ker_d2,  process_mean = mean_d2,
              white_noise = "mfBm", hurst = Hlogistic, 
              Hconst = Hconst, design = "d2")

simulate_data(Nmc = 1:30, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process_ker = ker_d2, process_mean = mean_d2,
              white_noise = "mfBm", hurst = Hlogistic,
              Hconst = Hconst, design = "d2")

# Data generation : FTS Model 2 with zero-mean ----

## Mean function
mean_d4 <- function(t) 0

## Autoregressive kernel
ker_d2 <- function(s,t, operator_norm = 0.5){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}


simulate_data(Nmc = 1:400, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process_ker = ker_d2,
              process_mean = mean_d4, white_noise = "mfBm",
              hurst = Hlogistic, Hconst = Hconst, design = "d4")

simulate_data(Nmc = 1:400, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process_ker = ker_d2,
              process_mean = mean_d4, white_noise = "mfBm",
              hurst = Hlogistic, Hconst = Hconst, design = "d4")




