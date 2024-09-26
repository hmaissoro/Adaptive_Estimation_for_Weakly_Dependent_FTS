################################################################################
# File:             05_estimate_true_autocov.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     23.01.2024
# 
# Contains R-script for numerical study: FTS Model 2:
#   FTS Model 2 with zero-mean \widetilde\{\gamma}_1(s,t) estimation
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

# Generate the data : FTS Model 2 with zero-mean ----
## Estimation points
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

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

## Design and N = number of curves of the FTS ----
design <- "d4"
N <- 5000

## Generate the data
dt_res <- data.table::rbindlist(parallel::mclapply(1:30, function(mc_i, Ni, t0, process_ker, process_mean, hurst){
  dt_gen <- simulate_far(N = Ni,
                         lambda = 40, # Not used if design = 'common'
                         tdesign = "common",
                         Mdistribution = NULL,
                         tdistribution = NULL,
                         tcommon = t0,
                         hurst_fun = hurst,
                         L = 4,
                         far_kernel = process_ker,
                         far_mean = process_mean,
                         int_grid = 100L,
                         burnin = 500L,
                         remove_burnin = TRUE)
  ## Change the name of the mean function
  data.table::setnames(dt_gen, "far_mean", "process_mean")
  dt_gen[, c("id_mc", "N") := .(mc_i, Ni)]
  data.table::setcolorder(x = dt_gen, neworder = c("id_mc", "N", "id_curve", "tobs", "ttag", "process_mean", "X"))
  
  return(dt_gen)
}, Ni = N, t0 = t0, process_ker = ker_d2, process_mean = mean_d4, hurst = Hlogistic, mc.cores = 30))

file_title <- paste0(script_path, "/FAR/autocov_estimates/dt_mc_common_FAR_mfBm_N=", N, "_", design,".RDS")

saveRDS(object = dt_res, file = file_title)
rm(dt_res, file_title) ; gc()

# Estimate \widetilde\{\gamma}_1(s,t) ----
t0 <- c(0.2, 0.4, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]
rm(dt_st) ; gc()

## Load data
data_file_name <- paste0(script_path, "/FAR/autocov_estimates/dt_mc_common_FAR_mfBm_N=", N, "_", design,".RDS")
dt <- readRDS(data_file_name)
index_mc <- dt[, unique(id_mc)]
  
## Estimate mean function
dt[, mutilde := mean(X), by = c("id_mc", "tobs")]

## Reshape the data
dt_autocovtilde_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, N, s0, t0, lag){
  dt_mc <- dt[id_mc == mc_i]
  ## Estimate \widetilde\gamma(s,t)
  dt_gammatilde_data <- data.table::rbindlist(lapply(dt[, unique(id_curve)], function(curve_index, data, s0, t0){
    dt_true_X <- data.table::merge.data.table(
      x = data.table::data.table("s" = s0, "t" = t0),
      y = data[id_curve == curve_index, .("t" = tobs, "Xt" = X)],
      by = "t"
    )
    dt_true_X <- data.table::merge.data.table(
      x = dt_true_X,
      y = data[id_curve == curve_index, .("s" = tobs, "Xs" = X)],
      by = "s"
    )
    return(data.table::data.table("id_curve" = curve_index, dt_true_X))
  }, data = dt_mc, s0 = s0, t0 = t0))
  
  ### Take into account the cross_lag
  ### The argument s is associated to the curves n = 1,..., N - lag
  dt_gammatilde_data_s <- dt_gammatilde_data[, list(id_curve, s, t, Xs)]
  dt_gammatilde_data_s <- dt_gammatilde_data_s[id_curve %in% 1:(N - lag)]
  dt_id_lag_s <- data.table::data.table(
    "id_curve" = sort(unique(dt_gammatilde_data_s[, id_curve])),
    "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
  dt_gammatilde_data_s <- data.table::merge.data.table(
    x = dt_id_lag_s,
    y = dt_gammatilde_data_s,
    by = "id_curve")
  
  ### The argument t is associated to the curves n = 1 + lag,..., N
  dt_gammatilde_data_t <- dt_gammatilde_data[, list(id_curve, s, t, Xt)]
  dt_gammatilde_data_t <- dt_gammatilde_data_t[id_curve %in% (1 + lag):N]
  dt_id_lag_t <- data.table::data.table(
    "id_curve" = sort(unique(dt_gammatilde_data_t[, id_curve])),
    "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
  dt_gammatilde_data_t <- data.table::merge.data.table(
    x = dt_id_lag_t,
    y = dt_gammatilde_data_t,
    by = "id_curve")
  
  ### Merge and clean
  dt_gammatilde_merge <- data.table::merge.data.table(
    x = dt_gammatilde_data_s,
    y = dt_gammatilde_data_t,
    by = c("id_lag", "s", "t"))
  dt_gammatilde_merge <- dt_gammatilde_merge[order(id_curve.x)]
  rm(dt_id_lag_s, dt_id_lag_t, dt_gammatilde_data_s, dt_gammatilde_data_t, dt_gammatilde_data) ; gc()
  
  ### The gamma function with true X's
  dt_gammatilde <- dt_gammatilde_merge[!(is.nan(Xs) | is.nan(Xt)), .("gammatilde" = mean(Xs * Xt)), by = c("s", "t")]
  rm(dt_gammatilde_merge) ; gc()
  
  ## Merge add mean function dt_gammatilde
  ### Add \mu(t)
  dt_autocovtilde <- data.table::merge.data.table(
    x = dt_gammatilde,
    y = unique(dt_mc[, .("t" = tobs, "mutilde_t" = mutilde)]),
    by = "t"
  )
  ### Add \mu(s)
  dt_autocovtilde <- data.table::merge.data.table(
    x = dt_autocovtilde,
    y = unique(dt_mc[, .("s" = tobs, "mutilde_s" = mutilde)]),
    by = "s"
  )
  dt_autocovtilde[, autocovtilde := gammatilde - mutilde_s * mutilde_t]
  
  dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lag" = lag, dt_autocovtilde)
  rm(dt_gammatilde, dt_autocovtilde) ; gc()
  return(dt_res)
}, mc.cores = 30, dt = dt, N = N, s0 = s0, t0 = t0, lag = 1))

## Estimate \widetilde\{\gamma}_1(s,t)
dt_autocovtilde_mc <- dt_autocovtilde_mc[
  ,  .("gammatilde" = mean(gammatilde),
       "mutilde_t" = mean(mutilde_t),
       "mutilde_s" = mean(mutilde_s),
       "autocovtilde" = mean(autocovtilde)),
  by = c("s", "t", "lag")]

## Save
file_name <- paste0(script_path, "/FAR/autocov_estimates/dt_autocovtilde_FAR_mfBm_N=", N, "_", design,".RDS")
saveRDS(object = dt_autocovtilde_mc, file = file_name)
rm(data_file_name, file_name, dt_autocovtilde_mc) ; gc()



