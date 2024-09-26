################################################################################
# File:             04_mean_estimates_normal_qq_plots.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     23.01.2024
# 
# Contains R-script for mean estimation for numerical study: FTS Model 2.
#         Normal Q-Q plots
# 
# For a detailed description see:
#   Maissoro, H., Patilea, V. and Vimond, M. (2024). Adaptive estimation for 
#     weakly dependent functional time series.
################################################################################

# Load packages and functions
library(data.table)
library(matrixStats)
library(ggplot2)

all_func <- list.files("../R/")
lapply(all_func, function(func){
  if (func != "zzz.R") source(file = file.path("../R/", func), echo = FALSE)
})

# Get the script path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# Estimate sandardized mean distribution  ----
t0 <- c(0.2, 0.4, 0.7, 0.8)
estim_standardised_mean_qqplot_fun <- function(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", t0, n_cores = 30){
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
    data_risk_file <- paste0(script_path, "/FAR/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda, 
                             "/dt_mean_risk_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", mc_i, "_", design, ".RDS")
    dt_mean_risk <- readRDS(data_risk_file)
    dt_optbw <- dt_mean_risk[, .("optbw" = h[which.min(mean_risk)] ** 1.1), by = c("id_mc", "t")]
    
    # Load data
    data_file_name <- paste0(script_path, "/FAR/data/", design, "_slice/N", N, "lambda", lambda,
                             "/dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
    dt <- readRDS(data_file_name)
  
    # Extract data
    dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
    id_curves <- dt_random_mc[, unique(id_curve)]
    optbw <- dt_optbw[id_mc == mc_i][order(t), optbw]
    t0 <- dt_optbw[id_mc == mc_i][order(t), t]

    # Estimate sigma
    dt_sig <- estimate_sigma(data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X", t = t0)
    dt_sig_mutrue <- data.table::merge.data.table(
      x = dt_sig,
      y = unique(dt[ttag == "tcommon"][id_mc == mc_i, .("t" = tobs, "mutrue" = process_mean)]),
      by = "t")

    # Smooth curves with optimal bandwidth parameters
    dt_Xhat <- data.table::rbindlist(lapply(id_curves, function(curve_index, t, data, optbw, smooth_ker){
      # \pi_n(t,h)
      Tn <- data[id_curve == curve_index, tobs]
      Yn <- data[id_curve == curve_index, X]
      pi_n <- sapply(X = 1:length(t), function(tidx, Tn, t, optbw_mean){
        as.numeric(abs(Tn - t[tidx]) <= optbw_mean[tidx])
      }, Tn = Tn, t = t, optbw_mean = optbw)
      pi_n <- t(pi_n)
      pi_n <- as.numeric(rowSums(pi_n, na.rm = TRUE) >= 1)

      # \sum{i=1}^{M_n} W_{n,i}(t,h)^2
      Wn_square <- sapply(X = 1:length(t), function(tidx, Tn, t, optbw){
        K_ni <- smooth_ker((Tn - t[tidx]) / optbw[tidx])
        kn_sum <- sum(K_ni)
        if (kn_sum == 0) {
          W_ni <- K_ni
        } else {
          W_ni <- {K_ni / kn_sum} ** 2
        }
        return(W_ni)
      }, Tn = Tn, t = t, optbw = optbw)
      Wn_square <- t(Wn_square)
      Wn_square <- rowSums(Wn_square)

      # \widehat X(t;h)
      Xhat <- mapply(function(t, h, Yn, Tn, ker){
        estimate_nw(y = Yn, t = Tn, tnew = t, h = h, smooth_ker = ker)$yhat
      }, t = t, h = optbw,
      MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

      dt_res <- data.table::data.table("id_curve" = curve_index, "t" = t, "Xhat" = Xhat, "pi_n" = pi_n, "Wn_square" = Wn_square, "sigma" = dt_sig$sig)
      return(dt_res)
    }, data = dt_random_mc, t = t0, optbw = optbw, smooth_ker = epanechnikov))

    # Estimate mean function \widehat \mu(t)
    dt_Xhat[is.nan(Xhat) & pi_n == 0, Xhat := 0]

    # Add sigma estimates and true mean function
    dt_Xhat <- data.table::merge.data.table(x = dt_Xhat, y = dt_sig_mutrue, by = "t")
    dt_Xhat[, PN := sum(pi_n), by = t]
    dt_Xhat <- dt_Xhat[!is.nan(Xhat)]
    dt_Xhat[, "muhat" := sum(pi_n * Xhat) / PN, by = "t"]
    dt_Xhat[, "Sigma_element" := sum(pi_n * (Wn_square ** 2) * (sigma ** 2)) / PN, by = "t"]
    dt_Xhat[, "SSmu_element" := sum(pi_n * (Xhat - mutrue) / sqrt(PN)), by = "t"]

    # Add the optimal bandwidth
    dt_muhat <- unique(dt_Xhat[, list(t, muhat, mutrue, PN, Sigma_element, SSmu_element)])
    dt_muhat[, "id_mc" := mc_i]
    dt_muhat <- data.table::merge.data.table(x = dt_optbw[, .(id_mc, t, "hN" = optbw)], y = dt_muhat, by = c("id_mc", "t"))

    return(dt_muhat)
  }, mc.cores = n_cores, Ni = N, lambdai = lambda))

  # Estimate the variance part : \mathbb S_mu(t)
  dt_mean_mc[, SSmu_t := var(SSmu_element), by = "t"]

  # Estimate the variance part : \Sigma(t)
  dt_mean_mc[, Sigma_t := mean(Sigma_element), by = "t"]

  # Estimate the standarize mean estimation
  dt_mean_mc[, mean_standardise := sqrt(PN) * (muhat - mutrue) / sqrt(Sigma_t + SSmu_t),  by = "t"]

  ## Save
  file_name <- paste0(script_path, "/FAR/mean_estimates/dt_mean_standardised_qqplot_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(file_name, dt_mean_mc) ; gc()

  return(paste0("Done : dt_mean_standardised_qqplot_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Q-Q plots :FTS Model 2
estim_standardised_mean_qqplot_fun(N = 150, lambda = 40, white_noise = "mfBm", design = "d2", n_cores = 30)
estim_standardised_mean_qqplot_fun(N = 1000, lambda = 40, white_noise = "mfBm", design = "d2", n_cores = 30)
estim_standardised_mean_qqplot_fun(N = 400, lambda = 300, white_noise = "mfBm", design = "d2", n_cores = 30)
estim_standardised_mean_qqplot_fun(N = 1000, lambda = 1000, white_noise = "mfBm", design = "d2", n_cores = 30)

# Q-Q plots
ggplot_mean_standardised_qq_by_t <- function(N = 400, lambda = 300, ti = 0.2, white_noise = "mfBm", design = "d2"){
  ## Load data
  file_name <- paste0(script_path, "/FAR/mean_estimates/dt_mean_standardised_qqplot_estimates_FAR_mfBm_N=",
                      N, "_lambda=", lambda, "_", design, ".RDS")
  dt_mean <- readRDS(file_name)
  dt_mean <- unique(dt_mean)

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))

  title_exp <- paste0("t=", ti, ", ", "N=", N, ", $\\lambda$=", lambda)
  dt_qq <- dt_mean[t == ti]
  dt_qq <- dt_qq[!is.nan(mean_standardise)]
  ggplt <-  ggplot(data = dt_qq, aes(sample = mean_standardise, group = 1)) +
    # ylim(-18, 18) +
    stat_qq(data = dt_qq, aes(sample = mean_standardise, group = 1)) +
    stat_qq_line(data = dt_qq, aes(sample = mean_standardise, group = 1)) +
    labs(y = "", x = "N(0,1)") +
    ggtitle(latex2exp::TeX(title_exp)) +
    geom_theme
  return(ggplt)
}

# Scenario 1 :
for(ti in c(0.2, 0.4, 0.7, 0.8)){
  g_mean_qq_sd_far_mfBm_d2  <- ggpubr::ggarrange(
    ggplot_mean_standardised_qq_by_t(N = 150, lambda = 40, ti = ti, white_noise = "mfBm", design = "d2"),
    ggplot_mean_standardised_qq_by_t(N = 1000, lambda = 40, ti = ti, white_noise = "mfBm", design = "d2"),
    ggplot_mean_standardised_qq_by_t(N = 400, lambda = 300, ti = ti, white_noise = "mfBm", design = "d2"),
    ggplot_mean_standardised_qq_by_t(N = 1000, lambda = 1000, ti = ti, white_noise = "mfBm", design = "d2"),
    nrow = 2, ncol = 2)

  ggsave(filename = paste0(script_path, "/../paper_graphs/mean_standardised_qq_far_mfBm_t=", ti,"_d2.png"),
         plot = g_mean_qq_sd_far_mfBm_d2, width = 9, height = 7, units = "in", bg = "white")
  rm(g_mean_qq_sd_far_mfBm_d2)
}
