################################################################################
# File:             07_generate_paper_graphs.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     23.01.2024
# 
# Contains R-script for figures generation for FTS Model 2.
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

# Local regularity graph function ----
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}
t0 <- c(0.2, 0.4, 0.7, 0.8)

ggplot_locreg <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht", Hfun = Hlogistic){
  ## Load data
  file_name <- paste0(script_path, "/FAR/locreg_estimates/dt_locreg_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(file_name)
  
  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          # axis.title.y = element_text(size = 11, margin = margin(t = 10, r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'),
          legend.position = "top")
  
  if (white_noise == "mfBm") {
    ## define segment and set scale label
    dt_pr <- unique(dt_locreg[, .("t" = as.factor(t), "x" = as.factor(t - 0.05),
                                  "xend" = as.factor(t + 0.049), "Htrue" =  Hfun(t))])
    scale_label <- c(dt_pr[, as.character(t)], dt_pr[, as.character(x)], dt_pr[, as.character(xend)])
    scale_label <- sort(as.character(scale_label))
    scale_label[- which(scale_label %in% as.character(dt_pr[, t]))] <- as.character(" ")
    
    ## set t as factor
    dt_locreg[, t := as.factor(t)]
    if (param == "Ht") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      title_exp <- paste0("$\\widehat{H}_t$ - ", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.1, 0.9)
      x_lab <- "t"
      y_lab <- ""
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = "twodash", size = 0.9)
    } else if (param == "Lt") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)*
      title_exp <- paste0("$\\widehat{L}_t^2 - $", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 12)
      x_lab <- "t"
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = "twodash", size = 0.9)
      y_lab <- ""
      scale_label <- scale_label[! scale_label == " "]
    }
    
    ggplt <- ggplot(data = dt_locreg, mapping = aes(x = t, y = get(param))) +
      geom_boxplot() +
      ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_true_param +
      scale_x_discrete(labels = scale_label) +
      geom_theme
    
  } else if (white_noise == "fBm") {
    ## define segment and set scale label
    dt_pr <- unique(dt_locreg[, .("t" = as.factor(t), "x" = as.factor(Htrue - 0.05),
                                  "xend" = as.factor(Htrue + 0.05), Htrue)])
    scale_label <- c(dt_pr[, as.character(Htrue)], dt_pr[, as.character(x)], dt_pr[, as.character(xend)])
    scale_label <- sort(as.character(unique(scale_label)))
    scale_label[- which(scale_label %in% as.character(dt_pr[, as.factor(Htrue)]))] <- ""
    
    ## set t and Htrue as factor
    dt_locreg[, t := as.factor(t)]
    dt_locreg[, Htrue := as.factor(Htrue)]
    
    if (param == "Ht") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      title_exp <- paste0("$\\widehat{H}_t$ - ", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.1, 0.9)
      x_lab <- latex2exp::TeX("True $H_t$")
      y_lab <- ""
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = "twodash", linewidth = 0.9)
    } else if (param == "Lt") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      title_exp <- paste0("$\\widehat{L}_t^2$ - ", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 12)
      x_lab <- latex2exp::TeX("True $H_t$")
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = "twodash", linewidth = 0.9)
      y_lab <- ""
      scale_label <- scale_label[! scale_label == ""]
    }
    
    ggplt <- ggplot(data = dt_locreg, mapping = aes(x = Htrue, y = get(param), fill = t)) +
      geom_boxplot() +
      ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_true_param +
      scale_x_discrete(labels = scale_label) +
      scale_fill_grey() +
      geom_theme
  }
  return(ggplt)
}

## Local regularity estimates : FTS Model 2
g_locreg_far_mfBm_d2  <- ggpubr::ggarrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt"),
  nrow = 2, ncol = 4)
g_locreg_far_mfBm_d2
ggsave(
  filename = paste0(script_path, "/../paper_graphs/locreg_far_mfBm_d2.png"), plot = g_locreg_far_mfBm_d2,
  width = 9, height = 6, units = "in", bg = "white")

# Estimate mean function ----
## Merge the mean risk
merge_mean_risk <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2"){
  data_path <- paste0(script_path, "/FAR/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda)
  data_file_done <- list.files(data_path)
  dt_res <- data.table::rbindlist(lapply(data_file_done, function(ff){
    dt <- readRDS(file.path(data_path, ff))
  }))
  saveRDS(object = dt_res,
          file = paste0(script_path, "/FAR/mean_estimates/dt_mean_risk_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS"))
}

merge_mean_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d2")
merge_mean_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d2")
merge_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2")
merge_mean_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d2")

## Risk plot the merged risk function
ggplot_mean_risk_by_t <- function(N_vec = c(150, 1000, 400, 1000), lambda_vec = c(40, 40, 300, 1000),
                                  ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d2"){

  ## Load data, remove NaN values and reshape data
  dt_risk <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
    file_name <- paste0(script_path, "/FAR/mean_estimates/dt_mean_risk_",
                        process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
    dt_mean_risk <- readRDS(file_name)
    dt_mean_risk <- dt_mean_risk[! is.nan(mean_risk)]
    dt_mean_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk / 2)), by = c("N", "lambda", "t", "h")]
    dt_mean_risk[, N_lambda := paste0("(", N, ",", lambda, ")")]
    return(dt_mean_risk)
  }))

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))

  if (white_noise == "mfBm") {
    if (design == "d2") {
      y_lim <- c(0, 5)
    } else if (design == "d3"){
      y_lim <- c(0, 7)
    } else {
      y_lim <- c(0, 0.4 / 2)
    }
    title_exp <- paste0("t=", ti)
    ggplt <- ggplot(data = dt_risk[t == ti], mapping = aes(x = h, y = mean_risk, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
      geom_line(size = 0.9) +
      ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      labs(y = "", x = "h") +
      scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                            values = c("(1000,1000)" = "solid", "(400,300)" = "dotted", "(1000,40)" = "dotdash", "(150,40)" = "dashed"),
                            labels = c("(1000,1000)" = "(1000,1000)  ", "(400,300)" = "(400,300)  ", "(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                          values = c("(1000,1000)" = "#273746", "(400,300)" = "#34495E", "(1000,40)" = "#707B7C", "(150,40)" = "#909497"),
                          labels = c("(1000,1000)" = "(1000,1000)  ", "(400,300)" = "(400,300)  ", "(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      geom_theme
    return(ggplt)
  }
}

## Mean risk : FTS Model 2
g_mean_risk_far_mfBm_d2  <- ggpubr::ggarrange(
  ggplot_mean_risk_by_t(ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d2"),
  ggplot_mean_risk_by_t(ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d2"),
  ggplot_mean_risk_by_t(ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d2"),
  ggplot_mean_risk_by_t(ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d2"),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
g_mean_risk_far_mfBm_d2

ggsave(filename = paste0(script_path, "/../paper_graphs/mean_risk_far_mfBm_d2.png"),
       plot = g_mean_risk_far_mfBm_d2, width = 9, height = 6, units = "in", bg = "white")


# Mean function table ----
table_mean <- function(N_vec = c(150, 1000, 400, 1000), lambda_vec = c(40, 40, 300, 1000),
                       process = "FAR", white_noise = "mfBm", design = "d1"){

  ## Load data, remove NaN values and reshape data
  if (white_noise == "mfBm") {
    dt_mean <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
      file_name <- paste0(script_path, "/FAR/mean_estimates/dt_mean_estimates_",
                          process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
      dt_mean <- readRDS(file_name)
      dt_mean <- dt_mean[! (is.nan(muhat) | is.nan(muhat))]
      dt_mean <- dt_mean[, .("bias" = as.character(round(mean(muhat - mutrue), 4)),
                             "sd" = as.character(round(sd(muhat), 4))), by = c("N", "lambda", "t")]
      dt_mean[, N_lambda := paste0("(", N, ",", lambda, ")")]
      return(dt_mean)
    }))
    dt_mean_dcast <- data.table::dcast(data = dt_mean, formula = N + lambda ~ t, value.var = c("bias", "sd"))
    data.table::setcolorder(x = dt_mean_dcast,
                            neworder = c("N", "lambda", "bias_0.2", "sd_0.2", "bias_0.4", "sd_0.4",
                                         "bias_0.7", "sd_0.7", "bias_0.8", "sd_0.8"))
    dt_mean_dcast <- dt_mean_dcast[order(lambda)]

    Hmisc::latex(
      object = dt_mean_dcast,
      file =  paste0(script_path, "/../paper_tables/latex_mean_estimates_", process,"_", white_noise, "_", design,".tex"),
      n.cgroup = c(1, 1, 2, 2, 2, 2),
      cgroup = c("", "", "t = 0.2", "t = 0.4", "t = 0.7", "t = 0.8"),
      colheads = c("$N$", "$\\lambda$", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd"),
      rowname = NULL)
  } else if (white_noise == "fBm") {
    # Not need for the paper.
  }
}

table_mean(N_vec = c(150, 1000, 400, 1000), lambda_vec = c(40, 40, 300, 1000),
           process = "FAR", white_noise = "mfBm", design = "d2")

# Estimate autocovariance function : FTS Model 3 with zero-mean ----
merge_autocov_risk <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  data_path <- paste0(script_path, "/FAR/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda)
  data_file_done <- list.files(data_path)
  dt_res <- data.table::rbindlist(lapply(data_file_done, function(ff){
    dt <- readRDS(file.path(data_path, ff))
  }))
  saveRDS(object = dt_res,
          file = paste0(script_path, "/FAR/autocov_estimates/dt_autocov_risk_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS"))
}

merge_autocov_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4")
merge_autocov_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4")

## Risk function of the autocovariance function
ggplot_autocov_risk_by_st_d4 <- function(N_vec = c(150, 1000), lambda_vec = c(40, 40),
                                         si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d4"){

  ## Load data, remove NaN values and reshape data
  dt_risk <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
    ## Load data, remove NaN values and reshape data
    file_name <- paste0(script_path, "/FAR/autocov_estimates/dt_autocov_risk_",
                        process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
    dt_autocov_risk <- readRDS(file_name)
    dt_autocov_risk <- dt_autocov_risk[! is.nan(autocov_risk)]
    dt_autocov_risk <- dt_autocov_risk[, .("autocov_risk" = mean(autocov_risk / 3)), by = c("N", "lambda", "s", "t", "h")]
    dt_autocov_risk[, N_lambda := paste0("(", N, ",", lambda, ")")]
    return(dt_autocov_risk)
  }))

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))
  if (white_noise == "mfBm") {

    title_exp <- paste0("(s,t)=(", si, ",", ti, ")")
    ggplt <- ggplot(data = dt_risk[s == si & t == ti], mapping = aes(x = h, y = autocov_risk, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
      geom_line(size = 0.9) +
      ylim(0, 30) +
      ggtitle(latex2exp::TeX(title_exp)) +
      labs(y = "", x = "h") +
      scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                            values = c("(1000,40)" = "solid", "(150,40)" = "dashed"),
                            labels = c("(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                          values = c("(1000,40)" = "#34495E", "(150,40)" = "#909497"),
                          labels = c("(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      geom_theme
  } else if (white_noise == "mfBm") {
    # comming ...
  }

  return(ggplt)
}
g_autocov_risk_far_mfBm_d4  <- ggpubr::ggarrange(
  ggplot_autocov_risk_by_st_d4(si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d4"),
  ggplot_autocov_risk_by_st_d4(si = 0.8, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d4"),
  ggplot_autocov_risk_by_st_d4(si = 0.4, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d4"),
  ggplot_autocov_risk_by_st_d4(si = 0.7, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d4"),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
g_autocov_risk_far_mfBm_d4

ggsave(filename = paste0(script_path, "/../paper_graphs/autocov_risk_far_mfBm_d4.png"),
       plot = g_autocov_risk_far_mfBm_d4, width = 9, height = 6, units = "in", bg = "white")

# Autocov function table ----
table_autocov <- function(N_vec = c(150, 1000), lambda_vec = c(40, 40),
                          process = "FAR", white_noise = "mfBm", design = "d1",
                          s_vec = c(0.2, 0.8, 0.4, 0.7), t_vec = c(0.4, 0.2, 0.7, 0.8)){
  st_to_select <- paste0("(", s_vec, ",", t_vec, ")")
  ## Load data, remove NaN values and reshape data
  if (white_noise == "mfBm") {
    dt_autocov <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
      file_name <- paste0(script_path, "/FAR/autocov_estimates/dt_autocov_estimates_",
                          process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
      dt_autocov <- readRDS(file_name)
      dt_autocov <- dt_autocov[! (is.nan(gammahat) | is.nan(gammahat))]
      dt_autocov[, st := paste0("(", s, ",", t, ")")]
      dt_autocov <- dt_autocov[, .("bias" = as.character(round(mean(gammahat - gammatilde), 4)),
                                   "sd" = as.character(round(sd(gammahat), 4))), by = c("N", "lambda", "st")]
      dt_autocov <- dt_autocov[st %in% st_to_select]
      return(dt_autocov)
    }))
    dt_autocov_dcast <- data.table::dcast(data = dt_autocov, formula = N + lambda ~ st, value.var = c("bias", "sd"))
    data.table::setcolorder(x = dt_autocov_dcast,
                            neworder = c("N", "lambda", "bias_(0.2,0.4)", "sd_(0.2,0.4)", "bias_(0.4,0.7)", "sd_(0.4,0.7)",
                                         "bias_(0.7,0.8)", "sd_(0.7,0.8)", "bias_(0.8,0.2)", "sd_(0.8,0.2)"))
    dt_autocov_dcast <- dt_autocov_dcast[order(lambda)]

    Hmisc::latex(
      object = dt_autocov_dcast,
      file =  paste0(script_path, "/../paper_tables/latex_autocov_estimates_", process,"_", white_noise, "_", design,".tex"),
      n.cgroup = c(1, 1, 2, 2, 2, 2),
      cgroup = c("", "", "(s,t) = (0.2,0.4)", "(s,t) = (0.4,0.7)", "(s,t) = (0.7,0.8)", "(s,t) = (0.8,0.2)"),
      colheads = c("$N$", "$\\lambda$", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd"),
      rowname = NULL)
  } else if (white_noise == "fBm") {
    # Comming ...
  }
}

table_autocov(N_vec = c(150, 1000), lambda_vec = c(40, 40),
              process = "FAR", white_noise = "mfBm", design = "d4",
              s_vec = c(0.2, 0.8, 0.4, 0.7), t_vec = c(0.4, 0.2, 0.7, 0.8))





