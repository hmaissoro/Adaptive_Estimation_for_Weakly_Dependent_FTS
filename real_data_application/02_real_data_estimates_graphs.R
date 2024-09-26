library(data.table)
library(ggplot2)

dt_real_data <- fread("./inst/12_mc_simulate_data/real_data/data_voltage.csv")
## Global empirical mean
dt_real_data <- data.table::fread("./inst/12_mc_simulate_data/real_data/data_voltage.csv")
dt_mean_emp <- dt_real_data[, .("muhat" = mean(X)), by = "tobs"]
saveRDS(dt_mean_emp, "./inst/12_mc_simulate_data/real_data/dt_emp_mean.RDS")

dt_mean <- readRDS("./inst/12_mc_simulate_data/real_data/dt_mean_real_data.RDS")
dt_mean[, c("N", "lambda") := .(1358, 1440)]
dt_mean_subset <- readRDS("./inst/12_mc_simulate_data/real_data/dt_mean_real_data_subset.RDS")
dt_mean_subset[, c("N", "lambda") := .(50, 1440)]

dt_mean_all <- rbind(dt_mean, dt_mean_subset)

# plot and save estimated curves : all data ----

# Smooth Local regularity estimates
tobs <- dt_mean[order(t)][, t]
Ht <- dt_mean[order(t)][, Ht]
Lt <- dt_mean[order(t)][, Lt]
K <- 20

## Define basis
mat_covariable <- splines::ns(x = tobs, df = (K + 1 + 1), intercept = TRUE)

### LASSO regression
cv_model_Ht <- glmnet::cv.glmnet(
  x = mat_covariable, y = Ht,
  intercept = TRUE, alpha = 1, nfolds = 20
)
plot(cv_model_Ht)

cv_model_Lt <- glmnet::cv.glmnet(
  x = mat_covariable, y = Lt,
  intercept = TRUE, alpha = 1, nfolds = 20
)
plot(cv_model_Lt)

best_lambda_Ht <- cv_model_Ht$lambda.min
best_lambda_Lt <- cv_model_Lt$lambda.min

Ht_model <- glmnet::glmnet(
  x = mat_covariable, y = Ht,
  intercept = FALSE, alpha = 1, lambda = best_lambda_Ht
)

Lt_model <- glmnet::glmnet(
  x = mat_covariable, y = Lt,
  intercept = FALSE, alpha = 1, lambda = best_lambda_Lt
)


Hthat <- predict(Ht_model, mat_covariable)
Lthat <- predict(Lt_model, mat_covariable)

dygraphs::dygraph(data = data.table::data.table(tobs, Hthat))
dygraphs::dygraph(data = data.table::data.table(tobs, Lthat))

dt_res_locreg <- data.table::data.table("N" = 1358, "lambda" = 1440, "t" = tobs, "Ht" = c(Hthat), "Lt" = c(Lthat))

# plot and save estimated curves : all data
## H_t
g_local_exponent <- ggplot(data = dt_res_locreg, aes(x = t, y = Ht)) +
  geom_line(size = 0.5)  +
  ylim(0, 0.6) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$H_t$  ")) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_local_exponent
ggsave(plot = g_local_exponent, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_Ht_real_data.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## Lt^2
g_local_constant <- ggplot(data = dt_res_locreg, aes(x = t, y = Lt)) +
  geom_line(size = 0.5)  +
  ylim(0, 60) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$L_t^2$  ")) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_local_constant
ggsave(plot = g_local_constant, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_Lt_real_data.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# Plot mean curves
## Real data mean
tobs_locreg <- tobs
tobs_all <- dt_real_data[, sort(unique(tobs))]

dt_real_data_silice <- dt_real_data[tobs %in% tobs_all[round(tobs_locreg * 1440)]]
dt_real_data_silice[, muhat_emp := mean(X), by = "tobs"]

dt_mean_emp <- unique(dt_real_data_silice[, .("t" = tobs, muhat_emp)])

g_mean_emp <- ggplot(data = dt_mean_emp, aes(x = t, y = muhat_emp)) +
  geom_line(size = 0.5)  +
  ylim(237, 244) +
  xlab(label = "t") +
  ylab(label = "") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_mean_emp
ggsave(plot = g_mean_emp, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_empirique_real_data.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## Our mean function estimates

g_muhat <- ggplot(data = dt_mean, aes(x = t, y = muhat)) +
  geom_line(size = 0.5)  +
  ylim(237, 244) +
  xlab(label = "t") +
  ylab(label ="") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_muhat
ggsave(plot = g_muhat, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_MPV_real_data.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

### N = 50 last curves
g_muhat_subset <- ggplot(data = dt_mean_subset, aes(x = t, y = muhat)) +
  geom_line(size = 0.5)  +
  ylim(236, 246) +
  xlab(label = "t") +
  ylab(label ="") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_muhat_subset
ggsave(plot = g_muhat_subset, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_MPV_real_data_subset_N=50.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## Rubin and Panateros mean function estimator
dt_mean_rp <- fread("../../matlab/rubin_mean_autocov_method/real_data_muhat_RP.csv")
dt_mean_rp[, c("N", "lambda") := .(1358, 1440)]

g_muhat_RP <- ggplot(data = dt_mean_rp, aes(x = t, y = muhat_RP)) +
  geom_line(size = 0.5)  +
  ylim(237, 244) +
  xlab(label = "t") +
  ylab(label ="") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_muhat_RP

ggsave(plot = g_muhat_RP, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_RP_real_data.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# plot and save estimated curves : subset of N=50 curves ----
dt_mean_subset_subset <- readRDS("./inst/12_mc_simulate_data/real_data/dt_mean_real_data_subset.RDS")

# Smooth Local regularity estimates
tobs <- dt_mean_subset[order(t)][, t]
Ht <- dt_mean_subset[order(t)][, Ht]
Lt <- dt_mean_subset[order(t)][, Lt]
K <- 20 / 2

## Define basis
mat_covariable <- splines::ns(x = tobs, df = (K + 1 + 1), intercept = TRUE)

### LASSO regression
cv_model_Ht <- glmnet::cv.glmnet(
  x = mat_covariable, y = Ht,
  intercept = TRUE, alpha = 1, nfolds = 20
)
plot(cv_model_Ht)

cv_model_Lt <- glmnet::cv.glmnet(
  x = mat_covariable, y = Lt,
  intercept = TRUE, alpha = 1, nfolds = 20
)
plot(cv_model_Lt)

best_lambda_Ht <- cv_model_Ht$lambda.min
best_lambda_Lt <- cv_model_Lt$lambda.min

Ht_model <- glmnet::glmnet(
  x = mat_covariable, y = Ht,
  intercept = FALSE, alpha = 1, lambda = best_lambda_Ht
)

Lt_model <- glmnet::glmnet(
  x = mat_covariable, y = Lt,
  intercept = FALSE, alpha = 1, lambda = best_lambda_Lt
)


Hthat <- predict(Ht_model, mat_covariable)
Lthat <- predict(Lt_model, mat_covariable)

dygraphs::dygraph(data = data.table::data.table(tobs, Hthat))
dygraphs::dygraph(data = data.table::data.table(tobs, Lthat))

dt_res_locreg_subset <- data.table::data.table("N" = 50, "lambda" = 1440, "t" = tobs, "Ht" = c(Hthat), "Lt" = c(Lthat))

# plot and save estimated curves
## H_t
g_local_exponent <- ggplot(data = dt_res_locreg_subset, aes(x = t, y = Ht)) +
  geom_line(size = 0.5)  +
  ylim(0, 0.45) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$H_t$  ")) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_local_exponent
ggsave(plot = g_local_exponent, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_Ht_real_data_subset_N=50.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## Lt^2
g_local_constant <- ggplot(data = dt_res_locreg_subset, aes(x = t, y = Lt)) +
  geom_line(size = 0.5)  +
  ylim(0, 50) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$L_t^2$  ")) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_local_constant
ggsave(plot = g_local_constant, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_Lt_real_data_subset_N=50.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## Our mean function estimates

### N = 50 last curves
g_muhat_subset <- ggplot(data = dt_mean_subset, aes(x = t, y = muhat)) +
  geom_line(size = 0.5)  +
  ylim(236, 246) +
  xlab(label = "t") +
  ylab(label ="") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_muhat_subset
ggsave(plot = g_muhat_subset, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_MPV_real_data_subset_N=50.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## Rubin and Panateros mean function estimator
dt_mean_rp_subset <- fread("../../matlab/rubin_mean_autocov_method/real_data_muhat_RP_subset.csv")
dt_mean_rp_subset[, c("N", "lambda") := .(50, 1440)]

g_muhat_RP <- ggplot(data = dt_mean_rp_subset, aes(x = t, y = muhat_RP)) +
  geom_line(size = 0.5)  +
  ylim(236, 246) +
  xlab(label = "t") +
  ylab(label ="") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_muhat_RP

ggsave(plot = g_muhat_RP, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_RP_real_data_subset_N=50.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# plot and save estimated curves : all curves AND subset of N=50 curves ----
dt_mean_RP_all <- rbind(dt_mean_rp, dt_mean_rp_subset)
dt_mean_RP_all[, "N_lambda" := paste0("(", N, ",", lambda,")")]
dt_res_locreg_all <- rbind(dt_res_locreg, dt_res_locreg_subset)
dt_res_locreg_all[, "N_lambda" := paste0("(", N, ",", lambda,")")]
dt_mean_all <- rbind(dt_mean, dt_mean_subset)
dt_mean_all[, "N_lambda" := paste0("(", N, ",", lambda,")")]


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

# Local regularity parameters
## H_t
ggplt_Ht <- ggplot(data = dt_res_locreg_all, mapping = aes(x = t, y = Ht, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
  geom_line(size = 0.9) +
  ylim(0, 0.6) +
  # ggtitle(latex2exp::TeX("$t \\to H_t$")) +
  labs(y = "", x = "h") +
  scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                        values = c("(1358,1440)" = "solid", "(50,1440)" = "dotdash"),
                        labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                      values = c("(1358,1440)" = "black", "(50,1440)" = "#A93226"),
                      labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  geom_theme
ggplt_Ht
ggsave(plot = ggplt_Ht, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_Ht_real_data_all.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## L_t^2
ggplt_Lt <- ggplot(data = dt_res_locreg_all, mapping = aes(x = t, y = Lt, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
  geom_line(size = 0.9) +
  ylim(0, 60) +
  labs(y = "", x = "h") +
  scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                        values = c("(1358,1440)" = "solid", "(50,1440)" = "dotdash"),
                        labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                      values = c("(1358,1440)" = "black", "(50,1440)" = "#A93226"),
                      labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  geom_theme
ggplt_Lt
ggsave(plot = ggplt_Lt, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_Lt_real_data_all.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# Mean function
## Our mean function
ggplt_mean_our <- ggplot(data = dt_mean_all, mapping = aes(x = t, y = muhat, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
  geom_line(size = 0.9) +
  ylim(236, 246) +
  labs(y = "", x = "h") +
  scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                        values = c("(1358,1440)" = "solid", "(50,1440)" = "dotdash"),
                        labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                      values = c("(1358,1440)" = "black", "(50,1440)" = "#A93226"),
                      labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  geom_theme
ggplt_mean_our
ggsave(plot = ggplt_mean_our, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_MPV_real_data_all.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

## RP mean function
ggplt_mean_rp <- ggplot(data = dt_mean_RP_all, mapping = aes(x = t, y = muhat_RP, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
  geom_line(size = 0.9) +
  ylim(236, 246) +
  labs(y = "", x = "h") +
  scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                        values = c("(1358,1440)" = "solid", "(50,1440)" = "dotdash"),
                        labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                      values = c("(1358,1440)" = "black", "(50,1440)" = "#A93226"),
                      labels = c("(1358,1440)" = "(1358,1440)  ", "(50,1440)" = "(50,1440)  ")) +
  geom_theme
ggplt_mean_rp
ggsave(plot = ggplt_mean_rp, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_estimate_RP_real_data_all.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")


