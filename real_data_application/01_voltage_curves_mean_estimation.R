################################################################################
# File:             01_voltage_curves_mean_estimation.R
# Created by:       Hassan Maissoro
# First released:   15.01.2024
# Last revised:     23.01.2024
# 
# Contains R-script for mean estimation for numerical study:
#     Voltage curves mean function estimation
# 
# For a detailed description see:
#   Maissoro, H., Patilea, V. and Vimond, M. (2024). Adaptive estimation for 
#     weakly dependent functional time series.
################################################################################

# Load packages and functions
library(data.table)
library(matrixStats)

all_func <- list.files("./package/adaptivefts/R/")
lapply(all_func, function(func){
  if (func != "zzz.R") source(file = file.path("./package/adaptivefts/R/", func), echo = FALSE)
})

# Get the script path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# Import data ----
dt_raw <- fread(file.path(script_path, "/household_power_consumption.txt"), sep = ";")

# Information about the data
names(dt_raw)
dt_raw[, length(unique(Date))]

# One observation per minute, expect "16/12/2006" and "26/11/2010" where some minutes are not observed
N_row_curves <- length(unique(dt_raw$Date))
dt_raw[, .N, by = "Date"][, summary(N)]
dt_raw[! Date %in% c("16/12/2006", "26/11/2010"), .N, by = "Date"]

# Remove "16/12/2006" and "26/11/2010"
dt_raw <- dt_raw[! Date %in% c("16/12/2006", "26/11/2010")]
dt_raw[, date := as.Date(Date, format = "%d/%m/%Y")]

# Normalize time
dt_raw[, t := (1:1440) / 1440, by = "Date"]
dt <- dt_raw[, .(date, "tobs" = t, "voltage" = Voltage)]
rm(dt_raw) ; gc()

# Convert voltage curve as nunmeric
dt[, voltage := as.numeric(voltage)]
date_na <- dt[is.na(voltage), unique(date)]
dt <- dt[! date %in% date_na]
dt[, length(unique(date))]

# Add curve index
date_vec <- unique(dt$date)
dt_id_curve <- data.table::data.table("id_curve" = 1:1358, "date" = date_vec)
dt_real_data <- data.table::merge.data.table(x = dt, dt_id_curve, by = "date")
rm(dt, dt_id_curve) ; gc()

# Estimate loacal regularity parameters ----
## Estimate bandwidth
N <- dt_real_data[, max(id_curve)]
(1 - N / N_row_curves) * 100 # %tage of missing days
lambda <- 1440
K <- 30
b0 <- 1 / lambda
bK <- lambda ** (- 1 / 3)
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(b0, bK, a, K) ; gc()

index_last30 <- dt_real_data[, tail(sort(unique(id_curve)), 30)]
presmooth_bw <- median(unlist(parallel::mclapply(index_last30, function(curve_index, data, bw_grid){
  dt_by_curve <- data[id_curve == curve_index]
  bw <- estimate_nw_bw(
    y = dt_by_curve[, voltage], t = dt_by_curve[, tobs],
    bw_grid = bw_grid, smooth_ker = epanechnikov
  )
  return(bw)
}, mc.cores = 30, data = dt_real_data, bw_grid)))

## Add presmoothing bandwidth 
dt_real_data[, presmooth_bw := presmooth_bw]

## Estimate local regularity parameters
lambda <- 1440
delta <- exp(-log(lambda) ** (1/3))
bw <- unique(dt_real_data[order(id_curve), .(id_curve, presmooth_bw)])[, presmooth_bw]
t0 <- seq(0.1, 0.99, len = 96)

dt_locreg <- estimate_locreg(
  data = dt_real_data[order(id_curve)], idcol = "id_curve", tcol = "tobs", ycol = "voltage",
  t = t0, Delta = delta, h = bw, smooth_ker = epanechnikov)
dt_locreg

# mean function estimation ----

## Define bandwidth grid
N <- dt_real_data[, length(unique(id_curve))]
lambda <- 1440
K <- 40
b0 <- 1 / {(N * lambda) ** (0.9)}
bK <- 0.5
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(K, b0, bK, a) ; gc()

## Get local regularity parameters
Ht <- dt_locreg[order(t), Ht]
Lt <- dt_locreg[order(t), Lt]
bw <- dt_locreg[, unique(locreg_bw)]
t0 <- seq(0.1, 0.99, len = 96)

## Estimate the mean function risk
dt_risk_muhat <- estimate_mean_risk(
  data = dt_real_data, idcol = "id_curve", tcol = 'tobs', ycol = "voltage",
  t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
)
saveRDS(dt_risk_muhat, file.path(script_path, "dt_mean_risk_real_data.RDS"))

## Estimate the mean function
dt_optbw <- dt_risk_muhat[, .("optbw" = h[which.min(mean_risk)]), by = "t"]

dt_muhat <- estimate_mean(
  data = dt_real_data, idcol = "id_curve", tcol = 'tobs', ycol = "voltage",
  t = t0, optbw = dt_optbw[, optbw], Ht = Ht, Lt = Lt, Delta = NULL, h = bw
)
Ht_vec <- Ht
Lt_vec <- Lt
dt_muhat[, "Ht" := Ht_vec]
dt_muhat[, "Lt" := Lt_vec]
dt_muhat[, "locreg_bw" := dt_locreg[order(t), locreg_bw]]
saveRDS(dt_muhat, file.path(script_path, "dt_mean_real_data.RDS"))

# Estimate mean on subset of curve function estimation ----
## Import data and take a subset
id_curve_subset <- tail(1:1358, 50)
dt_real_data_subset <- dt_real_data[id_curve %in% id_curve_subset]

## Estimate local regularity parameters
lambda <- 1440
delta <- exp(-log(lambda) ** (1/3))
bw <- unique(dt_real_data_subset[order(id_curve), .(id_curve, presmooth_bw)])[, presmooth_bw]
t0 <- seq(0.1, 0.99, len = 96)

dt_locreg <- estimate_locreg(
  data = dt_real_data_subset[order(id_curve)], idcol = "id_curve", tcol = "tobs", ycol = "voltage",
  t = t0, Delta = delta, h = bw, smooth_ker = epanechnikov)
dt_locreg

# mean function estimation
N <- dt_real_data_subset[, length(unique(id_curve))]
lambda <- 1440
K <- 40
# b0 <- 1 / lambda
b0 <- 1 / {(N * lambda) ** (0.9)}
# bK <- (N * lambda) ** (- 1 / (2 * max(dt_locreg$Ht) + 1))
# bK <- min(bK, 1 / 10)
bK <- 0.5
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(K, b0, bK, a) ; gc()

Ht <- dt_locreg[order(t), Ht]
Lt <- dt_locreg[order(t), Lt]
bw <- dt_locreg[, unique(locreg_bw)]
t0 <- seq(0.1, 0.99, len = 96)

start_time <- Sys.time()
dt_risk_muhat <- estimate_mean_risk(
  data = dt_real_data_subset, idcol = "id_curve", tcol = 'tobs', ycol = "voltage",
  t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
)
end_time <- Sys.time()
time_taken <- as.numeric(end_time - start_time)
# time_taken <- time_taken + dt_locreg[id_mc == mc_i, unique(time_taken)]
saveRDS(dt_risk_muhat, file.path(script_path, "dt_mean_risk_real_data_subset.RDS"))

# Compare optimal bandwidths
dt_optbw <- dt_risk_muhat[, .("optbw" = h[which.min(mean_risk)]), by = "t"]

dt_muhat <- estimate_mean(
  data = dt_real_data_subset, idcol = "id_curve", tcol = 'tobs', ycol = "voltage",
  t = t0, optbw = dt_optbw[, optbw], Ht = Ht, Lt = Lt, Delta = NULL, h = bw
)
Ht_vec <- Ht
Lt_vec <- Lt
dt_muhat[, "Ht" := Ht_vec]
dt_muhat[, "Lt" := Lt_vec]
dt_muhat[, "locreg_bw" := dt_locreg[order(t), locreg_bw]]
saveRDS(dt_muhat, file.path(script_path, "dt_mean_real_data_subset.RDS"))

