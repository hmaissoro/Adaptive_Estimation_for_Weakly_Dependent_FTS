#' Local Regularity Parameters Estimation
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the local regularity parameters of the underlying process.
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it will be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#' @param center \code{logical}. If \code{TRUE}, the curves are centered.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The points around which the local regularity parameters are estimated.}
#'            \item{locreg_bw :}{ The presmoothing bandwidth.}
#'            \item{Delta :}{ The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.}
#'            \item{Nused :}{ The number of curves that give non-degenerate estimates around \code{t}.}
#'            \item{Ht :}{ The local exponent estimates for each \code{t}. It corresponds to \eqn{H_t}}
#'            \item{Lt :}{ The Hölder constant estimates \code{t}. It corresponds to \eqn{L_t^2}.}
#'         }
#' @export
#'
#' @import data.table
#' @importFrom methods is
#' @importFrom stats median quantile
#'
#' @seealso [estimate_nw()], [estimate_nw_bw()], [simulate_far()], etc.
#'
#' @examples
#' \dontrun{
#'   # Generate a sample of FAR(1)
## Exponent H
#' Hfun <- function(t) {
#'   hurst_logistic(t = t, h_left = 0.4, h_right = 0.8, slope = 5)
#' }
#'
#' ## Hölder constant
#' L <- 4
#'
#' dt_far <- simulate_far(N = 200L, lambda = 100L,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = Hfun,
#'                        L = L,
#'                        far_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
#'                        far_mean = function(t) 4 * sin(1.5 * pi * t),
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Estimate local regularity at
#' t0 <- seq(0.2, 0.8, len = 8)
#'
#' ## If data is a data.table or a data. frame
#' dt_locreg <- estimate_locreg(data = dt_far,
#'                              idcol = "id_curve",
#'                              tcol = "tobs",
#'                              ycol = "X",
#'                              t = t0,
#'                              Delta = NULL,
#'                              h = NULL,
#'                              smooth_ker = epanechnikov)
#' DT::datatable(dt_locreg)
#'
#' ## If data is a list of data.table (or data. frame)
#' list_dt_far <- lapply(unique(dt_far[, id_curve]), function(idx){
#'   dt_far[id_curve == idx, list(tobs, X)]
#' })
#'
#' dt_locreg_2 <- estimate_locreg(data = list_dt_far,
#'                                idcol = NULL,
#'                                tcol = "tobs",
#'                                ycol = "X",
#'                                t = t0,
#'                                Delta = NULL,
#'                                h = NULL,
#'                                smooth_ker = epanechnikov)
#' DT::datatable(dt_locreg_2)
#'
#' ## If data is a list of list
#' list_list_far <- lapply(unique(dt_far[, id_curve]), function(idx){
#'   list("Obs_point" = dt_far[id_curve == idx, tobs],
#'        "Xobs" = dt_far[id_curve == idx, X])
#' })
#'
#' dt_locreg_3 <- estimate_locreg(data = list_list_far,
#'                                idcol = NULL,
#'                                tcol = "Obs_point",
#'                                ycol = "Xobs",
#'                                t = t0,
#'                                Delta = NULL,
#'                                h = NULL,
#'                                smooth_ker = epanechnikov)
#' DT::datatable(dt_locreg_2)
#'
#'}
#'
#'
#'
estimate_locreg <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                            t = 1/2, Delta = NULL, h = NULL,
                            smooth_ker = epanechnikov, center = TRUE){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (! methods::is(center, "logical"))
    stop("'center' must be a TRUE or FALSE.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Control and set arguments depending on the data
  if (! is.null(Delta)) {
    if (! (methods::is(Delta, "numeric") & data.table::between(Delta, 0, 1) & length(Delta) == 1))
      stop("'Delta' must be a numeric scalar value between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    Delta <- 2 * exp(-log(lambdahat) ** 0.72)
  }

  # Control on the pre-smoothing bandwidth
  if (! is.null(h)) {
    # h is a vector or a scalar
    if (! all(methods::is(h, "numeric") & data.table::between(h, 0, 1))){
      stop("'h' must be a numeric vector or scalar value(s) between 0 and 1.")
    } else if (length(h) > 1 & length(h) != N) {
      stop("If 'h' is given as a vector, its length must be equal to the number of curves in 'data'.")
    }
  } else {
    # If h = NULL, choose the bandwidth by CV
    if (N > 50) {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, smooth_ker = smooth_ker)
      h <- dt_optbw[, median(optbw)]
      rm(dt_optbw) ; gc()
    } else {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, smooth_ker = smooth_ker)
      h <- dt_optbw[, optbw]
      rm(dt_optbw) ; gc()
    }
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Step 1: estimate the curve at
  ## Take into account the cases where t-Delta/2 < 0 and t + Delta/2
  dt_t <- data.table::rbindlist(lapply(t, function(ti, D){
    if ((ti - D / 2 ) <= 0) {
      t1 <- ti
      t2 <- ti + D / 2
      t3 <- ti + D
    } else if ((ti + D / 2 ) >= 1) {
      t3 <- ti
      t2 <- ti - D / 2
      t1 <- ti - D
    } else {
      t1 <- ti - Delta / 2
      t2 <- ti
      t3 <- ti + Delta / 2
    }
    return(data.table::data.table("t1" = t1, "t2" = t2, "t3" = t3))
  }, D = Delta))
  t1 <- dt_t[, t1]
  t2 <- dt_t[, t2]
  t3 <- dt_t[, t3]
  rm(dt_t)

  dt_smooth <- data.table::rbindlist(lapply(data[, unique(id_curve)], function(i, data, h, t1, t2, t3){
    ## smooth an estimate
    dt1 <- estimate_nw(y = data[id_curve == i, X],
                       t = data[id_curve == i, tobs],
                       tnew = t1, h = h[i],
                       smooth_ker = smooth_ker)
    dt2 <- estimate_nw(y = data[id_curve == i, X],
                       t = data[id_curve == i, tobs],
                       tnew = t2, h = h[i],
                       smooth_ker = smooth_ker)
    dt3 <- estimate_nw(y = data[id_curve == i, X],
                       t = data[id_curve == i, tobs],
                       tnew = t3, h = h[i],
                       smooth_ker = smooth_ker)

    dt_out <- data.table::data.table(
      id_curve = i, t1 = t1, xt1 = dt1$yhat,
      t2 = t2, xt2 = dt2$yhat,
      t3 = t3, xt3 = dt3$yhat
    )
    # rm(xtilde, inKernelSupp)

    return(dt_out)
  }, h = h, data = data, t1 = t1, t2 = t2, t3 = t3))

  # Step 2 : estimate regularity parameters

  dt_reg <- data.table::rbindlist(lapply(1:length(t2), function(i, dt_smooth, t1, t2, t3, t, center){
    ## Extract X_1(g),...,X_N(g) where g = t1, t2 or t3
    xt1 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt1]
    xt2 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt2]
    xt3 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt3]

    ## Remove NA values
    any_na <- (is.na(xt1) | is.na(xt2) | is.na(xt3))
    xt1 <- xt1[! any_na]
    xt2 <- xt2[! any_na]
    xt3 <- xt3[! any_na]
    rm(any_na)
    ## Remove NaN values
    any_nan <- (is.nan(xt1) | is.nan(xt2) | is.nan(xt3))
    xt1 <- xt1[! any_nan]
    xt2 <- xt2[! any_nan]
    xt3 <- xt3[! any_nan]
    rm(any_nan)

    ## Remove extreme values
    rxt1 <- (xt1 >= quantile(xt1, 0.025, na.rm = TRUE, type = 1)) & (xt1 <= quantile(xt1, 0.975, na.rm = TRUE, type = 1))
    rxt2 <- (xt2 >= quantile(xt2, 0.025, na.rm = TRUE, type = 1)) & (xt2 <= quantile(xt2, 0.975, na.rm = TRUE, type = 1))
    rxt3 <- (xt3 >= quantile(xt3, 0.025, na.rm = TRUE, type = 1)) & (xt3 <= quantile(xt3, 0.975, na.rm = TRUE, type = 1))
    xt1 <- xt1[rxt1 & rxt2 & rxt3]
    xt2 <- xt2[rxt1 & rxt2 & rxt3]
    xt3 <- xt3[rxt1 & rxt2 & rxt3]
    Nused <- sum(rxt1 & rxt2 & rxt3)
    rm(rxt1, rxt2, rxt3)

    ## Center data if the argument center = TRUE
    xt1 <- xt1 - center * mean(xt1)
    xt2 <- xt2 - center * mean(xt2)
    xt3 <- xt3 - center * mean(xt3)

    ## Compute Unweighed local regularity parameters
    theta_t1_t3 <- mean((xt1 - xt3) ** 2, na.rm = TRUE)
    theta_t1_t2 <- mean((xt1 - xt2) ** 2, na.rm = TRUE)
    theta_t2_t3 <- mean((xt2 - xt3) ** 2, na.rm = TRUE)
    Ht <- (log(theta_t1_t3) - log(theta_t2_t3))  / (2 * log(2))
    Lt <- theta_t2_t3 / (abs(t2[i] - t3[i]) ** (2 * Ht))

    ## Return the result
    dt_out <- data.table(t = t[i], Ht = Ht, Lt = Lt, Nused = Nused)
    rm(theta_t1_t3, theta_t2_t3, theta_t1_t2, Ht, Lt, xt1, xt2, xt3)

    return(dt_out)
  }, dt_smooth = dt_smooth, t1 = t1, t2 = t2, t3 = t3, t = t, center = center))

  dt_reg[, c("locreg_bw", "Delta") := list(median(h), Delta)]

  data.table::setcolorder(x = dt_reg, neworder = c("t", "Delta", "Nused", "locreg_bw", "Ht", "Lt"))

  return(dt_reg)
}

#' Format data for local regularity estimation
#'
#' @param data \code{data.table (or data.frame)} or \code{list} of \code{data.table (or data.frame)} or \code{list} of \code{list}.
#' \itemize{
#'    \item{If \code{data.table}}{
#'        It must contain the raw binding of the curve observations with at least 3 columns.
#'        \itemize{
#'          \item{\code{idcol} :}{ The name of the column that contains the index of the curve in the sample.
#'                              Each index of a curve is repeated as many times as it has observation points.}
#'          \item{\code{tcol} :}{ The name of the column that contains the observation points associated to each curve index.}
#'          \item{\code{ycol} :}{ The name of the column that contains the observed value of a curve at each point of observation and for each index of the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{data.table}}{
#'         In this case, each element of the given \code{list} corresponds to the observation scheme of a curve, which is given as \code{data.table} or \code{data.frame}.
#'         The data.table contains at least 2 columns.
#'         \itemize{
#'          \item{\code{tcol} :}{ The name of the column that contains the observation points associated to the curve.}
#'          \item{\code{ycol} :}{ The name of the column that contains the observed value of the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{list}}{
#'      In the latter case, the \code{data} is a list \code{list} where each element is the observation scheme of a curve given as a \code{list} of 2 vectors.
#'      \itemize{
#'          \item{\code{tcol} :}{ The name of the vector that contains the observation points associated the curve.}
#'          \item{\code{ycol} :}{ The name of the vector that contains the observed value of the curve.}
#'        }
#'    }
#' }
#' @param idcol \code{character}. If \code{data} is given as \code{data.table} or \code{data.frame},
#' it is the name of the column that contains the index of the curve in the sample.
#' Each index of a curve is repeated as many times as it has observation points.
#' Opposite, if f \code{data} is given as \code{list} of \code{data.table (or data.frame)} of \code{list} of \code{list}, \code{idcol = NULL.}
#' @param tcol \code{character}. The name of the column (or vector) that contains the observation points associated to the curves.
#' @param ycol \code{character}. The name of the column that contains the observed value of the curves.
#'
#' @return A \code{data.table} containing 3 columns.
#'          \itemize{
#'            \item{id_curve :}{ The index of the curve.}
#'            \item{tobs :}{ The observation points associated to each curve \code{id_curve}.}
#'            \item{X :}{ The observed values of the curve associated to \code{id_curve} at \code{tobs} observation points.}
#'         }
#'
#' @import data.table
#' @importFrom methods is
#'
.format_data <- function(data, idcol = NULL, tcol = "tobs", ycol = "X"){
  # Check if data is a data.table or data.frame
  is_dt_or_df <- methods::is(data, "data.table") | methods::is(data, "data.frame")

  # Check if data is a list of data.table (or data.frame)
  is_list_dt_or_df <- methods::is(data, "list") &
    all(unlist(lapply(data, function(element){
      methods::is(element,"data.table") | methods::is(element,"data.frame")
    })))

  # Check if data is a list of list
  is_list_of_list <- methods::is(data, "list") &
    all(unlist(lapply(data, function(element){
      methods::is(element,"list")
    })))

  if (! (is_dt_or_df | is_list_dt_or_df | is_list_of_list))
    stop("'data' must of class data.table (or data.frame) or a list of data.table (or data.table) or a list of list.")

  if (is_dt_or_df) {
    if (is.null(idcol))
      stop("If the class of 'data' is data.table (or data.frame), 'idcol' need to be specifyed.")
    if (! all(c(idcol, tcol, ycol) %in% colnames(data))){
      stop("The specified column name 'idcol' or 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::as.data.table(data)
      data <- data[, .SD, .SDcols = c(idcol, tcol, ycol)]
      names(data) <- c("id_curve", "tobs", "X")
      data <- data[, list(id_curve, tobs, X)]
      Mn <- data[, .N, by = id_curve][, N]
      N <- length(Mn)
      id <- unlist(lapply(1:N, function(n, Mn){
        rep(n, Mn[n])
      }, Mn = Mn))
      data[, id_curve := id]
      rm(Mn, N, id)
    }
  } else if (is_list_dt_or_df) {
    if (! is.null(idcol))
      stop("If 'data' is a list of data.table (or data.table) or a list of list, 'idcol' must be NULL.")
    check_colname <- all(unlist(lapply(data, function(element, tcol, ycol){
      all(c(tcol, ycol) %in% colnames(element))
    }, tcol = tcol, ycol = ycol)))
    if (! check_colname) {
      stop("The specified column name 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::rbindlist(lapply(1:length(data), function(i, tcol, ycol){
        data[[i]] <- data.table::as.data.table(data[[i]])
        data.table::setnames(x = data[[i]], old = c(tcol, ycol), new = c("tobs", "X"))
        dt <- data.table::data.table("id_curve" = i, data[[i]][, list(tobs, X)])
      }, tcol = tcol, ycol = ycol))
    }

  } else if (is_list_of_list) {
    if (! is.null(idcol))
      stop("If 'data' is a list of data.table (or data.table) or a list of list, 'idcol' must be NULL.")
    check_vecname <- all(unlist(lapply(data, function(element, tcol, ycol){
      all(c(tcol, ycol) %in% names(element))
    }, tcol = tcol, ycol = ycol)))
    if (! check_vecname) {
      stop("The specified vector name 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::rbindlist(lapply(1:length(data), function(i, tcol, ycol){
        data[[i]] <- data.table::as.data.table(data[[i]])
        data.table::setnames(x = data[[i]], old = c(tcol, ycol), new = c("tobs", "X"))
        dt <- data.table::data.table("id_curve" = i, data[[i]][, list(tobs, X)])
      }, tcol = tcol, ycol = ycol))
    }
  } else {
    NA
  }
  return(data)
}
