#' Estimate the risk of the mean function
#'
#' The risk of the mean function is the function \eqn{R_\mu(t;h)} in
#' Section 4.1 of \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each \code{t}.
#' It can be \code{NULL} and that way it will be defined as an exponential grid of \eqn{N\times\lambda}.
#' @param Ht \code{vector (numeric)}. The estimates of the local exponent for each \code{t}.
#' Default \code{Ht = NULL} and thus it will be estimated.
#' @param Lt \code{vector (numeric)}. The estimates of the Hölder constant for each \code{t}.
#' It corresponds to \eqn{L_t^2}. Default \code{Lt = NULL} and thus it will be estimated.
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it will be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The points at which the risk function is estimated.}
#'            \item{h :}{ The candidate bandwidth.}
#'            \item{Ht :}{ The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}}
#'            \item{Lt :}{ The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
#'            \item{locreg_bw :}{ The bandwidth used to estimate the local regularity parameters.}
#'            \item{bias_term :}{ The bias term of the risk function.}
#'            \item{varriance_term :}{ The variance term of the risk function.}
#'            \item{dependence_term :}{ The dependence term of the risk function.}
#'            \item{mean_risk :}{ The estimates of the risk function of the mean.}
#'         }
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_autocov()].
#'
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist between
#' @importFrom Rdpack
#'
#' @references
#'  \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' # Generate a sample path of FTS
#' dt_far <- simulate_far(N = 50, lambda = 70,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        far_kernel = get_real_data_far_kenel,
#'                        far_mean = get_real_data_mean,
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Add noise
#' dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]
#'
#' # Estimate risk function
#' dt_mean_risk <- estimate_mean_risk(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
#'   Delta = NULL, h = NULL, smooth_ker = epanechnikov)
#'
#' # Plot mean risk
#' dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
#' manipulateWidget::combineWidgets(
#'   list = list(
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
#'       main = "t = 0.25", xlab = "h", ylab = "risk function"),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
#'       main = "t = 0.5", xlab = "h", ylab = "risk function"),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
#'       main = "t = 0.75", xlab = "h", ylab = "risk function")
#'   ),
#'   nrow = 3
#' )
#'
#' }
#'
#'
estimate_mean_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                               t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
                               Ht = NULL, Lt = NULL, Delta = NULL, h = NULL,
                               smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  # Control on local regularity parameters
  if ( ((!is.null(Ht)) & length(Ht) != length(t)) | ((!is.null(Lt)) & length(Lt) != length(t)))
    stop("If 'Ht' or 'Lt' is not NULL, it must be the same length as 't'.")
  if ((!is.null(Ht) & ! methods::is(Ht, "numeric")) |
      (!is.null(Lt) & ! methods::is(Lt, "numeric")))
    stop("If 'Ht' or 'Lt' is not NULL, it must be numeric.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (! is.null(bw_grid)) {
    if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
      stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    K <- 20
    b0 <- 4 * (N * lambdahat) ** (- 0.9)
    bK <- 4 * (N * lambdahat) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a, lambdahat) ; gc()
  }

  # Estimate local regularity parameters
  # This function controls the remaining arguments
  if (is.null(Ht) | is.null(Lt)) {
    dt_locreg <- estimate_locreg(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = t, Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    Ht <- dt_locreg[, Ht]
    Lt <- dt_locreg[, Lt]
    ht <- dt_locreg[, unique(locreg_bw)]
  } else {
    ht <- median(h)
  }

  # Estimation of the observation error standard deviation
  dt_sigma <- estimate_sigma(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = t)

  # Smooth curve with presmoothing bandwidth
  # Smooth curves
  # dt_smooth <- data.table::rbindlist(lapply(1:N, function(curve_index, hvec, t, kernel_smooth, data){
  #   x_smooth <- estimate_nw(y = data[id_curve == curve_index, X],
  #                           t = data[id_curve == curve_index, tobs],
  #                           tnew = t,
  #                           h = hvec[curve_index],
  #                           smooth_ker = kernel_smooth)
  #   xhat <- x_smooth$yhat
  #   dt_res <- data.table::data.table("id_curve" = curve_index, "t" = t, "x" = xhat)
  #   return(dt_res)
  # }, hvec = rep(ht, N), t = t, kernel_smooth = smooth_ker, data = data))

  # Estimation of the empirical autocovariance using the presmoothing bandwidth
  dt_autocov <- estimate_empirical_autocov(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = t, lag = 0:(N-1),
    h = ht, smooth_ker = smooth_ker)

  # If the variance is NaN, set it as 0
  dt_autocov[is.nan(autocov) & lag == 0, autocov := 0]

  # Estimate the risk function
  dt_mean_risk <- data.table::rbindlist(lapply(bw_grid, function(h, t, Ht, Lt, presmooth_bw, kernel_smooth, data, sig_error, N, dt_autocov){
    # compute quantities to estimate estimate the risk
    dt_risk <- data.table::rbindlist(lapply(1:N, function(curve_index, t, h, Ht, kernel_smooth, data){
      # Compute the weight of the estimator
      Tn <- data[id_curve == curve_index, tobs]
      ker <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, Ker) Ker((Tnm - ti) / h), Ker = kernel_smooth)
      ## By convention NaN = 0
      wker <- apply(ker, 1, function(r) r / ifelse(sum(r) != 0, sum(r), 1))
      wker <- t(wker)

      # Take the maximum and c_n(t,h), the sum of the weight
      wmax <- apply(wker, 1, max)
      cn <- apply(wker, 1, function(r) sum(abs(r)))

      Tn_t_2H <- outer(
        X = 1:length(t), Y = Tn,
        FUN = function(tid, Tnm, t, Ht) abs((Tnm - t[tid]) / h) ** (2 * Ht[tid]),
        t = t, Ht = Ht
      )

      bn2H <- diag(abs(wker) %*% t(Tn_t_2H))

      # \pi_n(t;h)
      pi_n <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, h) as.numeric(abs(Tnm - ti) <= h), h = h)
      pi_n <- as.numeric(rowSums(pi_n) >= 1)

      # data.table return
      dt_res <- data.table::data.table(
        "id_curve" = curve_index, "t" = t, "wmax" = wmax,
        "cn" = cn, "bn2H" = bn2H, "pi_n" = pi_n
      )
      return(dt_res)
    }, h = h, t = t, Ht = Ht, kernel_smooth = kernel_smooth, data = data))

    # Compute P_N(t, h)
    dt_risk[, PN := sum(pi_n), by = t]

    # compute \mathbb B(t,h, 2H) for each n = 1, ..., N
    dt_risk[, B := sum(pi_n * cn * bn2H) / PN, by = t]

    # Compute \mathbb V_mu(t,h) for each n = 1, ..., N
    dt_risk[, Vmu := sum(pi_n * cn * wmax) / (PN ** 2), by = t]

    # Compute the risk for each h
    dt_rk <- unique(dt_risk[, list(t, B, Vmu, PN)])

    ## Biais term
    bias_term <- 2 * Lt * (h ** (2 * Ht)) * dt_rk[, B]

    ## Variance term
    varriance_term <- 2 * (sig_error ** 2) * dt_rk[, Vmu]

    ## Dependence term
    ### Compute \rho_\ell(t;h)
    dt_rho_ell <- unique(dt_risk[, list(id_curve, t, pi_n, PN)])
    dt_rho <- data.table::rbindlist(lapply(1:(N-1), function(ell, t, dt_rho_ell){
      data.table::rbindlist(lapply(t, function(ti, ell, dt_rho_ell){
        PN <- unique(dt_rho_ell[t == ti, PN])
        pi_vec <- dt_rho_ell[t == ti][order(id_curve), pi_n]
        pi_i <- pi_vec[1:(N - ell)]
        pi_i_plus_ell <- pi_vec[(1 + ell):N]
        rho_ell <- sum(pi_i * pi_i_plus_ell) / PN
        dt_res <- data.table::data.table("t" = ti, "lag" = ell, "rho" = rho_ell)
        return(dt_res)
      }, ell = ell, dt_rho_ell = dt_rho_ell))
    }, t = t, dt_rho_ell =dt_rho_ell))

    dt_lr_var <- data.table::merge.data.table(
      x = dt_autocov[lag > 0],
      y = dt_rho,
      by = c("t", "lag")
    )
    dt_lr_var <- dt_lr_var[!is.nan(autocov), list("lr_var" = sum(2 * autocov * rho)), by = t]
    dependence_coef <- dt_autocov[lag == 0][order(t), autocov] + dt_lr_var[order(t), lr_var]
    ### Note that dependence_coef <= abs(dependence_coef), thus
    dependence_coef <- abs(dependence_coef)
    dependence_term <- 2 * dependence_coef  / dt_rk[, PN]

<<<<<<< HEAD

=======
>>>>>>> d8b5dd12dee2ab3dddb7b7a432b3c5de45fbdff4
    ## Final risk function
    mean_risk <- bias_term + varriance_term + dependence_term

    # Result to returned
    dt_res <- data.table::data.table("t" = t, "h" = h, "Ht" = Ht, "Lt" = Lt, "locreg_bw" = presmooth_bw,
                                     "bias_term" = bias_term, "varriance_term" = varriance_term,
                                     "dependence_term" = dependence_term, "mean_risk" = mean_risk)
    return(dt_res)
  }, t = t, Ht = Ht, Lt = Lt, presmooth_bw = ht, kernel_smooth = smooth_ker, data = data,
  sig_error = dt_sigma[, sig], N = N, dt_autocov = dt_autocov))

  return(dt_mean_risk)
}


#' Estimate mean function
#'
#' Mean function estimation using the adaptive estimator of \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
#'
#' @inheritParams estimate_mean_risk
#' @param optbw \code{vector (numeric)}. The optimal bandwidth parameter for mean function estimation for each \code{t}.
#' Default \code{optbw = NULL} and thus it will be estimated using \link{estimate_mean_risk} function.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The points at which the risk function is estimated.}
#'            \item{locreg_bw :}{ The bandwidth used to estimate the local regularity parameters.}
#'            \item{Ht :}{ The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}}
#'            \item{Lt :}{ The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
#'            \item{optbw :}{ The optimal bandwidth. That is the bandwidth which minimises the risk function.}
#'            \item{PN :}{ Number of selected curves for each \code{t}.}
#'            \item{muhat :}{ The estimates of the mean function.}
#'         }
#' @export
#'
#' @seealso [estimate_mean_risk()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_autocov()].
#'
#' @importFrom data.table data.table rbindlist setcolorder merge.data.table
#' @importFrom Rdpack reprompt
#'
#' @references
#' insertAllCited{}
#'
#' @examples
#' \dontrun{
#' # Generate a FAR A process
#' dt_far <- simulate_far(N = 50, lambda = 70,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        far_kernel = get_real_data_far_kenel,
#'                        far_mean = get_real_data_mean,
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Add noise
#' dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]
#'
#' # Estimate mean function
#' dt_mean <- estimate_mean(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
#'   Delta = NULL, h = NULL, smooth_ker = epanechnikov)
#'
#' # Table of the estimates of the mean function
#' DT::datatable(data = dt_mean[, lapply(.SD, function(X) round(X, 3))])
#' }
#'
#'
estimate_mean <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          t = c(1/4, 1/2, 3/4), optbw = NULL, bw_grid = seq(0.005, 0.15, len = 45),
                          Ht = NULL, Lt = NULL, Delta = NULL, h = NULL,
                          smooth_ker = epanechnikov){
  # Control on t, optbw and smooth_ker arguments
  # NB : The remaining arguments are controlled using the .format_data and estimate_mean_risk functions, if required.
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if ((!is.null(optbw)) & length(optbw) != length(t)) {
    stop("If 'optbw' is not NULL, it must be the same length as 't'.")
  } else if (is.null(optbw)) {
    if ((! is.null(bw_grid)) ) {
      if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
        stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
    } else {
      lambdahat <- mean(data[, .N, by = "id_curve"][, N])
      K <- 20
      b0 <- 4 * (N * lambdahat) ** (- 0.9)
      bK <- 4 * (N * lambdahat) ** (- 1 / 3)
      a <- exp((log(bK) - log(b0)) / K)
      bw_grid <- b0 * a ** (seq_len(K))
      rm(K, b0, bK, a) ; gc()
    }
  }

  if (is.null(optbw)) {
    # Estimate the risk function
    dt_mean_risk <- estimate_mean_risk(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = t, bw_grid = bw_grid,
      Ht = Ht, Lt = Lt, Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    # Take the optimum of the risk function
    dt_mean_risk[, optbw := h[which.min(mean_risk)], by = t]
    dt_mean_optbw <- unique(dt_mean_risk[, list(t, Ht, Lt, locreg_bw, optbw)])
  }else {
    dt_mean_optbw <- data.table::data.table("t" = t, "Ht" = NA, "Lt" = NA, "locreg_bw" = NA, "optbw" = optbw)
  }

  # Smooth curves with optimal bandwidth parameters
  dt_Xhat <- data.table::rbindlist(lapply(1:N, function(curve_index, t, data, dt_mean_optbw, smooth_ker){
    # \pi_n(t,h)
    Tn <- data[id_curve == curve_index, tobs]
    Yn <- data[id_curve == curve_index, X]
    pi_n <- sapply(X = 1:length(t), function(tidx, Tn, t, optbw_mean){
      as.numeric(abs(Tn - t[tidx]) <= optbw_mean[tidx])
    }, Tn = Tn, t = t, optbw_mean = dt_mean_optbw[order(t), optbw])
    pi_n <- t(pi_n)
    pi_n <- as.numeric(rowSums(pi_n, na.rm = TRUE) >= 1)

    # \widehat X(t;h)
    Xhat <- mapply(function(t, h, Yn, Tn, ker){
      estimate_nw(y = Yn, t = Tn, tnew = t, h = h, smooth_ker = ker)$yhat
    }, t = dt_mean_optbw[, t], h = dt_mean_optbw[, optbw],
    MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

    dt_res <- data.table::data.table("id_curve" = curve_index, "t" = t, "Xhat" = Xhat, "pi_n" = pi_n)
    return(dt_res)
  }, data = data, t = t, dt_mean_optbw = dt_mean_optbw, smooth_ker = smooth_ker))

  # Estimate mean function \widehat \mu(t)
  dt_Xhat[is.nan(Xhat) & pi_n == 0, Xhat := 0]
  dt_Xhat[, PN := sum(pi_n), by = t]
  dt_Xhat <- dt_Xhat[!is.nan(Xhat)]
  dt_Xhat[, "muhat" := sum(pi_n * Xhat) / PN, by = t]

  # Add the optimal bandwidth
  dt_muhat <- unique(dt_Xhat[, list(t, muhat, PN)])
  dt_muhat <- data.table::merge.data.table(x = dt_muhat, y = dt_mean_optbw, by = "t")
  data.table::setcolorder(x = dt_muhat, neworder = c("t", "locreg_bw", "Ht", "Lt", "optbw", "PN", "muhat"))
  return(dt_muhat)
}

# mean function estimator : Rubìn et Paranaretos ----
#' Estimate mean function using \insertCite{rubin2020;textual}{adaptiveFTS} method.
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @import data.table data.table
#' @importFrom fastmatrix kronecker.prod
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The Observation points at which the mean function is estimated.}
#'            \item{h :}{ The bandwidth parameter.}
#'            \item{muhat_RP :}{ The estimates of the mean function using Rubìn and Panaretos (2020) method.}
#'         }
#' @export
#'
#' @seealso [estimate_mean_bw_rp()]
#' @importFrom Rdpack reprompt
#' @references
#' insertAllcited{}
#'
#' @examples
#' \dontrun{
#' # Generate a FAR A process
#' dt_far <- simulate_far(N = 50, lambda = 70,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        far_kernel = get_real_data_far_kenel,
#'                        far_mean = get_real_data_mean,
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Add noise
#' dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]
#'
#' # Estimate mean function using Rubìn and Panaretos (2020) method
#' dt_mean_rp <- estimate_mean_rp(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), h = 5/70, smooth_ker = epanechnikov)
#'
#' DT::datatable(data = dt_mean_rp[, lapply(.SD, function(X) round(X, 5))])
#'
#' }
#'
estimate_mean_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             t = c(1/4, 1/2, 3/4), h, smooth_ker = epanechnikov){
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]
  lambdahat <- mean(data[, .N, by = "id_curve"][, N])

  if (N <= 250 && lambdahat < 200) {
    # Estimate mean function
    Tn <- data[order(tobs), tobs]
    Yn <- data[order(tobs), X]
    data_curve <- fastmatrix::kronecker.prod(
      x = matrix(data = rep(1, length(t)), ncol = 1),
      y = cbind(Tn, Yn)
    )
    colnames(data_curve) <- c("Tn", "Yn")
    data_curve <- data.table::as.data.table(data_curve)
    tvec <- rep(t, each = length(Tn))
    data_curve[, t := tvec]
    rm(Yn, Tn, tvec, data) ; gc() ; gc()

    data_curve[, Tn_minus_t := (Tn - t)]
    data_curve[, Tn_minus_t_over_h := Tn_minus_t / h]

    # Compute mean Q and S function
    dt_res_by_curve <- data_curve[
      ,
      .("Q0" = sum((Tn_minus_t ** 0) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "Q1" = sum((Tn_minus_t ** 1) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "S0" = sum((Tn_minus_t ** 0) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "S1" = sum((Tn_minus_t ** 1) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "S2" = sum((Tn_minus_t ** 2) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N),
      by = "t"
    ]
    rm(data_curve) ; gc() ; gc()

    # Estimate mean
    dt_res <- dt_res_by_curve[, .("muhat_RP" = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)), by = "t"]
    dt_res[, "h" := h]
    data.table::setcolorder(x = dt_res, neworder = c("t", "h", "muhat_RP"))
    rm(dt_res_by_curve) ; gc() ; gc()

  } else {
    # Split t if N x \lambda >> 0
    if (N <= 450 && lambdahat <= 300) {
      N_t_by_list <- 40 * 300 / 10
      t_list <- split(t, ceiling(seq_along(t) / N_t_by_list))
    } else if (N <= 1000 & lambdahat <= 50){
      N_t_by_list <- 1000 * 40 / 50
      t_list <- split(t, ceiling(seq_along(t) / N_t_by_list))
    } else {
      t_list <- t
    }

    # Estimate mean function
    dt_res <- data.table::rbindlist(lapply(t_list, function(t_list_i, data, N, h){
      Tn <- data[order(tobs), tobs]
      Yn <- data[order(tobs), X]
      if (length(t_list_i) > 1) {
        data_curve <- fastmatrix::kronecker.prod(
          x = matrix(data = rep(1, length(t_list_i)), ncol = 1),
          y = cbind(Tn, Yn)
        )
        colnames(data_curve) <- c("Tn", "Yn")
        data_curve <- data.table::as.data.table(data_curve)
        tvec <- rep(t_list_i, each = length(Tn))
        data_curve[, t := tvec]
        rm(Yn, Tn, tvec) ; gc() ; gc()
      } else {
        data_curve <- data.table::data.table("Tn" = Tn, "Yn" = Yn, "t" = t_list_i)
      }

      data_curve[, Tn_minus_t := (Tn - t)]
      data_curve[, Tn_minus_t_over_h := Tn_minus_t / h]

      # Compute mean Q and S function
      dt_res_by_t_list_i <- data_curve[
        ,
        .("Q0" = sum((Tn_minus_t ** 0) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "Q1" = sum((Tn_minus_t ** 1) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "S0" = sum((Tn_minus_t ** 0) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "S1" = sum((Tn_minus_t ** 1) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "S2" = sum((Tn_minus_t ** 2) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N),
        by = "t"
      ]
      rm(data_curve) ; gc() ; gc()
      # Estimate mean
      dt_res_by_t_list_i <- dt_res_by_t_list_i[, .("muhat_RP" = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)), by = "t"]
      dt_res_by_t_list_i[, "h" := h]
      data.table::setcolorder(x = dt_res_by_t_list_i, neworder = c("t", "h", "muhat_RP"))

      # Return
      return(dt_res_by_t_list_i)
    }, data = data, N = N, h = h))
  }
  return(dt_res)
}

#' Bandwidth estimation using cross-validation for the \insertCite{rubin2020;textual}{adaptiveFTS} mean function estimator.
#'
#' @inheritParams .format_data
#' @param Kfold \code{integer (positive)}. Number of fold for the cross-validation.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{h :}{ The candidate bandwidth.}
#'            \item{cv_error :}{ The estimates of the Cross-Validation error for each \code{h}.}
#'         }
#' @export
#' @importFrom caret createFolds
#' @importFrom data.table data.table rbindlist
#' @seealso [estimate_mean_rp()]
#'
#' @importFrom Rdpack reprompt
#' @references
#' insertAllcited{}
#'
#' @examples
#' \dontrun{
#' # Generate a FAR A process
#' dt_far <- simulate_far(N = 50, lambda = 70,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        far_kernel = get_real_data_far_kenel,
#'                        far_mean = get_real_data_mean,
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Add noise
#' dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]
#'
#' ## Estimate the bandwidth by Cross-Validation
#' dt_bw_mean_rp <- estimate_mean_bw_rp(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
#'   smooth_ker = epanechnikov)
#'
#' ## Plot the Cross-Validation error
#' dygraphs::dygraph(dt_bw_mean_rp)
#'
#' ## Select the best bandwidth
#' optbw <- dt_bw_mean_rp[, h[which.min(cv_error)]]
#'
#' ## Estimate the mean function
#' dt_mean_rp <- estimate_mean_rp(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), h = optbw, smooth_ker = epanechnikov)
#'
#' DT::datatable(data = dt_mean_rp[, lapply(.SD, function(X) round(X, 5))])
#'
#' }
#'
estimate_mean_bw_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
                                smooth_ker = epanechnikov){
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  # Create Kfold folds
  fold <- caret::createFolds(y = unique(data[, id_curve]), k = Kfold, list = TRUE)

  # Get risk for each bandwidth in the grid
  dt_bw <- data.table::rbindlist(parallel::mclapply(bw_grid, function(Bmu0, data, fold, kernel_smooth){

    # Compute the cross-validation error for each f in fold
    err_fold <- tryCatch(
      expr = sapply(fold, function(f, data, Bmu0, kernel_smooth){
        # split train - test
        dt_test <- data[id_curve %in% unlist(f)]
        dt_test <- dt_test[order(tobs)]
        dt_train <- data[id_curve %in% setdiff(unlist(fold), unlist(f))]
        dt_train <- dt_train[order(tobs)]

        # Estimation of mean on fold\f and test on f
        dt_mu <- estimate_mean_rp(
          data = dt_train, idcol = "id_curve", tcol = "tobs", ycol = "X",
          t = dt_test[, tobs], h = Bmu0, smooth_ker = kernel_smooth)

        Sqerror <- (dt_test[, X] - dt_mu[, muhat_RP]) ** 2
        err <- sum(Sqerror)
        return(err)
      }, data = data, Bmu0 = Bmu0, kernel_smooth = kernel_smooth, simplify = TRUE),
      error = function(e){
        message("Error in estimating the mean function:")
        print(e)
        return(NA)

      })

    # Cross-validaiton error
    cv_err <- mean(err_fold[!is.nan(err_fold)], na.rm = TRUE)

    # Return the result
    dt_res <- data.table::data.table("h" = Bmu0, "cv_error" = cv_err)
    return(dt_res)

  }, mc.cores = 20, data = data, fold = fold, kernel_smooth = smooth_ker))
  rm(data, fold) ; gc() ; gc()

  return(dt_bw)
}
