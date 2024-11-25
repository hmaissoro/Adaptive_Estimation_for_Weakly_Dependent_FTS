#' Arctan Hurst function
#'
#' Arctan Hurst function that can be used to generate multifractional Brownian motion (mfBm).
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to evaluate the function.
#'
#' @return A \code{vector (float)} corresponding to the value of the function evaluated at \code{t}.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom data.table between
#'
#' @seealso [hurst_linear()], [hurst_logistic()].
#'
#' @examples
#'
#' t0 <- seq(0.2, 0.8, len = 10)
#' htan <- hurst_arctan(t = t0)
#' plot(x = t0, y = htan, type = "b", col = "red")
#'
#'
hurst_arctan <- function(t = seq(0.2, 0.8, len = 10)){
  if (! methods::is(t, "numeric") && all(t >= 0 & t <= 1))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  hval <- atan(t) / pi + 1/2
  return(hval)
}

#' Linear Hurst function
#'
#' Linear Hurst function that can be used to generate multifractional Brownian motion (mfBm).
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to evaluate the function.
#' @param h_left \code{Float}. A scalar value in the interval between 0 and 1 indicating the minimum of the function.
#' @param h_right \code{Float}. A scalar value in the interval between 0 and 1 indicating the maximum of the function.
#'
#' @return A \code{vector (float)} corresponding to the value of the function evaluated at \code{t}.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom data.table between
#'
#' @seealso [hurst_arctan()], [hurst_logistic()].
#'
#' @examples
#' t0 <- seq(0.2, 0.8, len = 10)
#' hlinear <- hurst_linear(t = t0)
#' plot(x = t0, y = hlinear, type = "b", col = "red")
#'
#'
hurst_linear <- function(t = seq(0.2, 0.8, len = 10), h_left = 0.2, h_right = 0.8) {

  if (! (methods::is(t, "numeric") && all(t >= 0 & t <= 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  if (!(methods::is(h_left, "numeric") && methods::is(h_right, "numeric") &&
        (h_left > 0 && h_left < 1) && (h_right > 0 && h_right < 1))) {
    stop("'h_left' and 'h_right' must be scalar values between 0 and 1")
  }

  t1 <- 1
  t0 <- 0
  a <- (h_right - h_left) / (t1 - t0)
  b <- h_right - a * t1
  hval <- pmin(a * t + b, 1)
  return(hval)
}

#' Logistic Hurst function
#'
#' Logistic Hurst function that can be used to generate multifractional Brownian motion (mfBm).
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to evaluate the function.
#' @param h_left \code{Float}. A scalar value in the interval between 0 and 1 indicating the minimum of the function.
#' @param h_right \code{Float}. A scalar value in the interval between 0 and 1 indicating the maximum of the function.
#' @param slope \code{Float (positive)}. A scalar positive value corresponding to the slope of the logistic function.
#' @param change_point_position \code{Float}. A scalar value in the interval between 0 and 1 corresponding ti the change point position.
#'
#' @return A \code{vector (float)} corresponding to the value of the function evaluated at \code{t}.
#'
#' @export
#'
#' @importFrom methods is
#'
#' @seealso [hurst_arctan()], [hurst_linear()].
#'
#' @examples
#' t0 <- seq(0.2, 0.8, len = 10)
#' hlogistic <- hurst_logistic(t = t0, h_left = 0.2,
#'                             h_right = 0.8, slope = 30,
#'                             change_point_position = 0.5)
#' plot(x = t0, y = hlogistic, type = "b", col = "red")
#'
#'
hurst_logistic <- function(t, h_left = 0.2, h_right = 0.8, slope = 30,
                           change_point_position = 0.5) {
  if (! (methods::is(t, "numeric") && all(t >= 0 & t <= 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  if (!(methods::is(h_left, "numeric") && methods::is(h_right, "numeric") &&
        (h_left > 0 && h_left < 1) && (h_right > 0 && h_right < 1) &&
        (length(h_left) == 1) && (length(h_right) == 1))) {
    stop("'h_left' and 'h_right' must be scalar values between 0 and 1.")
  }

  if (!(methods::is(change_point_position, "numeric") &&
        (change_point_position > 0 && change_point_position < 1) &&
        length(change_point_position) == 1)) {
    stop("'change_point_position' must be scalar value between 0 and 1.")
  }

  if (! (methods::is(slope, "numeric") && slope > 0 && length(slope) == 1)) {
    stop("'slope' must be a positive scalar value.")
  }

  u <- (t - change_point_position) / (1 - 0)
  hval <- (h_right - h_left) / (1 + exp(- slope * u)) + h_left
  return(hval)
}

#' Constant D(x,y) function
#'
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param x \code{Float (positive)}. First argument of the function.
#' @param y \code{Float (positive)}. Second argument of the function.
#'
#' @return A positive \code{Float} corresponding to the value the function evaluate at (x,y).
#'
.constant_d <- function(x, y) {
  a <- gamma(2 * x + 1) * gamma(2 * y + 1) * sin(pi * x) * sin(pi * y)
  b <- 2 * gamma(x + y + 1) * sin(pi * (x + y) / 2)
  val <- sqrt(a) / b
  return(val)
}

#' Covariance matrix of the multi-fractional Brownian Motion
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to compute the covariance function.
#' @param hurst_fun \code{function}. Hurst function. It can be \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}.
#' @param ... Hurst function additional arguments.
#'
#' @return a \code{matrix} of \code{t} x \code{t} covariance.
#'
.covariance_mfBm <- function(t = seq(0.2, 0.8, len = 10), hurst_fun = hurst_logistic, ...) {
  tmp <- expand.grid(u = t, v = t)
  u <- tmp$u
  v <- tmp$v
  hu <- hurst_fun(u, ...)
  hv <- hurst_fun(v, ...)
  hh <- hu + hv
  values <- .constant_d(hu, hv) *
    (u**hh + v**hh - abs(u - v) ** hh)
  mat <- matrix(values, ncol = length(t))
  return(mat)
}


#' Draw a multifractional Brownian motion sample path.
#'
#' This function generates a sample path of a multifractional Brownian motion (mfBm) based on the provided Hurst function and other parameters.
#'
#' @param t \code{vector (float)}. Grid of points between 0 and 1 where the sample path will be generated.
#' @param hurst_fun \code{function}. Hurst function. It can be \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}, or any custom Hurst function.
#' @param L \code{float (positive)}. Hölder constant.
#' @param shift_var \code{float (positive)}. The variance of the shift Gaussian random variable. Default is \code{shift_var = 1}, meaning a normal random variable with mean 0 and variance 1 is added.
#' @param tied \code{boolean}. If \code{TRUE}, the sample path is tied down.
#' @param ... Additional arguments for the Hurst function.
#'
#' @return A \code{data.table} containing 2 columns: \code{t} and \code{mfBm}, representing the grid points and the corresponding values of the mfBm sample path.
#'
#' @importFrom MASS mvrnorm
#' @importFrom data.table data.table
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' t0 <- seq(0.2, 0.8, len = 20)
#' dt_mfBm <- simulate_mfBm(t = t0, hurst_fun = hurst_logistic, L = 1, tied = TRUE)
#' plot(x = dt_mfBm$t, y = dt_mfBm$mfBm, type = "l", col = "red")
#'
simulate_mfBm <- function(t = seq(0.2, 0.8, len = 50), hurst_fun = hurst_logistic, L = 1, shift_var = 1, tied = TRUE, ...) {
  if (! (methods::is(t, "numeric") && all(t >= 0 & t <= 1))) {
    stop("'t' must be a numeric vector with values between 0 and 1.")
  }
  if (!methods::is(hurst_fun, "function")) {
    stop("'hurst_fun' must be a function.")
  }
  if (! (methods::is(L, "numeric") && L > 0 && length(L) == 1)) {
    stop("'L' must be a positive scalar value.")
  }
  if (! (methods::is(shift_var, "numeric") && shift_var > 0 && length(shift_var) == 1)) {
    stop("'shift_var' must be a positive scalar value.")
  }

  if (!methods::is(tied, "logical")) {
    stop("'tied' must be TRUE or FALSE.")
  }

  t <- sort(t)
  cov_mat <- .covariance_mfBm(t = t, hurst_fun = hurst_fun, ...) + shift_var
  out <- MASS::mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = L * cov_mat)
  mfBm_path <- out - tied * t * out[length(out)]
  dt <- data.table::data.table("t" = t, mfBm = mfBm_path)

  return(dt)
}



#' Draw a fractional Brownian motion sample path.
#'
#' @param t \code{vector (float)}. Grid of points between 0 and 1 where we want to generate the sample path.
#' @param hurst \code{float (positive)}. The Hurst exponent scalar value between 0 and 1.
#' @param L \code{float (positive)}. Hölder constant.
#' @param tied \code{boolean}. If \code{TRUE}, the sample path is tied-down.
#'
#' @return A \code{data.table} containing 2 column : \code{t} and \code{mfBm}, the sample path.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom data.table data.table between
#' @importFrom methods is
#'
#' @examples
#'
#' t0 <- seq(0.2, 0.8, len = 20)
#' dt_fBm <- simulate_fBm(t = t0, hurst = 0.6, L = 1, tied = TRUE)
#' plot(x = dt_fBm$t, y = dt_fBm$fBm, type = "l", col = "red")
#'
simulate_fBm <- function(t = seq(0.2, 0.8, len = 20), hurst = 0.6, L = 1, tied = TRUE) {
  if (! methods::is(t, "numeric") && all(t >= 0 & t <= 1))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(hurst, "numeric") & (hurst >= 0 & hurst <=1) & length(hurst) == 1))
    stop("'hurst' must be a positive scalar value between 0 and 1.")
  if (! (methods::is(L, "numeric") & L > 0 & length(L) == 1))
    stop("'L' must be a positive scalar value.")

  tmp <- expand.grid(u = t, v = t)
  u <- tmp$u
  v <- tmp$v
  values <- (1 / 2) *
    (u ** hurst + v ** hurst - abs(u - v) ** hurst)
  cov_mat <- matrix(values, ncol = length(t))

  out <- MASS::mvrnorm(1,
                       mu = rep(0, ncol(cov_mat)),
                       Sigma = L * cov_mat)
  fBm_path <- out - tied * t * out[length(out)]
  dt <- data.table::data.table("t" = t, "fBm" = fBm_path)
  return(dt)
}

#' Generate a random design of the place where the observation is made.
#'
#' @param N \code{integer}. Number of curves.
#' @param lambda \code{integer}. Mean of the number of observations per curve.
#' @param Mdistribution \code{function}. Distribution of the number of observation points per curve.
#' The first argument of the function must correspond to \code{N} and the second to \code{lambda}.
#' Default \code{Mdistribution = rpois}.
#' @param tdistribution \code{function}. Distribution of the observation point in the domain.
#' Currently only \code{runif} is accepted.
#' @param ... Additional argument of \code{tdistribution}.
#'
#' @return A \code{data.table} containing 3 column :
#' \itemize{
#'    \item{id_curve :}{ Index of the curve. It goes from 1 to N.}
#'    \item{Mn :}{ Number of sampled observation location.}
#'    \item{Tn :}{ Sampled observation location.}
#' }
#'
#' @import data.table
#' @importFrom methods is
#' @importFrom stats rpois runif
#'
#'
.random_design <- function(N, lambda, Mdistribution = rpois, tdistribution = runif, ...) {
  if (! (N - floor(N) == 0) & N > 1)
    stop("'N' must be an integer greater than 1.")
  if (! (lambda - floor(lambda) == 0) & lambda > 1)
    stop("'lambda' must be an integer greater than 1.")
  if (! methods::is(Mdistribution, "function"))
    stop("'Mdistribution' must be a function.")
  if (! (methods::is(tdistribution, "function") & identical(tdistribution, runif)))
    stop("'tdistribution' must be a function, and currently only 'runif' is accepted.")

  M <- Mdistribution(N, lambda)
  data.table::rbindlist(lapply(1:N, function(n, M = M, ...){
    data.table::data.table("id_curve" = n, "Mn" = M[n], "Tn" =  sort(tdistribution(M[n], ...)))
  }, M = M, ... = ...))
}

#' Functional Autoregressive process of order 1 (FAR(1)) simulation
#'
#' @param N \code{integer}. Number of curves.
#' @param lambda \code{integer}. Mean of the number of observations per curve.
#' @param tdesign \code{character}. Type of the design. It is either 'random' or 'common'.
#' @param Mdistribution \code{function}. Distribution of the number of observation points per curve.
#' The first argument of the function must correspond to \code{N} and the second to \code{lambda}.
#' Default \code{Mdistribution = rpois}.
#' @param tdistribution \code{function (or NULL)}. Observation point distribution if \code{tdesign = 'random'} and \code{NULL} otherwise.
#' @param tcommon \code{vector (float)}. Observation point vector if \code{tdesign = 'common'}.
#' If \code{tdesign = 'random'} and if we want to run some tests at a particular observation position, this can also be specified.
#' @param hurst_fun \code{function}. Hurst function. It can be \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}.
#' @param L \code{float (positive)}. Hölder constant.
#' @param far_kernel \code{function}. Kernel function of the operator of the FAR(1).
#' @param far_mean \code{function}. Mean function of the FAR(1).
#' @param int_grid \code{integer}. Length of the grid used to approximate the integral.
#' @param burnin \code{integer}. Burnin period of the FAR(1).
#' @param remove_burnin \code{boolean}. If \code{TRUE}, burnin period is removed.
#'
#' @return A \code{data.table} containing 3 column :
#' \itemize{
#'    \item{id_curve :}{ Index of the curve. It goes from 1 to N.}
#'    \item{tobs :}{ Sampled observation points, for each \code{id_curve}.}
#'    \item{ttag :}{ Tag on the observations points, for each \code{id_curve}. It is either \code{tcommon} for common design grid or \code{tcommon} pour random design.}
#'    \item{far_mean :}{ The mean of the process evaluate at \code{tobs}, for each \code{id_curve}.}
#'    \item{X :}{ The process observed at tobs, for each \code{id_curve}.}
#' }
#'
#' @importFrom data.table data.table rbindlist setnames
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#'
#'\dontrun{
#' dt_far <- simulate_far(N = 2L, lambda = 70L,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = seq(0.2, 0.8, len = 50),
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        far_kernel = function(s,t) 9/4 * exp(- (t + 2 * s) ** 2),
#'                        far_mean = function(t) 4 * sin(1.5 * pi * t),
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#'}
#'
simulate_far <- function(N = 2L, lambda = 70L,
                         tdesign = "random",
                         Mdistribution = rpois,
                         tdistribution = runif,
                         tcommon = seq(0.2, 0.8, len = 50),
                         hurst_fun = hurst_logistic,
                         L = 4,
                         far_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
                         far_mean = function(t) 4 * sin(1.5 * pi * t),
                         int_grid = 100L,
                         burnin = 100L,
                         remove_burnin = TRUE) {
  #TODO : Ajouter une description car grosse fonction
  #TODO : try sur farkernel ou alors une erreur spécifique si c'es lui qui fait péter, pareil pour far_mean
  # C'est quoi la burnin period ?
  if (! (N - floor(N) == 0) & N > 1)
    stop("'N' must be an integer greater than 1.")
  if (! methods::is(tdesign, "character")){
    stop("'tdesign' must be a character.")
  }else{
    tdesign <- match.arg(arg = tdesign, choices = c("random", "common"))
  }
  if ((tdesign == "random") & (! (lambda - floor(lambda) == 0) & lambda > 1))
    stop("'lambda' must be an integer greater than 1.")
  if (( ! (methods::is(Mdistribution, "function") & methods::is(tdistribution, "function"))) & tdesign == "random")
    stop("If tdesign = 'random', then 'Mdistribution' and 'tdistribution' must be functions.")
  if ((! (is.null(Mdistribution) & is.null(tdistribution))) & tdesign == "common")
    stop("If tdesign = 'common', then 'Mdistribution' and 'tdistribution' must be NULL.")
  if (tdesign == "common"){
    if (is.null(tcommon) | ! (any(tcommon > 0 & tcommon <= 1) & length(tcommon) > 2))
      stop("'tcommon' must be of minimum length 2 with values between 0 and 1.")
  }else{
    if (! is.null(tcommon) & ! (any(tcommon > 0 & tcommon <= 1) & length(tcommon) > 2))
      stop("If tdesign = 'random', 'tcommon' must be either NULL or of minimum length 2 with values between 0 and 1.")
  }
  if (! methods::is(hurst_fun, "function"))
    stop("'hurst_fun' must be a function.")
  if (! (methods::is(L, "numeric") & L > 0 & length(L) == 1))
    stop("'L' must be a positive scalar value.")
  if (! methods::is(far_kernel, "function"))
    stop("'far_kernel' must be bevariate function")
  if (! methods::is(far_mean, "function"))
    stop("'far_mean' must be a function")
  if (! (is.integer(int_grid) & int_grid > 50))
    stop("'int_grid' must be an integer greater than 30.")
  if (! (is.integer(burnin) & burnin > 30))
    stop("'burnin' must be an integer greater than 30.")
  if (! methods::is(remove_burnin, "logical"))
    stop("'remove_burnin' must be boolean.")
  n <- N + burnin
  grid <- (1:int_grid) / int_grid

  # If random design
  if (tdesign == "random"){

    dt_rdesign <- .random_design(N = n, lambda = lambda, Mdistribution = Mdistribution, tdistribution = tdistribution)
    M <- dt_rdesign[, unique(Mn), by = "id_curve"][, V1]

    dt_far <- data.table::rbindlist(lapply(1:n, function(i, dt_rdesign, grid, tcommon, M, hurst_fun, L){
      # Combine design + integration grid + tcommon
      tall <- c(dt_rdesign[id_curve == i, Tn], grid, tcommon)
      ttag <- c(rep("trandom", M[i]), rep("int_grid", length(grid)), rep("tcommon", length(tcommon)))
      dt <- data.table::data.table("id_curve" = i, "tall" = tall, "ttag" = ttag)
      dt <- dt[order(tall)]

      # Generate and add mfBm
      dt_eps <- simulate_mfBm(t = dt[, tall], hurst_fun = hurst_fun, L = L, tied = FALSE)
      dt[, eps := dt_eps[, mfBm]]

      # Add mean function
      dt[, far_mean := far_mean(tall)]
    }, dt_rdesign = dt_rdesign, grid = grid, tcommon = tcommon, M = M, hurst_fun = hurst_fun, L = L))
  } else {
    # Common design case
    dt_far <- data.table::rbindlist(lapply(1:n, function(i, tcommon, grid, hurst_fun, L){
      # Combine design + integration grid
      tall <- c(tcommon, grid)
      ttag <- c(rep("tcommon", length(tcommon)), rep("int_grid", length(grid)))
      dt <- data.table::data.table("id_curve" = i, "tall" = tall, "ttag" = ttag)
      dt <- dt[order(tall)]

      # Generate and add mfBm
      dt_eps <- simulate_mfBm(t = dt[, tall], hurst_fun = hurst_fun, L = L, tied = FALSE)
      dt[, eps := dt_eps[, mfBm]]

      # Add mean function
      dt[, far_mean := far_mean(tall)]
    }, tcommon = tcommon, grid = grid, hurst_fun = hurst_fun, L = L))
  }

  # Generate FAR(1)
  dt_far[id_curve == 1, X := far_mean + eps]
  for(i in 2:n){
    tall <- dt_far[id_curve == i, tall]
    Xold_centred <- dt_far[id_curve == i - 1 & ttag == "int_grid", X - far_mean]
    Xold_centred <- matrix(Xold_centred, ncol = 1)
    Enew <- dt_far[id_curve == i, eps]
    far_mean_new <- dt_far[id_curve == i, far_mean]

    tmp <- expand.grid(u = tall, v = grid)
    u <- tmp$u
    v <- tmp$v
    beta <- matrix(far_kernel(u,v), ncol = int_grid, byrow = FALSE)
    Xi <- far_mean_new + as.numeric((1/int_grid) * beta %*% Xold_centred + Enew)
    dt_far[id_curve == i, X := Xi]
  }

  # Remove the data for integral approximation
  dt_far <- dt_far[ttag != "int_grid"]
  dt_far[, eps := NULL]
  if(remove_burnin){
    dt_far <- dt_far[! id_curve %in% 1:burnin]
    dt_far[, id_curve := id_curve - burnin]
  }
  data.table::setnames(x = dt_far, old = "tall", new = "tobs")
  return(dt_far)
}

#' Functional Moving Average process of order 1 (FMA(1)) simulation
#'
#'@param N \code{integer}. Number of curves.
#' @param lambda \code{integer}. Mean of the number of observations per curve.
#' @param tdesign \code{character}. Type of the design. It is either 'random' or 'common'.
#' @param Mdistribution \code{function}. Distribution of the number of observation points per curve.
#' The first argument of the function must correspond to \code{N} and the second to \code{lambda}.
#' Default \code{Mdistribution = rpois}.
#' @param tdistribution \code{function (or NULL)}. Observation point distribution if \code{tdesign = 'random'} and \code{NULL} otherwise.
#' @param tcommon \code{vector (float)}. Observation point vector if \code{tdesign = 'common'}.
#' If \code{tdesign = 'random'} and if we want to run some tests at a particular observation position, this can also be specified.
#' @param hurst_fun \code{function}. Hurst function. It can be \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}.
#' @param L \code{float (positive)}. Hölder constant.
#' @param fma_kernel \code{function}. Kernel function of the operator of the FMA(1).
#' @param fma_mean \code{function}. Mean function of the FMA(1).
#' @param int_grid \code{integer}. Length of the grid used to approximate the integral.
#' @param burnin \code{integer}. Burnin period of the FMA(1).
#' @param remove_burnin \code{boolean}. If \code{TRUE}, burnin period is removed.
#'
#' @return A \code{data.table} containing 3 column :
#' \itemize{
#'    \item{id_curve :}{ Index of the curve. It goes from 1 to N.}
#'    \item{tobs :}{ Sampled observation points, for each \code{id_curve}.}
#'    \item{ttag :}{ Tag on the observations points, for each \code{id_curve}. It is either \code{tcommon} for common design grid or \code{tcommon} pour random design.}
#'    \item{fma_mean :}{ The mean of the process evaluate at \code{tobs}, for each \code{id_curve}.}
#'    \item{X :}{ The process observed at tobs, for each \code{id_curve}.}
#' }
#'
#' @importFrom data.table data.table rbindlist setnames
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#'
#'\dontrun{
#' dt_fma <- simulate_fma(N = 2L, lambda = 70L,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = seq(0.2, 0.8, len = 50),
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        fma_kernel = function(s,t) 9/4 * exp(- (t + 2 * s) ** 2),
#'                        fma_mean = function(t) 4 * sin(1.5 * pi * t),
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#' # plot simulated curve
#' library(ggplot2)
#'
#' ggplot(data = dt_fma[ttag == "trandom", .("id_curve" = as.factor(id_curve), tobs, X)],
#'        mapping = aes(x = tobs, y = X, group = id_curve, color = id_curve)) +
#'   geom_line() +
#'   scale_colour_grey() +
#'   theme_minimal()
#'
#'}
#'
simulate_fma <- function(N = 2L, lambda = 70L,
                         tdesign = "random",
                         Mdistribution = rpois,
                         tdistribution = runif,
                         tcommon = seq(0.2, 0.8, len = 50),
                         hurst_fun = hurst_logistic,
                         L = 4,
                         fma_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
                         fma_mean = function(t) 4 * sin(1.5 * pi * t),
                         int_grid = 100L,
                         burnin = 100L,
                         remove_burnin = TRUE) {
  if (! (N - floor(N) == 0) & N > 1)
    stop("'N' must be an integer greater than 1.")
  if (! (lambda - floor(lambda) == 0) & lambda > 1)
    stop("'lambda' must be an integer greater than 1.")
  if (! methods::is(tdesign, "character")){
    stop("'tdesign' must be a character.")
  }else{
    tdesign <- match.arg(arg = tdesign, choices = c("random", "common"))
  }
  if (( ! (methods::is(Mdistribution, "function") & methods::is(tdistribution, "function"))) & tdesign == "random")
    stop("If tdesign = 'random', then 'Mdistribution' and 'tdistribution' must be functions.")
  if ((! (is.null(Mdistribution) & is.null(tdistribution))) & tdesign == "common")
    stop("If tdesign = 'common', then 'Mdistribution' and 'tdistribution' must be NULL.")
  if (tdesign == "common"){
    if (is.null(tcommon) | ! (any(tcommon > 0 & tcommon <= 1) & length(tcommon) > 2))
      stop("'tcommon' must be of minimum length 2 with values between 0 and 1.")
  }else{
    if (! is.null(tcommon) & ! (any(tcommon > 0 & tcommon <= 1) & length(tcommon) > 2))
      stop("If tdesign = 'random', 'tcommon' must be either NULL or of minimum length 2 with values between 0 and 1.")
  }
  if (! methods::is(hurst_fun, "function"))
    stop("'hurst_fun' must be a function.")
  if (! (methods::is(L, "numeric") & L > 0 & length(L) == 1))
    stop("'L' must be a positive scalar value.")
  if (! methods::is(fma_kernel, "function"))
    stop("'fma_kernel' must be bevariate function")
  if (! methods::is(fma_mean, "function"))
    stop("'fma_mean' must be a function")
  if (! (is.integer(int_grid) & int_grid > 50))
    stop("'int_grid' must be an integer greater than 30.")
  if (! (is.integer(burnin) & burnin > 30))
    stop("'burnin' must be an integer greater than 30.")
  if (! methods::is(remove_burnin, "logical"))
    stop("'remove_burnin' must be boolean.")
  n <- N + burnin
  grid <- (1:int_grid) / int_grid

  # If random design
  if (tdesign == "random"){
    dt_rdesign <- .random_design(N = n, lambda = lambda, Mdistribution = Mdistribution, tdistribution = tdistribution)
    M <- dt_rdesign[, unique(Mn), by = "id_curve"][, V1]

    dt_fma <- data.table::rbindlist(lapply(1:n, function(i, dt_rdesign, grid, tcommon, M, hurst_fun, L){
      # Combine design + integration grid + tcommon
      tall <- c(dt_rdesign[id_curve == i, Tn], grid, tcommon)
      ttag <- c(rep("trandom", M[i]), rep("int_grid", length(grid)), rep("tcommon", length(tcommon)))
      dt <- data.table::data.table("id_curve" = i, "tall" = tall, "ttag" = ttag)
      dt <- dt[order(tall)]

      # Generate and add mfBm
      dt_eps <- simulate_mfBm(t = dt[, tall], hurst_fun = hurst_fun, L = L, tied = FALSE)
      dt[, eps := dt_eps[, mfBm]]

      # Add mean function
      dt[, fma_mean := fma_mean(tall)]
    }, dt_rdesign = dt_rdesign, grid = grid, tcommon = tcommon, M = M, hurst_fun = hurst_fun, L = L))
  } else {
    # Common design case
    dt_fma <- data.table::rbindlist(lapply(1:n, function(i, tcommon, grid, hurst_fun, L){
      # Combine design + integration grid
      tall <- c(tcommon, grid)
      ttag <- c(rep("tcommon", length(tcommon)), rep("int_grid", length(grid)))
      dt <- data.table::data.table("id_curve" = i, "tall" = tall, "ttag" = ttag)
      dt <- dt[order(tall)]

      # Generate and add mfBm
      dt_eps <- simulate_mfBm(t = dt[, tall], hurst_fun = hurst_fun, L = L, tied = FALSE)
      dt[, eps := dt_eps[, mfBm]]

      # Add mean function
      dt[, fma_mean := fma_mean(tall)]
    }, tcommon = tcommon, grid = grid, hurst_fun = hurst_fun, L = L))
  }

  # Generate FAR(1)
  dt_fma[id_curve == 1, X := fma_mean + eps]
  for(i in 2:n){
    tall <- dt_fma[id_curve == i, tall]
    Eold <- dt_fma[id_curve == i - 1 & ttag == "int_grid", eps]
    Enew <- dt_fma[id_curve == i, eps]
    fma_mean_new <- dt_fma[id_curve == i, fma_mean]

    tmp <- expand.grid(u = tall, v = grid)
    u <- tmp$u
    v <- tmp$v
    beta <- matrix(fma_kernel(u,v), ncol = int_grid, byrow = FALSE)
    Xi <- fma_mean_new + Enew + as.numeric((1 / int_grid) * beta %*% matrix(Eold, ncol = 1))
    dt_fma[id_curve == i, X := Xi]
  }

  # Remove the data for integral approximation
  dt_fma <- dt_fma[ttag != "int_grid"]
  dt_fma[, eps := NULL]
  if(remove_burnin){
    dt_fma <- dt_fma[! id_curve %in% 1:burnin]
    dt_fma[, id_curve := id_curve - burnin]
  }
  data.table::setnames(x = dt_fma, old = "tall", new = "tobs")
  return(dt_fma)
}

