#' Biweight kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [triweight()], [tricube()], [epanechnikov()], [triangular()], and [uniform()].
#'
biweight <- function(u){
  ifelse(abs(u) <= 1, (15 / 16) * (1 - u ** 2) ** 2, 0)
}

#' Triweight kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [tricube()], [epanechnikov()], [triangular()], and [uniform()].
#'
triweight <- function(u){
  ifelse(abs(u) <= 1, (35 / 32) * (1 - u ** 2) ** 3, 0)
}

#' Tricube kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [epanechnikov()], [triangular()], and [uniform()].
#'
tricube <- function(u){
  ifelse(abs(u) <= 1, (70 / 81) * (1 - abs(u) ** 3) ** 3, 0)
}

#' Epanechnikov kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [tricube()], [triangular()], and [uniform()].
#'
epanechnikov <- function(u){
  ifelse(abs(u) <= 1, (3 / 4) * (1 - u ** 2), 0)
}

#' Triangular kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [tricube()], [epanechnikov()], and [uniform()].
#'
triangular <- function(u){
  ifelse(abs(u) <= 1, (1 - abs(u)), 0)
}

#' Uniform kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [tricube()], [epanechnikov()], and [triangular()].
#'
uniform <- function(u){
  ifelse(abs(u) <= 1, 1/2, 0)
}

#' Nadaraya-Watson estimator
#'
#' @param y \code{vector (numeric)}. A numeric vector containing the observed values of the independent variable corresponding to the observation points \code{t}.
#' @param t \code{vector (numeric)}. A numeric vector containing the observed values of the dependent variable.
#' @param tnew \code{vector (numeric)}. New \code{t} values at which we want to estimate the regression function.
#' @param h \code{numeric (positive)}. The bandwidth parameter such that \code{h > (2 * length(x))}.
#' Default \code{h = NULL} and such it will be computed automatically.
#' @param smooth_ker \code{function}. The kernel function of the estimator.
#'
#' @return A \code{data.table} containing
#' \itemize{
#'    \item{h :}{ The bandwidth used to estimate the regression function.}
#'    \item{inKernelSupp :}{ For each \code{t_i} in \code{tnew}, it is the number of  \code{t} between \code{t_i - h } and \code{t_i + h }.}
#'    \item{tnew :}{ The vector \code{new}.}
#'    \item{yhat :}{ The regression function's vector of estimates at \code{tnew}.}
#' }
#'
#' @importFrom methods is
#' @importFrom data.table data.table
#' @importFrom matrixStats rowSums2
#'
#' @export
#'
#' @seealso [estimate_nw_bw()], [epanechnikov()], [biweight()], [triweight()], [tricube()], [uniform()], etc.
#'
#' @examples
#' \dontrun{
#' # The model
#' ## Let
#' m <- function(t) 4 * sin(1.5 * pi * t)
#'
#' ## Observation points
#' t <- runif(n = 200, min = 0, max = 1)
#' t <- sort(t)
#'
#' ## Measure error
#' e <- rnorm(n = 200, mean = 0, sd = 0.2)
#'
#' ## Regression model
#' y <- m(t) + e
#'
#' plot(x = t, y = y, main = "Observed points and true regression function")
#' lines(x = t, y = m(t), type = "l", col = "red")
#'
#' ## Estimate the best bandwidth
#' bw_grid <- seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100)
#' hbest <- estimate_nw_bw(y = y, t = t,
#'                         bw_grid = bw_grid,
#'                         smooth_ker = epanechnikov)
#'
#' ## Estimate the regression function
#' dt_nw <- estimate_nw(y = y, t = t,
#'                      tnew = seq(0.01, 0.99, len = 100),
#'                      h = hbest, smooth_ker = epanechnikov)
#'
#' plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
#'      main = "Estimated and true regression function.")
#' lines(x = dt_nw[, tnew], y = m(dt_nw[, yhat]), type = "l", col = "red")
#' legend(x = 0.64, y = 4.1, fill = c("blue", "red"),legend = c("Estimated m", "True m"))
#'
#' }
#'
estimate_nw <- function(y, t, tnew, h = NULL, smooth_ker = epanechnikov){
  if (! is.numeric(y) | ! is.numeric(t) | ! is.numeric(tnew) |! is.numeric(h))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (! is.null(h) & length(h) > 1)
    stop("The bandwidth 'h' must be either NULL or saclar.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  m <- length(tnew)
  n <- length(y)
  A <- matrix(0, nrow = m, ncol = n)
  if (is.null(h)) {
    hcv <- estimate_nw_bw(y = y, t = t,
                          bw_grid = seq(2 / n, n ** (-1/3), len = 100),
                          smooth_ker = smooth_ker)
  }
  h <- ifelse(is.null(h), hcv, h)
  A <- outer(tnew, t, function(u, v) smooth_ker((u - v) / h))
  A <- A / matrixStats::rowSums2(A)
  yhat <- A %*% y
  yhat <- c(yhat)

  # ## Get the number of points used to estimate y for each tnew
  # inKernelSupp <- outer(tnew, t, function(u, v) abs(u - v) <= h)
  # inKernelSupp <- matrixStats::rowSums2(inKernelSupp)
  # inKernelSupp <- c(inKernelSupp)
  # inKernelSupp <- unlist(lapply(tnew, function(tnewi, t, h){
  #   sum(abs(tnewi - t) <= h)
  # }, t = t, h = h))

  dt <- data.table::data.table("h" = h,
                               # "inKernelSupp" = inKernelSupp,
                               "tnew" = tnew, "yhat" = yhat)
  return(dt)
}

#' Nadaraya-Watson Bandwidth Selection using cross validation.
#'
#' @param y \code{vector (numeric)}. A numeric vector containing the observed values of the independent variable corresponding to the observation points \code{t}.
#' @param t \code{vector (numeric)}. A numeric vector containing the observed values of the dependent variable.
#' @param bw_grid \code{vector (numeric)}. A grid of bandwidth to test. Default \code{bw_grid = NULL}, so it will be set as an exponential grid of \code{length(t)}.
#' @param smooth_ker \code{function}. The kernel function of the estimator.
#'
#' @return A \code{numeric} value corresponding to the best bandwidth.
#' @export
#' @seealso [estimate_nw()]
#'
#' @examples
#' \dontrun{
#' # The model
#' ## Let
#' m <- function(t) 4 * sin(1.5 * pi * t)
#'
#' ## Observation points
#' t <- runif(n = 200, min = 0, max = 1)
#' t <- sort(t)
#'
#' ## Measure error
#' e <- rnorm(n = 200, mean = 0, sd = 0.2)
#'
#' ## Regression model
#' y <- m(t) + e
#'
#' plot(x = t, y = y, main = "Observed points and true regression function")
#' lines(x = t, y = m(t), type = "l", col = "red")
#'
#' ## Estimate the best bandwidth
#' bw_grid <- seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100)
#' hbest <- estimate_nw_bw(y = y, t = t,
#'                         bw_grid = bw_grid,
#'                         smooth_ker = epanechnikov)
#'
#' ## Estimate the regression function
#' dt_nw <- estimate_nw(y = y, t = t,
#'                      tnew = seq(0.01, 0.99, len = 100),
#'                      h = hbest, smooth_ker = epanechnikov)
#'
#' plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
#'      main = "Estimated and true regression function.")
#' lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), type = "l", col = "red")
#' legend(x = 0.64, y = 4.1, fill = c("blue", "red"),legend = c("Estimated m", "True m"))
#'
#' }
#'
estimate_nw_bw <- function(y, t, bw_grid = NULL,
                           smooth_ker = epanechnikov) {
  if (! is.numeric(y) | ! is.numeric(t) |! is.numeric(bw_grid))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (! is.numeric(bw_grid) & ! is.numeric(bw_grid))
    stop("If the bandwidth grid 'bw_grid' is not NULL, so it must be a scalar or vector of numeric.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  if (is.null(bw_grid)) {
    K <- 100
    b0 <- 1 / length(t)
    bK <- length(t) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(b0, bK, a, K) ; gc()
  }

  cv_error <- sapply(bw_grid, function(hi, y, t, K){
    yhat <- estimate_nw(y = y, t = t, h = hi, tnew = t, smooth_ker = K)$yhat
    wmat <- outer(X = t, Y = t, function(u, v) K((u-v) / hi))
    metric <- (y - yhat) / (1 - K(0) / rowSums(wmat))

    # If there is only one value in the kernel support, it return NaN.
    error_hi <- mean(metric[! is.nan(metric)] ** 2)
  }, y = y, t = t, K = smooth_ker)

  # If cv_error is NaN, do take it into account
  if (any(is.nan(cv_error))) {
    bw_grid <- bw_grid[-which(is.nan(cv_error))]
    cv_error <- cv_error[! is.nan(cv_error)]
    hcv <- bw_grid[which.min(cv_error)]
  } else {
    hcv <- bw_grid[which.min(cv_error)]
  }
  return(hcv)
}

#' Estimate Nadayara-Watson optimal bandwidth on all or a subset of curves
#'
#' @inheritParams .format_data
#' @param bw_grid \code{vector (numeric)}. The cross-validation bandwidth grid.
#' Default \code{bw_grid = NULL} and so it will be set as an exponential grid using the average of the number of observation points per curve.
#' @param nsubset \code{integer (positive integer)}. The number of subset curves to be randomly and uniformly selected.
#' Default \code{nsubset = NULL} and thus an optimal bandwidth is calculated for each curve.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator.
#' Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{id_curve :}{The index of the curve.}
#'            \item{optbw :}{ The optimal bandwidth obtained by Cross-Validation.}
#'         }
#' @export
#'
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist
#'
get_nw_optimal_bw <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                              bw_grid = NULL, nsubset = NULL, smooth_ker = epanechnikov){
  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Control parameters
  if ((!is.null(bw_grid)) & (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1)))
    stop("If'bw_grid' is not NULL, then it must be a vector of positive values between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (!is.null(nsubset))
    if (any(nsubset < 0) | (length(nsubset) > 1) | any(nsubset - floor(nsubset) > 0) | any(N <= nsubset))
      stop("If 'nsubset' is not NULL, then if must be a positive integer lower than the number of curves.")

  # Define the set of curve
  if (! is.null(nsubset)) {
    sample_curves <- sample(x = 1:N, size = nsubset)
  } else {
    sample_curves <- 1:N
  }

  # Define the grid
  if (is.null(bw_grid)) {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    K <- 100
    b0 <- 2 / lambdahat
    bK <- lambdahat ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))

    rm(K, b0, bK, a) ; gc()
  } else {
    NA
  }

  dt_optbw <- data.table::rbindlist(lapply(sample_curves, function(i, bw_grid, data, smooth_ker){
    K <- 100
    b0 <- 2 / lambdahat
    bK <- lambdahat ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    hgrid <- b0 * a ** (seq_len(K))
    hbest <- estimate_nw_bw(y = data[id_curve == i, X],
                            t = data[id_curve ==i, tobs],
                            bw_grid = bw_grid,
                            smooth_ker = smooth_ker)
    dt_res <- data.table::data.table("id_curve" = i, "optbw" = hbest)
    return(dt_res)
  }, bw_grid = bw_grid, data = data, smooth_ker = smooth_ker))

  return(dt_optbw)
}
