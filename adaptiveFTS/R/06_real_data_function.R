#' Mean function learned from the voltage curves of the electricity
#'
#' For more details see the vignette:
#' \code{vignette("hybrid-simulation-setup", package = "adaptiveFTS")}
#'
#' @param t \code{vector (numeric)}. Points at which we want to return the mean function.
#' It can be a scalar.
#'
#' @return A \code{data.table} containing 2 columns.
#'          \itemize{
#'            \item{t :}{ The vector or scalar \code{t}.}
#'            \item{mean :}{ The values of the mean function evaluated at \code{t}.}
#'         }
#' @export
#'
#' @importFrom data.table data.table
#'
#' @examples
#'
#' t0 <- seq(0.1, 0.9, len = 10)
#' m <- get_real_data_mean(t = t0)
#' plot(x = t0, y = m, type = "b", col = "red",
#'      xlab = "t", ylab = "mean", main = "Mean function")
#'
get_real_data_mean <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:5, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:5, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    240.851203112216, 0.509378236915314, 0.0666785737279956, 0.402943145860831,
    0.161933581079031, 0.112863126651063, 0.420525704902966, 1.00346816098248,
    -0.242895339672357, -0.259141006436404, 0.00114630490804474
  )

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(cost_mat, sint_mat, eta, basis_coef)
  gc()

  return(muhat[, 1])
}

#' FAR kernel learned from the voltage curves of the electricity
#'
#' For more details see the vignette:
#' \code{vignette("hybrid-simulation-setup", package = "adaptiveFTS")}
#'
#' @param s \code{numeric (positive)}. A vector or scalar value(s) between 0 and 1.
#' @param t \code{numeric (positive)}. A vector or scalar value(s) between 0 and 1.
#' @param operator_norm \code{numeric (positive)}. A scalar corresponding to the norm of the integral operator associated with this kernel function.
#'
#' @return A vector (or scalar) of \code{numeric} values corresponding to the value of the kernel function evaluated at (\code{s}, \code{t}).
#' @export
#'
#' @examples
#'
#' # get the value of the kernel at (s,t) = (0.2, 0.3)
#' kerval <- get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
#' kerval
#'
get_real_data_far_kenel <- function(s = 0.2, t = 0.3, operator_norm = 0.5){
  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    0.887265486496153, -0.158284777828367, -0.433123270896265, -0.383368407909871,
    0.145655492369033, -0.00932791858596785, 0.25405721049976, 0.0360507006945973,
    0.0389539855934984, 0, -0.0133553863644848, 0.0177582032888235, 0.189761421268642,
    0.0195864450427664, 0.0887495150023169, 0, 0.0347257788913602, 0, 0.298938773778208,
    0.360062724244617, 0.00694075505838772, 0.0383993219719295, 0.0889742879270508,
    0.108124616829882, 0.597015339786177
  )
  # Transform to (K, L) matrix
  basis_coef_mat <- t(matrix(data = basis_coef, ncol = 5))

  ker_values <- mapply(function(s,t, coef_mat){
    # \eta(s)
    etas <- c(1, sqrt(2) * cos(2 * pi * 1:2 * s), sqrt(2) * sin(2 * pi * 1:2 * s))

    # \theta(t)
    thetat <- c(1, sqrt(2) * cos(2 * pi * 1:2 * t), sqrt(2) * sin(2 * pi * 1:2 * t))

    # Basis function
    ker_val <- matrix(etas, nrow = 1) %*% coef_mat %*% matrix(thetat, ncol = 1)
    return(c(ker_val))
  }, s = s, t = t, MoreArgs = list(coef_mat = basis_coef_mat))

  # Normalize values using operator norm
  op_norm <- 0.9128311
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef)
  gc()

  return(ker_values)
}

