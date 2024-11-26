% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_estimate_mean.R
\name{estimate_mean_risk}
\alias{estimate_mean_risk}
\title{Estimate the risk of the mean function}
\usage{
estimate_mean_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw_grid = seq(0.005, 0.15, len = 45),
  Ht = NULL,
  Lt = NULL,
  Delta = NULL,
  h = NULL,
  smooth_ker = epanechnikov
)
}
\arguments{
\item{data}{\code{data.table (or data.frame)} or \code{list} of \code{data.table (or data.frame)} or \code{list} of \code{list}.
\itemize{
   \item{If \code{data.table}}{
       It must contain the raw binding of the curve observations with at least 3 columns.
       \itemize{
         \item{\code{idcol} :}{ The name of the column that contains the index of the curve in the sample.
                             Each index of a curve is repeated as many times as it has observation points.}
         \item{\code{tcol} :}{ The name of the column that contains the observation points associated to each curve index.}
         \item{\code{ycol} :}{ The name of the column that contains the observed value of a curve at each point of observation and for each index of the curve.}
       }
   }
   \item{If \code{list} of \code{data.table}}{
        In this case, each element of the given \code{list} corresponds to the observation scheme of a curve, which is given as \code{data.table} or \code{data.frame}.
        The data.table contains at least 2 columns.
        \itemize{
         \item{\code{tcol} :}{ The name of the column that contains the observation points associated to the curve.}
         \item{\code{ycol} :}{ The name of the column that contains the observed value of the curve.}
       }
   }
   \item{If \code{list} of \code{list}}{
     In the latter case, the \code{data} is a list \code{list} where each element is the observation scheme of a curve given as a \code{list} of 2 vectors.
     \itemize{
         \item{\code{tcol} :}{ The name of the vector that contains the observation points associated the curve.}
         \item{\code{ycol} :}{ The name of the vector that contains the observed value of the curve.}
       }
   }
}}

\item{idcol}{\code{character}. If \code{data} is given as \code{data.table} or \code{data.frame},
it is the name of the column that contains the index of the curve in the sample.
Each index of a curve is repeated as many times as it has observation points.
Opposite, if f \code{data} is given as \code{list} of \code{data.table (or data.frame)} of \code{list} of \code{list}, \code{idcol = NULL.}}

\item{tcol}{\code{character}. The name of the column (or vector) that contains the observation points associated to the curves.}

\item{ycol}{\code{character}. The name of the column that contains the observed value of the curves.}

\item{t}{\code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.}

\item{bw_grid}{\code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each \code{t}.
It can be \code{NULL} and that way it will be defined as an exponential grid of \eqn{N\times\lambda}.}

\item{Ht}{\code{vector (numeric)}. The estimates of the local exponent for each \code{t}.
Default \code{Ht = NULL} and thus it will be estimated.}

\item{Lt}{\code{vector (numeric)}. The estimates of the Hölder constant for each \code{t}.
It corresponds to \eqn{L_t^2}. Default \code{Lt = NULL} and thus it will be estimated.}

\item{Delta}{\code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
Default \code{Delta = NULL} and thus it will be estimated from the data.}

\item{h}{\code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
and each element of the vector must correspond to a curve given in the same order as in \code{data}.}

\item{smooth_ker}{\code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.}
}
\value{
A \code{data.table} containing the following columns.
         \itemize{
           \item{t :}{ The points at which the risk function is estimated.}
           \item{h :}{ The candidate bandwidth.}
           \item{Ht :}{ The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}}
           \item{Lt :}{ The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
           \item{locreg_bw :}{ The bandwidth used to estimate the local regularity parameters.}
           \item{bias_term :}{ The bias term of the risk function.}
           \item{varriance_term :}{ The variance term of the risk function.}
           \item{dependence_term :}{ The dependence term of the risk function.}
           \item{mean_risk :}{ The estimates of the risk function of the mean.}
        }
}
\description{
The risk of the mean function is the function \eqn{R_\mu(t;h)} in
Section 4.1 of \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
}
\examples{
\dontrun{
# Generate a sample path of FTS
dt_far <- simulate_far(N = 50, lambda = 70,
                       tdesign = "random",
                       Mdistribution = rpois,
                       tdistribution = runif,
                       tcommon = NULL,
                       hurst_fun = hurst_logistic,
                       L = 4,
                       far_kernel = get_real_data_far_kenel,
                       far_mean = get_real_data_mean,
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)

# Add noise
dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]

# Estimate risk function
dt_mean_risk <- estimate_mean_risk(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
  Delta = NULL, h = NULL, smooth_ker = epanechnikov)

# Plot mean risk
dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
      main = "t = 0.25", xlab = "h", ylab = "risk function"),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
      main = "t = 0.5", xlab = "h", ylab = "risk function"),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
      main = "t = 0.75", xlab = "h", ylab = "risk function")
  ),
  nrow = 3
)

}


}
\references{
\insertAllCited{}
}
\seealso{
[estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_autocov()].
}
