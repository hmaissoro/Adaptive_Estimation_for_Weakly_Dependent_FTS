% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_estimate_regularity.R
\name{estimate_locreg}
\alias{estimate_locreg}
\title{Local Regularity Parameters Estimation}
\usage{
estimate_locreg(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = 1/2,
  Delta = NULL,
  h = NULL,
  smooth_ker = epanechnikov,
  center = TRUE
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

\item{t}{\code{vector (numeric)}. Observation points at which we want to estimate the local regularity parameters of the underlying process.}

\item{Delta}{\code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
Default \code{Delta = NULL} and thus it will be estimated from the data.}

\item{h}{\code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
and each element of the vector must correspond to a curve given in the same order as in \code{data}.}

\item{smooth_ker}{\code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.}

\item{center}{\code{logical}. If \code{TRUE}, the curves are centered.}
}
\value{
A \code{data.table} containing the following columns.
         \itemize{
           \item{t :}{ The points around which the local regularity parameters are estimated.}
           \item{locreg_bw :}{ The presmoothing bandwidth.}
           \item{Delta :}{ The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.}
           \item{Nused :}{ The number of curves that give non-degenerate estimates around \code{t}.}
           \item{Ht :}{ The local exponent estimates for each \code{t}. It corresponds to \eqn{H_t}}
           \item{Lt :}{ The Hölder constant estimates \code{t}. It corresponds to \eqn{L_t^2}.}
        }
}
\description{
Local Regularity Parameters Estimation
}
\examples{
\dontrun{
  # Generate a sample of FAR(1)
Hfun <- function(t) {
  hurst_logistic(t = t, h_left = 0.4, h_right = 0.8, slope = 5)
}

## Hölder constant
L <- 4

dt_far <- simulate_far(N = 200L, lambda = 100L,
                       tdesign = "random",
                       Mdistribution = rpois,
                       tdistribution = runif,
                       tcommon = NULL,
                       hurst_fun = Hfun,
                       L = L,
                       far_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
                       far_mean = function(t) 4 * sin(1.5 * pi * t),
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)

# Estimate local regularity at
t0 <- seq(0.2, 0.8, len = 8)

## If data is a data.table or a data. frame
dt_locreg <- estimate_locreg(data = dt_far,
                             idcol = "id_curve",
                             tcol = "tobs",
                             ycol = "X",
                             t = t0,
                             Delta = NULL,
                             h = NULL,
                             smooth_ker = epanechnikov)
DT::datatable(dt_locreg)

## If data is a list of data.table (or data. frame)
list_dt_far <- lapply(unique(dt_far[, id_curve]), function(idx){
  dt_far[id_curve == idx, list(tobs, X)]
})

dt_locreg_2 <- estimate_locreg(data = list_dt_far,
                               idcol = NULL,
                               tcol = "tobs",
                               ycol = "X",
                               t = t0,
                               Delta = NULL,
                               h = NULL,
                               smooth_ker = epanechnikov)
DT::datatable(dt_locreg_2)

## If data is a list of list
list_list_far <- lapply(unique(dt_far[, id_curve]), function(idx){
  list("Obs_point" = dt_far[id_curve == idx, tobs],
       "Xobs" = dt_far[id_curve == idx, X])
})

dt_locreg_3 <- estimate_locreg(data = list_list_far,
                               idcol = NULL,
                               tcol = "Obs_point",
                               ycol = "Xobs",
                               t = t0,
                               Delta = NULL,
                               h = NULL,
                               smooth_ker = epanechnikov)
DT::datatable(dt_locreg_2)

}



}
\seealso{
[estimate_nw()], [estimate_nw_bw()], [simulate_far()], etc.
}
