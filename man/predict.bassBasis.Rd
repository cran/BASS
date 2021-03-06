% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bassPCA.R
\name{predict.bassBasis}
\alias{predict.bassBasis}
\title{BASS Prediction}
\usage{
\method{predict}{bassBasis}(
  object,
  newdata,
  mcmc.use = NULL,
  trunc.error = FALSE,
  nugget = T,
  n.cores = 1,
  parType = "fork",
  ...
)
}
\arguments{
\item{object}{a fitted model, output from the \code{bass} function.}

\item{newdata}{a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bassBasis} or \code{bassPCA}  functions.}

\item{mcmc.use}{a vector indexing which MCMC iterations to use for prediction.}

\item{trunc.error}{logical, use basis truncation error when predicting?}

\item{nugget}{logical, use individual \code{bass} nugget variances when predicting?}

\item{n.cores}{number of cores, though 1 is often the fastest.}

\item{parType}{either "fork" or "socket".  Forking is typically faster, but not compatible with Windows. If \code{n.cores==1}, \code{parType} is ignored.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An array with first dimension corresponding to MCMC iteration, second dimension corresponding to the rows of \code{newdata}, and third dimension corresponding to the multivariate/functional response.
}
\description{
Predict function for BASS.  Outputs the posterior predictive samples based on the specified MCMC iterations.
}
\details{
Prediction combined across \code{bass} models.
}
\examples{
# See examples in bass documentation.

}
\seealso{
\link{bassPCA} and \link{bassBasis} for model fitting and \link{sobolBasis} for sensitivity analysis.
}
