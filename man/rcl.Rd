% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcl.R
\name{rcl}
\alias{rcl}
\title{Regression-Calibrated Lasso}
\usage{
rcl(
  me_list,
  SCov_x,
  SCov_u,
  SCov_z = NULL,
  SCov_xz = NULL,
  Lambda = NULL,
  fun = NULL,
  ...
)
}
\arguments{
\item{me_list}{List containing, at the very least, an n-dimensional vector indicating the number of replicates per subject, k, the response y, a list of n matrices of error-prone observations W -- one for each subject and with dimension k\link{i} x p, and a matrix of error-free predictors Z if one desires to include it in the model.}

\item{SCov_x}{p x p dimensional matrix for the estimated covariance matrix for true predictors X.}

\item{SCov_u}{p x p dimensional matrix for the estimated covariance matrix for the measurement errors U.}

\item{SCov_z}{q x q dimensional matrix for the estimated covariance matrix for the error-free predictors Z (optional).}

\item{SCov_xz}{p x q dimensional matrix for the estimated cross-covariance matrix between X and Z (optional).}

\item{Lambda}{(p + q) x p dimensional matrix for the attenuation matrix. If provided, it is used instead of any of the other provided covariance matrix estimates. Only implemented for equal number of replicates per subject. Default NULL.}

\item{fun}{function used to estimate precision matrices. Must take as the first input the matrix to invert, followed by additional arguments.}

\item{...}{Additional inputs to "fun" and to the glmnet function.}
}
\value{
An object of class "glmnet".
}
\description{
Lasso for linear or logistic regression with regression calibration correction for measurement error.
}
