#' Generate Error-Prone Data
#' 
#' @description
#' Generates a design matrix as well as its error-prone version with the specified number of replicates, along with the response varaible, based on specified parameters.
#' @param n Number of observations.
#' @param p Number of predictors.
#' @param k Number of replicates.
#' @param mu Population intercept (scalar).
#' @param mu_x p-dimensional mean vector for error-prone predictors.
#' @param beta p-dimensional vector of true coefficients (used to generate the response variable).
#' @param Sig_x Covariance matrix for true predictors X.
#' @param Sig_u Covariance matrix for measurement error U.
#' @param sig2e Variance of response (or, equivalently, of residual noise) in linear regression (default 0.1).
#' @param family Distribution of response conditional on true predictors X ("gaussian" for linear regression, or "binomial" for logistic regression).
#' @return A list with true predictor matrix X, list of k error-prone replicates W, response variable y, and n, p, and k.
#' @export

data_gen <- function(n, p, k, mu=0, mu_x=NULL, beta, Sig_x, Sig_u, sig2e=0.1, family="gaussian"){
  if (is.null(mu_x)){
    mu_x <- rep(0, p)
  }
  X <- MASS::mvrnorm(n, mu=mu_x, Sigma=Sig_x)
  
  if (family == "binomial"){
    z <- mu + X %*% beta
    pr <- 1/(1+exp(-z))
    y <- rbinom(n, 1, pr)
  } else if (family == "gaussian"){
    e <- rnorm(n, 0, sqrt(sig2e))
    y <- mu + X %*% beta + e
  } else {
    stop(paste0("cannot generate data for family ", family))
  }
  
  W <- list()
  if (k == 1){
    W[[1]] <- X + MASS::mvrnorm(n, mu=rep(0,p), Sigma=Sig_u)
  } else {
    for (r in 1:k){
      W[[r]] <- X + MASS::mvrnorm(n, mu=rep(0,p), Sigma=Sig_u)
    }
  }
  
  return (
    list(X=X, W=W, y=y, n=n, p=p, k=k)
  )
}
