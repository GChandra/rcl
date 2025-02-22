#' Generate Error-Prone Data
#' 
#' @description
#' Generates a design matrix as well as its error-prone version with the specified number of replicates, along with the response variable, based on specified parameters. If any of the parameters required to generate error-free predictors Z is missing, Z will not be generated.
#' @param n Number of observations.
#' @param p Number of error-prone predictors.
#' @param q Number of error-free predictors (can be omitted if no Z desired).
#' @param k Number of replicates.
#' @param mu Population intercept (scalar).
#' @param mu_x p-dimensional mean vector for error-prone predictors.
#' @param mu_z q-dimensional mean vector for error-free predictors.
#' @param beta p-dimensional vector of true coefficients for error-prone predictors (used to generate the response variable).
#' @param gamma q-dimensional vector of true coefficients for error-free predictors.
#' @param Sig_x Covariance matrix for true error-prone predictors X.
#' @param Sig_u Covariance matrix for measurement error U.
#' @param Sig_z Covariance matrix for error-free predictors Z (can be omitted if no Z desired).
#' @param Sig_xz (p x q)-dimensional covariance matrix for the covariance between X and Z (can be omitted if no Z desired).
#' @param sig2e Variance of response (or, equivalently, of residual noise) in linear regression (default 0.1).
#' @param family Distribution of response conditional on true predictors X ("gaussian" for linear regression, or "binomial" for logistic regression).
#' @return A list with true predictor matrix X, list of k error-prone replicates W, response variable y, and n, p, and k.
#' @export

data_gen <- function(n, p, q=NULL, k,
                     mu=0, mu_x=NULL, mu_z=NULL,
                     beta, gamma=NULL,
                     Sig_x, Sig_u, Sig_z=NULL, Sig_xz=NULL,
                     sig2e=0.1, family="gaussian"){
  if (is.null(mu_x)){
    mu_x <- rep(0, p)
  }
  
  if ( any(sapply(list(q, gamma, Sig_z, Sig_xz), is.null)) ){
    no_Z <- TRUE
    X <- MASS::mvrnorm(n, mu=mu_x, Sigma=Sig_x)
    q <- 0
    gamma <- matrix(0, nrow=1, ncol=1)
    Z <- rep(0,n)
    Sig_z <- NULL
    Sig_xz <- NULL
  } else {
    no_Z <- FALSE
    if (is.null(mu_z)){
      mu_x <- rep(0, q)
    }
    pred <- MASS::mvrnorm(n, mu=c(mu_x, mu_z),
                          Sigma=cbind(
                            rbind(Sig_x, t(Sig_xz)),
                            rbind(Sig_xz, Sig_z)
                          ))
    X <- pred[,1:p]
    Z <- pred[,(p+1):(p+q)]
  }
  
  if (family == "binomial"){
    z <- mu + X %*% beta + Z %*% gamma
    pr <- 1/(1+exp(-z))
    y <- rbinom(n, 1, pr)
  } else if (family == "gaussian"){
    e <- rnorm(n, 0, sqrt(sig2e))
    y <- mu + X %*% beta + Z %*% gamma + e
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
  
  if (no_Z){
    Z <- NULL
  }
  
  return (
    list(X=X, W=W, Z=Z, y=y, n=n, p=p, q=q, k=k)
  )
}
