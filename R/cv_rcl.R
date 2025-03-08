#' Cross-Validated Regression-Calibrated Lasso
#' 
#' @description
#' Cross-validation for RC Lasso for linear or logistic regression.
#' @param me_list List containing, at the very least, an n-dimensional vector indicating the number of replicates per subject, k, the response y, a list of n matrices of error-prone observations W -- one for each subject and with dimension k[i] x p, and a matrix of error-free predictors Z if one desires to include it in the model.
#' @param SCov_x p x p dimensional matrix for the estimated covariance matrix for true predictors X.
#' @param SCov_u p x p dimensional matrix for the estimated covariance matrix for the measurement errors U.
#' @param SCov_z q x q dimensional matrix for the estimated covariance matrix for the error-free predictors Z (optional).
#' @param SCov_xz p x q dimensional matrix for the estimated cross-covariance matrix between X and Z (optional).
#' @param fun function used to estimate precision matrices. Must take as the first input the matrix to invert, followed by additional arguments.
#' @param ... Additional inputs to "fun" and to the glmnet function.
#' @return An object of class "cv.glmnet".
#' @examples
#' if(!require("glasso")){
#'   stop("this example requires package 'glasso'")
#' }
#' 
#' set.seed(1424)
#' 
#' inv_fun <- function(M, rho, ...){
#'   return( glasso::glasso(M, rho=rho)$wi )
#' }
#' 
#' n <- 300
#' p <- 500
#' q <- 4
#' s <- 10 # number of non-zero coefficients
#' k <- sample(1:3, size=n, replace=TRUE)
#' 
#' sig2x <- 1
#' sig2z <- 1
#' sig2u <- 1
#' Sig_x <- choose_matrix(sig2x, rho=0.7, p=p, blk_sz=20, structure="block")
#' Sig_u <- choose_matrix(min(k)*sig2u, p=p, structure="diag")
#' Sig_z <- choose_matrix(sig2z, p=q, structure="diag")
#' Sig_xz <- matrix(0, nrow=p, ncol=q)
#' 
#' mu <- 100
#' mu_x <- rep(5, p)
#' mu_z <- rep(5, q)
#' 
#' # randomly sample sparsity index set.
#' S <- sample(1:p, s, replace=FALSE)
#' beta <- rep(0, p); beta[S] <- 2
#' gamma <- sample(1:3, size=q, replace=TRUE)
#' data <- data_gen(n=n, p=p, q=q, k=k,
#'                  mu=mu, mu_x=mu_x, mu_z=mu_z,
#'                  beta=beta, gamma=gamma,
#'                  Sig_x=Sig_x, Sig_u=Sig_u, Sig_z=Sig_z, Sig_xz=Sig_xz)
#' 
#' # average of replicates
#' W_bar <- t(sapply(data$W, colMeans, simplify=TRUE))
#' 
#' # cross-validated RC Lasso
#' cv.rcl.fit <- cv_rcl(data, Sig_x, Sig_u, Sig_z, Sig_xz,
#'                      rho=0.2,
#'                      family="gaussian", intercept=TRUE, standardize=TRUE, nfolds=10)
#' # for comparison, cross-validated Naive Lasso (i.e. no measurement error correction)
#' cv.nl.fit <- glmnet::cv.glmnet(cbind(W_bar, data$Z), data$y, family="gaussian",
#'                                intercept=TRUE, standardize=TRUE, nfolds=10)
#' # get estimated intercept from cross-validated 1se lambda
#' muRC_hat <- as.matrix(coef(cv.rcl.fit, s=cv.rcl.fit$lambda.1se))[1]
#' muNL_hat <- as.matrix(coef(cv.nl.fit, s=cv.nl.fit$lambda.1se))[1]
#' # get estimated beta coefficients from cross-validated 1se lambda
#' betaRC_hat <- as.matrix(coef(cv.rcl.fit, s=cv.rcl.fit$lambda.1se))[2:(p+1)]
#' betaNL_hat <- as.matrix(coef(cv.nl.fit, s=cv.nl.fit$lambda.1se))[2:(p+1)]
#' # get estimated gamma coefficients from cross-validated 1se lambda
#' gammaRC_hat <- as.matrix(coef(cv.rcl.fit, s=cv.rcl.fit$lambda.1se))[(p+2):(p+q+1)]
#' gammaNL_hat <- as.matrix(coef(cv.nl.fit, s=cv.nl.fit$lambda.1se))[(p+2):(p+q+1)]
#' @export

cv_rcl <- function(me_list, SCov_x, SCov_u, SCov_z=NULL, SCov_xz=NULL,
                   fun=NULL, ...){
  n <- length(me_list$y)
  p <- ncol(me_list$W[[1]])
  if (!is.null(me_list$Z)){
    q <- ncol(me_list$Z)
    mu_z_hat <- colMeans(me_list$Z)
  } else {
    q <- 0
    mu_z_hat <- NULL
  }
  
  unique_k <- sort(unique(me_list$k))
  
  # average of replicates
  W_bar <- t(sapply(me_list$W, colMeans, na.rm=TRUE, simplify=TRUE))
  
  # estimate other quantities
  mu_w_hat <- colSums( diag(me_list$k) %*% W_bar ) / sum(me_list$k)
  
  if(is.null(fun)) {
    if (n < p){
      warning("n < p, subsequent estimate may be nonsingular")
    }
    
    if (is.null(me_list$Z)){
      Lambda <- lapply(unique_k,
                       function(k){
                         solve(SCov_x + SCov_u/k) %*% SCov_x
                       }
      )
    } else {
      Lambda <- lapply(unique_k,
                       function(k){
                         solve(cbind(
                           rbind(SCov_x + SCov_u/k, t(SCov_xz)),
                           rbind(SCov_xz, SCov_z)
                         )) %*% rbind( SCov_x, t(SCov_xz) )
                       }
      )
    }
  } else {
    if (is.null(me_list$Z)){
      Lambda <- lapply(unique_k,
                       function(k){
                         fun(SCov_x + SCov_u/k, ...) %*% SCov_x
                       }
      )
    } else {
      Lambda <- lapply(unique_k,
                       function(k){
                         fun(cbind(
                           rbind(SCov_x + SCov_u/k, t(SCov_xz)),
                           rbind(SCov_xz, SCov_z)
                         ), ...) %*% rbind( SCov_x, t(SCov_xz) )
                       }
      )
    }
  }
  
  X_hat <- matrix(0, nrow=n, ncol=p)
  for (i in 1:n){
    X_hat[i,] <- ( t(mu_w_hat) - t(c(mu_w_hat, mu_z_hat)) %*%
                     Lambda[[ which(unique_k==me_list$k[i]) ]] ) +
      c(W_bar[i,], me_list$Z[i,]) %*% Lambda[[ which(unique_k==me_list$k[i]) ]]
  }
  
  return (
    glmnet::cv.glmnet(cbind(X_hat, me_list$Z), me_list$y, ...)
  )
}
