#' Regression-Calibrated Lasso
#' 
#' @description
#' Lasso for linear or logistic regression with regression calibration correction for measurement error.
#' @param me_list List containing, at the very least, an n-dimensional vector indicating the number of replicates per subject, k, the response y, a list of n matrices of error-prone observations W -- one for each subject and with dimension k[i] x p, and a matrix of error-free predictors Z if one desires to include it in the model.
#' @param SCov_x p x p dimensional matrix for the estimated covariance matrix for true predictors X.
#' @param SCov_u p x p dimensional matrix for the estimated covariance matrix for the measurement errors U.
#' @param SCov_z q x q dimensional matrix for the estimated covariance matrix for the error-free predictors Z (optional).
#' @param SCov_xz p x q dimensional matrix for the estimated cross-covariance matrix between X and Z (optional).
#' @param fun function used to estimate precision matrices. Must take as the first input the matrix to invert, followed by additional arguments.
#' @param ... Additional inputs to "fun" and to the glmnet function.
#' @return An object of class "glmnet".
#' @export

rcl <- function(me_list, SCov_x, SCov_u, SCov_z=NULL, SCov_xz=NULL,
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
  W_bar <- t(sapply(me_list$W, colMeans, simplify=TRUE))
  
  # estimate other quantities
  mu_w_hat <- colSums( diag(me_list$k) %*% W_bar ) / sum(me_list$k)
  mu_z_hat <- colMeans(me_list$Z)
  
  if(is.null(fun)) {
    if (n < p){
      warning("n < p, subsequent estimate may be nonsingular")
    }
    
    Lambda <- lapply(unique_k,
                     function(k){
                       solve(cbind(
                         rbind(SCov_x + SCov_u/k, t(SCov_xz)),
                         rbind(SCov_xz, SCov_z)
                       )) %*% rbind( SCov_x, t(SCov_xz) )
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
  
  X_hat <- matrix(0, nrow=n, ncol=p)
  for (i in 1:n){
    X_hat[i,] <- ( t(mu_w_hat) - t(c(mu_w_hat, mu_z_hat)) %*% Lambda[[ which(unique_k==k[i]) ]] ) +
      c(W_bar[i,], me_list$Z[i,]) %*% Lambda[[ which(unique_k==k[i]) ]]
  }
  
  return (
    glmnet::glmnet(cbind(X_hat, me_list$Z), me_list$y, ...)
  )
}