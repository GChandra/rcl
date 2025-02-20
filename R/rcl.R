#' Regression-Calibrated Lasso
#' 
#' @description
#' Lasso for linear or logistic regression with regression calibration correction for measurement error.
#' @param me_list List containing, at the very least, the number of replicates, k, the response y, a list of k matrices of error-prone observations W -- one for each replicate, and a matrix of error-free predictors Z if one desires to include it in the model.
#' @param Lambda Reliability matrix (or its estimate, if unknown).
#' @param ... Inputs to the glmnet function.
#' @return An object of class "glmnet".
#' @export

rcl <- function(me_list, Lambda, ...){
  n <- length(me_list$y)
  p <- ncol(me_list$W[[1]])
  if (!is.null(me_list$Z)){
    q <- ncol(me_list$Z)
    mu_z_hat <- colMeans(me_list$Z)
  } else {
    q <- 0
    mu_z_hat <- NULL
  }
  
  # average of replicates
  W_bar <- Reduce('+', me_list$W) / me_list$k
  
  mu_x_hat <- colMeans(W_bar)
  
  X_hat <- rep(1,n) %*% t(c(mu_x_hat, mu_z_hat)) %*% (diag(p+q) - Lambda) +
    cbind(W_bar, me_list$Z) %*% Lambda
  
  return (
    glmnet::glmnet(cbind(X_hat, me_list$Z), me_list$y, ...)
  )
}