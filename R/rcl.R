#' Regression-Calibrated Lasso
#' 
#' @description
#' Lasso for linear or logistic regression with regression calibration correction for measurement error.
#' @param me_list List containing, at the very least, the number of replicates, k, the response y, and a list of k matrices of error-prone observations W -- one for each replicate.
#' @param Lambda Reliability matrix (or its estimate, if unknown).
#' @param ... Inputs to the glmnet function.
#' @return An object of class "glmnet".
#' @export

rcl <- function(me_list, Lambda, ...){
  # average of replicates
  W_bar <- Reduce('+', me_list$W) / me_list$k
  W_bar <- scale(W_bar, center=TRUE, scale=FALSE)
  
  X_hat <- W_bar %*% Lambda
  
  return (
    glmnet::glmnet(X_hat, me_list$y, ...)
  )
}