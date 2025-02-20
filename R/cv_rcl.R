#' Cross-Validated Regression-Calibrated Lasso
#' 
#' @description
#' Cross-validation for RC Lasso for linear or logistic regression.
#' @param me_list List containing, at the very least, the number of replicates, k, the response y, a list of k matrices of error-prone observations W -- one for each replicate, and a matrix of error-free predictors Z if one desires to include it in the model.
#' @param Lambda Reliability matrix (or its estimate, if unknown).
#' @param ... Inputs to the cv.glmnet function.
#' @return An object of class "cv.glmnet".
#' @examples
#' if(!require("glasso")){
#' stop("this example requires package 'glasso'")
#' }
#' 
#' set.seed(1424)
#' 
#' n <- 100
#' p <- 300
#' s <- 10 # number of non-zero coefficients
#' k <- 5
#' 
#' sig2x <- 1
#' sig2u <- 0.2
#' Sig_x <- choose_matrix(sig2x, rho=0.7, p=p, blk_sz=20, structure="block")
#' Sig_u <- choose_matrix(k*sig2u, rho=0, p=p, structure="diag")
#' 
#' # randomly sample sparsity index set.
#' S <- sample(1:p, s, replace=FALSE)
#' beta <- rep(0, p); beta[S] <- 2
#' data <- data_gen(n=n, p=p, k=k, beta=beta, Sig_x=Sig_x, Sig_u=Sig_u)
#' 
#' # average of replicates
#' W_bar <- Reduce('+', data$W) / k
#' W_bar <- scale(W_bar, center=TRUE, scale=FALSE)
#' # use Graphical Lasso (Friedman et al., 2008)
#' gl <- glasso::glasso(cov(W_bar), rho=0.2)
#' 
#' # estimate Sig_u, a diagonal matrix
#' svar2u <- rep(0, p)
#' for (i in 1:p){
#'   for (j in 1:k){
#'     svar2u[i] <- svar2u[i] + sum( (data$W[[j]][,i] - W_bar[,i])^2 )
#'   }
#' }
#' svar2u <- svar2u / (n*(k-1))
#' SCov_uu <- diag(svar2u)
#' 
#' # estimate reliability matrix
#' Lambda <- gl$wi %*% (gl$w - SCov_uu/k)
#' 
#' # cross-validated RC Lasso
#' cv.rcl.fit <- cv_rcl(data, Lambda, family="gaussian", intercept=FALSE,
#'                      standardize=TRUE, nfolds=10)
#' # for comparison, cross-validated Naive Lasso (i.e. no measurement error correction)
#' cv.nl.fit <- glmnet::cv.glmnet(W_bar, data$y, family="gaussian", intercept=FALSE,
#'                        standardize=TRUE, nfolds=10)
#' # get estimated coefficients (modulo intercept) from cross-validated 1se lambda
#' betaRC_hat <- as.matrix(coef(cv.rcl.fit, s=cv.rcl.fit$lambda.1se))[2:(p+1)]
#' betaNL_hat <- as.matrix(coef(cv.nl.fit, s=cv.nl.fit$lambda.1se))[2:(p+1)]
#' @export

cv_rcl <- function(me_list, Lambda, ...){
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
  W_bar <- scale(W_bar, center=TRUE, scale=FALSE)
  
  mu_x_hat <- colMeans(W_bar)
  
  X_hat <- rep(1,n) %*% t(c(mu_x_hat, mu_z_hat)) %*% (diag(p+q) - Lambda) +
    cbind(W_bar, me_list$Z) %*% Lambda
  
  return (
    glmnet::cv.glmnet(cbind(X_hat, me_list$Z), me_list$y, ...)
  )
}