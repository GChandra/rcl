#' Choose Matrix
#' 
#' @description
#' Creates a covariance matrix, either diagonal or block-diagonal, with the given parameters. The covariance matrix has identical variance across the diagonal.
#' @param sig2 Variance (identical for all predictors).
#' @param rho Correlation used for blocks in block-diagonal structure.
#' @param p Number of features.
#' @param blk_sz Size of blocks in block-diagonal structure (default 20x20).
#' @param structure Covariance matrix structure (either "diag" for diagonal matrix, or "block" for block-diagonal matrix).
#' @return A p x p matrix.
#' @export

choose_matrix <- function(sig2, rho=0.7, p, blk_sz = 20, structure){
  if(structure == "diag") {
    mat <- matrix(0, p, p)
    diag(mat) <- sig2
  } else if(structure == "block") {
    blk <- matrix(NA, blk_sz, blk_sz)
    for (i in 1:blk_sz){
      for (j in 1:blk_sz){
        blk[i,j] <- sig2 * rho^(abs(i-j))
      }
    }
    mat <- as.matrix( Matrix::bdiag(replicate(p/blk_sz, blk, simplify=FALSE)) )
  } else {
    mat <- matrix(NA, p, p)
  }
  
  return(mat)
}