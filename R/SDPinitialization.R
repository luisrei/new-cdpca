#' SDPinitialization
#'
#' Generates a binary and row stochastic matrix using SDP-based approach
#' 
#' @param dim1 number of rows
#' @param dim2 number of columns
#' 
#' @return The random matrix (dim1xdim2) with only one nonzero element per row.
#' 
#' @export
#' 
#' @examples
#' SDPinitialization(1, 2, data)

SDPinitialization <- function(dim1,dim2,D){
  em <- cbind(rep(1,dim1))
  matIem <- diag(dim1)-(1/dim1)*em%*%t(em)
  if (dim1 == dim(data)[1]){
    data <- D
  }else{
    data <- t(D)
  }  
  So <- matIem%*%data%*%t(data)%*%matIem
  matF <- matrix(svd(So)$v[, 1:(dim2-1)], ncol = dim2-1)
  matFFT <- matF%*%t(matF)
  # approximate solution to the SDP-based model for clustering
  Z.bar <- matFFT + 1/dim1 *em%*%t(em)
  
  # refine solutions
  cent <- Z.bar%*%data
  centunique <- unique(cent)
  t <- dim(centunique)[1]
  
  # initial kmeans step
  indcentroids <- sort(sample(1:t, dim2, replace=FALSE), decreasing=FALSE)
  centroids <- centunique[indcentroids, ] 
  solkmeans <- kmeans(data, centroids, iter.max = 100000)
  M <- matrix(0, dim1, dim2)
  for (i in 1:dim1){
    M[i, solkmeans$cluster[i]] <- 1
  }
  M
}

