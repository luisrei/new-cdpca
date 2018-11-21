#' RandMat
#'
#' Generates a random binary and row stochastic matrix.
#' 
#' @param dim1 number of rows
#' @param dim2 number of columns
#' 
#' @return The random matrix (dim1xdim2) with only one nonzero element per row.
#' 
#' @export
#' 
#' @examples
#' RandMat(3, 5)

RandMat <- function(dim1, dim2) {
  U0 <- matrix(0, dim1, dim2)
  U0[1:dim2, 1:dim2] <- diag(dim2)
  for (j in (dim2+1):dim1){
    p <- sample(1:dim2, 1)
    U0[j, p] <- 1
  }
  U0
}  # end RandMat function