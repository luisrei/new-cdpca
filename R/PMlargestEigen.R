#' PMlargestEigen
#'
#' The Power Method simple function
#'
#' @param M covariance matrix restricted to some attributes (symmetric)
#'
#' @export

######################################
#
# The Power Method simple function
#
######################################
PMlargestEigen <- function(M){
  maxiteig <- 10000
  #x0 <- cbind(rep(1,nrow(M)))
  x0 <- cbind(runif(nrow(M))) #; x0
  it <- 0
  repeat
  {
    x1 = M %*% x0 ;
    x1 = x1 / norm(x1, "F") #; x1; x1[1] <- -1 ; x1 = x1 / norm(x1, "F"); x1
    #
    #NEW nonnegative loadings
    #
    if (nnloads ==1){
      nna <- which(x1 >= 0) #; nna # which indices are >=0 ?
      if ( length(nna) < length(x1) ){
        # se existem elementos <0 vamos calcular os elementos para os índices com >=0
        x1rednna <- as.matrix(x1[nna]) #; x1rednna
        x1nna <- M[nna, nna] %*% x1rednna #; x1nna
        x1[nna] <- x1nna %*% diag(1/sqrt(t(x1nna)%*% x1nna)) # normalizar vetor
        # posições com índices <0 ficam com 0
        na <- which( x1 < 0 )
        x1[na] <- 0
      }
    }
    if (norm(abs(x1 - x0), type=c("F")) <= tol ) break
    x0 <- x1
    it <- it + 1
    if (it == maxiteig) break
  }
  lambda = sum((M %*% x1) * x1)
  list(vector = x1, value = lambda, iterations = it)
}
