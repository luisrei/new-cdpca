#' NestedCDpca
#'
#' NestedCDpca performs (sequentially) two clustering and disjoint principal components analysis
#' on the given numeric data matrix and returns a list of results.
#'
#' @param data data frame (numeric).
#' @param method method to apply. CDpca or RCDpca.
#' @param class vector (numeric) or 0, if classes of objects are unknown.
#' @param fixAtt vector (numeric) or 0, for a selection of attributes.
#' @param nnloads 1 or 0, for nonnegative loadings
#' @param Q integer, number of clusters of variables.
#' @param P integer, number of clusters of objects.
#' @param K vector (numeric), representing number of sub-clusters.
#' @param ent integer, fuzzifier parameter (default 1).
#' @param tol real number, small positive tolerance.
#' @param maxit integer, maximum of iterations.
#' @param r number of runs of the cdpca algoritm for the final solution.
#' @param stand integer, 1 to standardize data before applying relaxed FKM (0 otherwise).
#'
#' @keywords cluster kmeans pca cdpca
#'
#' @return
#'  firstCDpca: CDpca object.
#'  secondCDpca: CDpca objects of second run of CDpca.
#'
#' @export
#'
#' @examples
#' example = NestedCDpca(data, method, class, fixAtt, nnloads, Q, P, K, ent, tol, maxit, r)

NestedCDpca <- function(data, method = "CDpca", class, fixAtt, nnloads = 0, Q, P, K, ent = 1, tol, maxit, r, stand) {


  # check arguments
  if(dim(array(K)) != P)    stop("Sub-clusters dimensions are not properly defined.")

  # initialize vector of lists carrying CDpca information
  secCdpca <- vector("list", P+1)

  I <- nrow(data)  # number of objects I
  #J <- ncol(data)  # number of variables J


  #
  #
  # CDPCA
  #
  #
  #
  #
  if (method == "CDpca"){

    # build upon first run of algorithm and analyse each object cluster
    for (subIter in 0:P){

      cat(sprintf("NestedCDPCA run #%d\n", subIter+1))

      if (class == 0 && length(class) == 1) {
        class0 <- numeric(I)
      }else{
        class0 <- class
      }  # end if

      # perform first CDpca run
      if (subIter == 0){
        if(P != 1){
          firstCdpca <- CDpca(data, class0, fixAtt, nnloads, Q, P, tol, maxit, r)
          cat("\n")
        }else{
          stop ("No clustering specified (P = 1).")
        }
      }else{
        if(K[subIter] != 1){
          secCdpca[[subIter]] <- CDpca(data[which(firstCdpca$U[, subIter] == 1), ], class0[which(firstCdpca$U[, subIter] == 1)], fixAtt, nnloads, Q, K[subIter], tol, maxit, r)
          if(subIter != P) cat("\n")
        }else {
          cat(sprintf("No clustering specified (P = 1).\n"))
          cat("-------------------\n")
          if(subIter != P) cat("\n")
          tab <- data.frame( class0[which(firstCdpca$U[, subIter] == 1)], rep.int(subIter, sum(firstCdpca$U[,subIter])))
          colnames(tab)<-c("Real Class", "CDPCA class")
          secCdpca[[subIter]] <- list(timeallloops = NULL,
                                      timebestloop = NULL,
                                      loop = NULL,
                                      iter = NULL,
                                      bcdev = NULL ,
                                      bcdevTotal = NULL,
                                      V = NULL,
                                      U = NULL,
                                      A = NULL,
                                      Fobj = NULL,
                                      Enorm = NULL,
                                      tableclass = tab,
                                      pseudocm = NULL,
                                      originalYbar = NULL,
                                      Y = NULL,
                                      Ybar = NULL,
                                      dorder = NULL,
                                      Dscale = NULL)
        }
      } # end if
    } # end for

    #
    #
    # Build second layer's classification vector
    #
    classcdpca2 <- numeric(I)
    for(subIter in 1:P){
      if(K[subIter] == 1){
        classcdpca2[which(firstCdpca$U[ , subIter] == 1 )] <-  max(classcdpca2)+1
      }else{
        for(subclust in 1:K[subIter]){
          classcdpca2[which(firstCdpca$U[ , subIter] == 1 )[which(secCdpca[[subIter]]$U[ , subclust] == 1 )]] <- max(classcdpca2)+1
        }
      }
    }

    secCdpca[[P+1]] <- classcdpca2

    #
    # OUTPUT
    #output <- new.env()
    #output$
    list(firstCDpca = firstCdpca,
         secondCDpca = secCdpca)
    #t <- as.list(output)




    #
    #
    #
    # RCDPCA
    #
    #
    #
    #
  }else if(method == "RCDpca"){
    ##CASE RELAXED CDPCA

    # build upon first run of algorithm and analyse each object cluster
    for (subIter in 0:P){

      cat(sprintf("NestedCDPCA run #%d\n", subIter+1))

      if (class == 0 && length(class) == 1) {
        class0 <- numeric(I)
      }else{
        class0 <- class
      }  # end if

      # perform first RCDpca run
      if (subIter == 0){
        if(P != 1){
          firstCdpca <- RelaxedCDpca(data, class0, fixAtt, nnloads, Q, P, ent, tol, maxit, r, stand)

          #
          #
          #Round fuzzy assignments
          U <- round(firstCdpca$U, 2)
          rownames(U) <- NULL
          colnames(U) <- NULL
          U2 <- matrix(0, I, P)
          for (i in 1:I) {
            U2[i, which.is.max(U[i,])] <- 1
          }  # end for

          cat("\n")
        }else{
          stop ("No clustering specified (P = 1).")
        }
      }else{

        if(K[subIter] != 1){
          secCdpca[[subIter]] <- RelaxedCDpca(data[which(U2[, subIter] == 1), ], class0[which(U2[, subIter] == 1)], fixAtt, nnloads, Q, K[subIter], ent, tol, maxit, r, stand)
          if(subIter != P) cat("\n")
        }else {
          cat(sprintf("No clustering specified (P = 1).\n"))
          cat("-------------------\n")
          if(subIter != P) cat("\n")

          tab <- data.frame( class0[which(U2[, subIter] == 1)], rep.int(subIter, sum(U2[,subIter])))
          colnames(tab)<-c("Real Class", "RCDPCA class")
          secCdpca[[subIter]] <- list(timeallloops = NULL,
                                      timebestloop = NULL,
                                      loop = NULL,
                                      iter = NULL,
                                      bcdev = NULL ,
                                      bcdevTotal = NULL,
                                      V = NULL,
                                      U = NULL,
                                      A = NULL,
                                      Fobj = NULL,
                                      Enorm = NULL,
                                      tableclass = tab,
                                      pseudocm = NULL,
                                      originalYbar = NULL,
                                      Y = NULL,
                                      Ybar = NULL,
                                      dorder = NULL,
                                      Dscale = NULL)
        }
      } # end if
    } # end for

    #
    #
    # Build second layer's classification vector
    #
    classcdpca2 <- numeric(I)
    for(subIter in 1:P){
      if(K[subIter] == 1){
        classcdpca2[which(U2[ , subIter] == 1)] <- max(classcdpca2)+1
      }else{
        for(subclust in 1:K[subIter]){
          classcdpca2[which(U2[ , subIter] == 1)[which(secCdpca[[subIter]]$U[ , subclust] == max(secCdpca[[subIter]]$U[ , subclust]) )]] <- max(classcdpca2)+1
        }
      }
    }

    secCdpca[[P+1]] <- classcdpca2

    #
    # OUTPUT
    #output <- new.env()
    #output$
    list(firstCDpca = firstCdpca,
         secondCDpca = secCdpca)
    #t <- as.list(output)

  }else{
    stop("Invalid method chosen.")
  }
}  # end NestedCDpca function
