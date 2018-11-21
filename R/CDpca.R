#' CDpca
#'
#' CDpca performs a clustering and disjoint principal components analysis
#' on the given numeric data matrix and returns a list of results.
#'
#' @param data data frame (numeric).
#' @param class vector (numeric) or 0, if classes of objects are unknown.
#' @param fixAtt vector (numeric) or 0, for a selection of attributes.
#' @param nnloads 1 or 0, for nonnegative loadings.
#' @param Q integer, number of clusters of variables.
#' @param P integer, number of clusters of objects.
#' @param tol real number, small positive tolerance.
#' @param maxit integer, maximum of iterations.
#' @param r number of runs of the cdpca algoritm for the final solution.
#'
#' @keywords cluster kmeans pca cdpca
#'
#' @return
#'  iter: iterations used in the best loop for computing the best solution
#'  loop: best loop number
#'  timebestloop: computation time on the best loop
#'  timeallloops: computation time for all loops
#'  Y: the component score matrix
#'  Ybar: the object centroids matrix in the reduced space
#'  A: the component loading matrix
#'  U: the partition of objects
#'  V: the partition of variables
#'  Fobj: function to maximize
#'  bcdev: between cluster deviance
#'  bcdevTotal: between cluster deviance over the total variability
#'  tableclass: cdpca classification
#'  pseudocm: pseudo confusion matrix of the real and cdpca classifications
#'  Enorm: error norm for the obtained cdpca model
#'
#' @export
#'
#' @examples
#' exampleCDpca = CDpca(data, class, fixAtt, nnloads, Q, P, tol, maxit, r)

CDpca <- function(data, class, fixAtt, nnloads = 0, Q, P, tol, maxit, r) {

  if(P == 1)    stop("No clustering specified (P = 1).")

  data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))  # normalized data
  #data.sd <- data.frame(apply(data, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))

  # in the paper: standardized to have mean zero and unit variance
  I <- nrow(data)  # number of objects I
  J <- ncol(data)  # number of variables J
  Xs <- as.matrix(data.sd*sqrt(I/(I-1)))  # matrix of normlized data (with var divided by I)
  # matriz X (I x J) (obs x vars), numeric matrix argument
  tfbest <- matrix(0, r, 1)  # computational time at each loop
  nnloads <- nnloads

  cat(sprintf("CDPCA running...\n"))
  # Run CDPCA ALS algorithm r times
  for (loop in 1:r) {
    t1 <- proc.time()

    # Initialization
    iter <- 0
    U <- RandMat(I, P)
    V <- RandMat(J, Q)

    #
    # NEW
    #
    # if there is no attribute selection else, with selection
    if (fixAtt == 0 && length(fixAtt) == 1){
      indexSelectFree <- 1:J
    }else {
      indexSelect <- which(fixAtt>0)
      for (i in 1:length(indexSelect)){
        V[indexSelect[i], ] <- 0
        V[indexSelect[i], fixAtt[indexSelect[i]]] <- 1
      }
      indexSelectFree <- which(fixAtt==0)
    }
    #
    #
    #
    X.bar <- diag(1/colSums(U))%*%t(U)%*%Xs  # (PxJ) object centroid matrix. It identifies the P centroids in the variable space
    X.group <- U%*%X.bar  # (IxJ) matrix. Identify the I objects by their centroids.
    vJ <- matrix(1:J, ncol = 1)
    zJ <- matrix(0, J, 1)
    zQ <- matrix(0, 1, Q)
    # matrix A
    A <- matrix(0, J, Q)
    for (j in 1:Q) {
      Jxx <- which(V[, j] == 1)
      len.Jxx <- length(Jxx)
      if (sum(V[, j])>1) {
        if (I >= len.Jxx) {
          # CASE 1: I >= J (more observations than variables)
          S.group <- t(X.group[, Jxx])%*%X.group[, Jxx]
          #
          #NEW
          #
          #         A[Jxx, j] <- eigen(S.group)$vectors[, 1]
          A[Jxx, j] <- PMlargestEigen(S.group)$vector
        }else{
          # CASE 2: I < J ( more variables than observations)
          SS.group <- X.group[, Jxx]%*%t(X.group[, Jxx])
          ##          A[Jxx, j] <- t(X.group[, Jxx])%*%eigen(SS.group)$vectors[, 1]/sqrt(eigen(SS.group)$values[1])
          PMinitial <- PMlargestEigen(SS.group)
          A[Jxx, j] <- t(X.group[, Jxx])%*%PMinitial$vector/sqrt(PMinitial$value)
        }  # end if
      }else{
        A[Jxx, j] <- 1
      }  # end if
    }  # end for
    A0 <- A
    Y.bar <- X.bar%*%A  # (PxQ) object centroid matrix. It identifies the P centroids in the reduced space of the principal components.
    F0 <- sum(diag(t(U%*%Y.bar)%*%U%*%Y.bar))  # Objective function to maximize.
    Fmax <- F0  ##########################################################################AQUI
    conv <- 2*tol
    while (conv > tol) {
      iter <- iter+1
      Y <- Xs%*%A  # (IxQ) component score matrix. It identifies the I objects in the reduced space.
      U <- matrix(0, I, P)
      # update U
      for (i in 1:I) {
        dist <- rep(0, P)
        for (p in 1:P) {
          dist[p] <- sum((Y[i, ]-Y.bar[p, ])^2)
        } # end for
        min.dist <- which.min(dist)
        U[i, min.dist] <- 1
      }  # end for
      su <- colSums(U)
      while (sum(su == 0) > 0) {
        ind.max <- which.max(su)  # p2
        su.max <- max(su)  # m2
        ind.min <- which.min(su)  # p1
        su.min <- min(su)  # m1
        ind.nzU <- which(U[, ind.max] == 1)
        ind.sel <- ind.nzU[1:floor(su.max)/2]
        U[ind.sel, ind.min] <- 1
        U[ind.sel, ind.max] <- 0
        su <- colSums(U)
      }  # end while
      # Given U and A compute X.bar
      X.bar <- diag(1/colSums(U))%*%t(U)%*%Xs
      Y.bar <- X.bar%*%A
      X.group <- U%*%X.bar
      # Update V and A
      #
      # NEW
      #
      ##indexSelectFree <- which(fixAtt==0)
      for (j in indexSelectFree ){

        # for (j in 1:J) {
        posmax <- which(V[j, ] == 1)
        for (g in 1:Q) {
          #################################
          V[j, ] <- diag(Q)[g, ]
          if (g!=posmax) {  ### if 1
            for (gg in c(g,posmax)) {
              xx <- V[, gg]
              Jxx <- which(xx == 1)
              len.Jxx <- length(Jxx)
              A[, gg] <- zJ
              if (sum(xx) > 1) {
                if (I >= len.Jxx) {
                  # CASE 1: I >= J (more observations than variables)
                  S.group <- t(X.group[, Jxx])%*%X.group[, Jxx]
                  #              	A[Jxx, gg] <- matrix(eigen(S.group)$vectors[, 1], nrow = len.Jxx)
                  A[Jxx, gg] <- matrix(PMlargestEigen(S.group)$vector, nrow = len.Jxx)
                }else{
                  # CASE 2: I < J (more variables than observations)
                  SS.group <- X.group[, Jxx]%*%t(X.group[, Jxx])
                  ##              	A[Jxx, gg] <- matrix((t(X.group[, Jxx])%*%eigen(SS.group)$vectors[, 1]/sqrt(eigen(SS.group)$values[1]))[, 1], nrow = len.Jxx)
                  PMgeneral <- PMlargestEigen(SS.group)
                  A[Jxx, gg] <- matrix((t(X.group[, Jxx])%*%PMgeneral$vector/sqrt(PMgeneral$value))[, 1], nrow = len.Jxx)
                }  # end if
              }else{
                if (sum(xx)==1) {
                  A[Jxx, gg] <- 1
                }  # end if
              }  # end if
            }  # end for
          }  # end ### if 1
          ##################################
          Y.bar <- X.bar%*%A
          Fobj <- sum(diag(t(U%*%Y.bar)%*%U%*%Y.bar))
          if (Fobj > Fmax) {
            Fmax <- Fobj
            posmax <- g
            A0 <- A
          }else{
            A <- A0
          }  # end if
        }  # end for
        V[j, ]=diag(Q)[posmax, ]
      }  # end for
      Y <- Xs%*%A
      Y.bar <- X.bar%*%A
      Fobj <- sum(diag(t(U%*%Y.bar)%*%U%*%Y.bar))
      conv <- Fobj-F0
      if (conv > tol) {
        F0 <- Fobj
        A0 <- A
      }else{
        break
      }  # end if
      if (iter == maxit) {
        print("Maximum iterations reached.")
        break
      }  # end if
    }  # end while
    t2 <- proc.time()
    tfinal <- t2-t1
    tfbest[loop] <- tfinal[3]  # Computation time for each loop
    # BetClusDevTotal <- Fobj/(I*J)*100  # between cluster deviance
    BetClusDev <- (Fobj/sum(diag(t(Y)%*%Y)))*100  # between cluster deviance
    tabperloop <- data.frame(cbind(loop, iter, tfinal[3], BetClusDev, Fobj, conv))
    rownames(tabperloop) <- c(" ")
    colnames(tabperloop) <- c("Loop", "Iter", "Loop time", "Between cluster deviance(%):", "F", "Convergence")
    print(tabperloop)  # Loop display
    if (loop == 1) {
      Vcdpca <- V
      Ucdpca <- U
      X.barcdpca <- diag(1/colSums(Ucdpca))%*%t(Ucdpca)%*%Xs
      Acdpca <- A
      Ycdpca <- Xs%*%Acdpca
      Y.barcdpca <- X.barcdpca%*%Acdpca
      Fcdpca <- Fobj
      loopcdpca <- 1
      itercdpca <- iter
      convcdpca <- conv
    }  # end if
    if (Fobj > Fcdpca) {
      Vcdpca <- V
      Ucdpca <- U
      X.barcdpca <- diag(1/colSums(Ucdpca))%*%t(Ucdpca)%*%Xs
      Acdpca <- A
      Ycdpca <- Xs%*%Acdpca
      Y.barcdpca <- X.barcdpca%*%Acdpca
      Fcdpca <- Fobj
      loopcdpca <- loop
      itercdpca <- iter
      convcdpca <- conv
    }  # end if
  }  # end for loop

  cat("-------------------\n")

  tftotal=sum(tfbest)  # Computation time for all loops
  # maximum between cluster deviance
  BetClusDevTotal <- Fcdpca/(I*J)*100  # (sum(diag(var(Y)))*(I-1))*100
  BetClusDev <- (sum(diag(t(Ucdpca%*%Y.barcdpca)%*%(Ucdpca%*%Y.barcdpca)))/sum(diag(t(Ycdpca)%*%Ycdpca)))*100
  # error in cdpca model
  Enormcdpca <- 1/I*norm(Xs-Ucdpca%*%Y.barcdpca%*%t(Acdpca), type = c("F"))
  # Table Real vs CDPCA classification
  classcdpca <- Ucdpca%*%as.matrix(1:ncol(Ucdpca))
  class <- data.frame(class)
  maxclass <- max(class)
  # Variability Ycdpca and decreasing order of PC
  varYcdpca <- var(Ycdpca)
  d <- round(diag(varYcdpca)*100/J, 2)
  dorder <- d[order(d, decreasing = TRUE)]
  tablevar <- data.frame(1:Q, d)
  colnames(tablevar) <- c("Dim", "Var (%)")
  #
  # Presentation of the matrices associated to the CDPCA component sortedby their variances.
  #
  Yorder <- t(t(rbind(d, Ycdpca))[order(d, decreasing = TRUE), ])[-1, ]
  Ybarorder <- t(t(rbind(d, Y.barcdpca))[order(d, decreasing=TRUE),])[-1, ]
  Aorder <- t(t(rbind(d, Acdpca))[order(d, decreasing = TRUE), ])[-1, ]
  Vorder <- t(t(rbind(d, Vcdpca))[order(d, decreasing = TRUE), ])[-1, ]
  #
  # We can check the model using this colunm sort matrices and observe that
  #
  # Ucdpca%*%Y.barcdpca%*%t(Acdpca) = Ucdpca%*%Ybarorder%*%t(Aorder)

  if (class == 0 && length(class) == 1) {
    realclasstrue <- 2
    pseudocm <- NULL
    tabrealcdpca <- data.frame(classcdpca)
    colnames(tabrealcdpca) <- c("CDPCA Class")
  }else{
    realclasstrue <- 3
    tabrealcdpca <- data.frame(class, classcdpca)
    colnames(tabrealcdpca) <- c("Real Class", "CDPCA Class")
    # pseudo-confusion matrix
    pseudocm <- table(tabrealcdpca)
  }  # end if

  # #
  # # PLOTS: CDPCA classification
  # #
  # data.graphic <- data.frame(Yorder, class, classcdpca, as.matrix(1:I))
  # if (cdpcaplot == 1) {
  #   displaygraphics <- par(no.readonly = TRUE)
  #   par(mfrow = c(1, realclasstrue))
  #   if (realclasstrue == 3) {
  #     # plot1
  #     matplot(data.graphic[, 1], data.graphic[, 2], xlab = "Dim 1", ylab = "Dim 2", type = "n", main = "Real classification")
  #     for (i in 1:maxclass) {
  #       points(data.graphic[, 1][class == i], data.graphic[, 2][class == i], pch = 1, col = i+1)
  #     }
  #     # plot2
  #     matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "CDPCA classification")
  #     for (i in 1:P) {
  #       points(data.graphic[, 1][classcdpca == i], data.graphic[, 2][classcdpca == i], pch = 2, col = i+1)
  #     }
  #     points(Ybarorder[, 1], Ybarorder[, 2],pch = 15)  # introduce the centroids into the plot
  #     # plot3 legend
  #     matplot(Y.barcdpca, type = "n", axes = FALSE, xlab = "", ylab = "")
  #     legend("topleft", pch = 1, col=c(2:(maxclass+1)), legend = c(1:maxclass), ncol = maxclass)
  #     legend("bottomleft", pch = 2, col=c(2:(P+1)), legend = 1:P, ncol = P)
  #   }else{
  #     # plot2
  #     matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"),type = "n", main = "CDPCA classification")
  #     for (i in 1:P) {
  #       points(data.graphic[, 1][classcdpca == i], data.graphic[, 2][classcdpca == i], pch = 2, col = i+1)
  #     }
  #     points(Ybarorder[, 1],Ybarorder[, 2], pch = 15)  # introduce the centroids into the plot
  #     # plot3
  #     matplot(Y.barcdpca, type = "n", axes = FALSE, xlab = "", ylab = "")
  #     legend("bottomleft", pch = 2, col=c(2:(P+1)), legend = 1:P, ncol = P)
  #   }	 # end if
  #   if (save == 1){
  #     # identify file path
  #     if (class == 0 && length(class) == 1) {
  #       path <- file.path(getwd(),"CDpcaImages", paste(filename, ".jpeg", sep = ""))
  #     }else{
  #       path <- file.path(getwd(),"CDpcaImages", paste(filename, "_noClass.jpeg", sep = ""))
  #     }
  #     # create desired directory
  #     dir.create(dirname(path), showWarnings = FALSE)
  #     # save dev to jpeg
  #     dev.copy(jpeg, path)
  #     dev.off()
  #   } #end if save
  # } #end if cdpcaplot


  # OUTPUT

  #output <- new.env()
  #output$
  list(timeallloops = tftotal,
       timebestloop = tfbest[loopcdpca],
       loop = loopcdpca,
       iter = itercdpca,
       bcdevTotal = BetClusDevTotal ,
       bcdev = BetClusDev ,
       V = Vorder,
       U = Ucdpca,
       A = Aorder,
       Fobj = Fcdpca,
       Enorm = Enormcdpca,
       tableclass = tabrealcdpca,
       pseudocm = pseudocm,
       originalYbar = Y.barcdpca,
       Y = Yorder,
       Ybar = Ybarorder,
       dorder = dorder,
       Dscale = Xs)
  #   t <- as.list(output)
}  # end CDpca function
