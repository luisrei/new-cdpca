#' CDpcaplot
#'
#' CDpcaplot plots CDpca object.
#'
#' @param data data frame (numeric).
#' @param obj CDpca object.
#' @param class Real classification for comparison (optional).
#'
#' @keywords cluster kmeans pca cdpca
#'
#' @export
#'
#' @examples
#' exampleplot = CDpcaplot(data, obj, class)

CDpcaplot <- function(data, obj, class = 0, save = 0, filename) {


  #verify save/filename relationship
  if (save == 1 && exists("filename")){
    if (!is.character(filename))   stop("Variable filename not of valid type.")
  }else if (save == 1) {
    stop("Filename not provided despite save option enabled.")
  }

  data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))  # normalized data, with variance divided by (I-1)
  I <- nrow(data)  # number of objects I

  if (class == 0 && length(class) == 1) {
    realclasstrue <- 2
    nb <- 10
    class0 <- numeric(nrow(data))
  }else{
    realclasstrue <- 3
    nb <- 11
    class0 <- data.frame(class)
    maxclass <- max(class0)
  }  # end if

  Ybar <- obj$originalYbar
  Yorder <- obj$Y
  Ybarorder <- obj$Ybar
  dorder <- obj$dorder
  Ucdpca <- obj$U

  # Table Real vs CDPCA classification
  classcdpca <- Ucdpca%*%as.matrix(1:ncol(Ucdpca))
  data.graphic <- data.frame(Yorder, class0, classcdpca, as.matrix(1:I))

  #
  # PLOTS: CDPCA classification
  #
  displaygraphics <- par(no.readonly = TRUE)
  par(mfrow = c(1, realclasstrue))
  if (realclasstrue == 3) {
    # plot1
    cl <- rainbow(maxclass)
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = "Dim 1", ylab = "Dim 2", type = "n", main = "Real classification")
    for (i in 1:maxclass) {
      points(data.graphic[, 1][class == i], data.graphic[, 2][class == i], pch = 1, col = cl[i])
    }
    # plot2
    cl <- rainbow(P)
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "CDPCA classification")
    for (i in 1:P) {
      points(data.graphic[, 1][classcdpca == i], data.graphic[, 2][classcdpca == i], pch = 2, col = cl[i])
    }
    points(Ybarorder[, 1], Ybarorder[, 2],pch = 15)  # introduce the centroids into the plot
    # plot3 legend
    matplot(Ybar, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("topleft",  title="Real class", pch = 1, col=c(cl[1:maxclass]), legend = c(1:maxclass), ncol = round(maxclass/2))
    legend("bottomleft",  title="CDPCA class", pch = 2, col=c(cl[1:P]), legend = 1:P, ncol = round(P/2))
  }else{
    # plot2
    cl <- rainbow(P)
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"),type = "n", main = "CDPCA classification")
    for (i in 1:P) {
      points(data.graphic[, 1][classcdpca == i], data.graphic[, 2][classcdpca == i], pch = 2, col = i+1)
    }
    points(Ybarorder[, 1],Ybarorder[, 2], pch = 15)  # introduce the centroids into the plot
    # plot3
    matplot(Ybar, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("bottomleft",  title="CDPCA class", pch = 2, col=c(cl[1:P]), legend = 1:P, ncol = round(P/2))
  }	 # end if



  if (save == 1){
    # identify file path
    if(substring(nb, 2) == 0){
      path <- file.path(getwd(),"CDpcaImages", paste(filename, "_noClass.jpeg", sep = ""))
    }else{
      path <- file.path(getwd(),"CDpcaImages", paste(filename, ".jpeg", sep=""))
    }
    # create desired directory
    dir.create(dirname(path), showWarnings = FALSE)
    # save dev to jpeg
    dev.copy(jpeg, path)
    dev.off()
  } #end if save


}  # end CDpcaplot function
