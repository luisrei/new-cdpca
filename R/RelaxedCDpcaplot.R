#' RelaxedCDpcaplot
#'
#' RelaxedCDpcaplot plots RelaxedCDpca object.
#'
#' @param data data frame (numeric).
#' @param obj RelaxedCDpca object.
#' @param class Real classification for comparison (optional).
#'
#' @keywords cluster kmeans pca cdpca
#'
#' @export
#'
#' @examples
#' exampleplot = RelaxedCDpcaplot(data, obj, class)

RelaxedCDpcaplot <- function(data, obj, class = 0, save = 0, filename) {


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
  classcdpca <- obj$tableclass[, ncol(obj$tableclass)]
  data.graphic <- data.frame(Yorder, class0, classcdpca, as.matrix(1:I))

  #
  # PLOTS: CDPCA classification
  #
  displaygraphics <- par(no.readonly = TRUE)
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('orange', 'red', 'dark blue'))
  maxU <- apply(obj$U[, 1:ncol(obj$U)], 1, max)
  #This adds a column of color values based on the y values
  Cols <- rbPal(10)[as.numeric(cut(maxU,breaks = 10))]


  par(mfrow = c(1, realclasstrue))
  if (realclasstrue == 3) {
    # plot1
    cl <- rainbow(maxclass)
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = "Dim 1", ylab = "Dim 2", type = "n", main = "Real classification")
    for (i in 1:maxclass) {
      points(data.graphic[, 1][class == i], data.graphic[, 2][class == i], pch = 1, col = cl[i])
    }
    # plot2
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "CDPCA classification")
    for (i in 1:I) {
      points(data.graphic[i, 1], data.graphic[i, 2], pch = 1+classcdpca[i], col=Cols[i])
    }
    points(Ybarorder[, 1], Ybarorder[, 2],pch = 15)  # introduce the centroids into the plot
    # plot3 legend
    matplot(Ybar, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("topleft",  title="Real class", pch = 1, col=c(cl[1:maxclass]), legend = c(1:maxclass), ncol = round(maxclass/2))
    legend("bottomleft",  title="RCDPCA class", pch = c(2:(P+1)), col=c("black"), legend = 1:P, ncol = round(P/2))

  }else{
    # plot2
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "CDPCA classification")
    for (i in 1:I) {
      points(data.graphic[i, 1], data.graphic[i, 2], pch = 1+classcdpca[i], col=Cols[i])
    }
    points(Ybarorder[, 1], Ybarorder[, 2],pch = 15)  # introduce the centroids into the plot
    # plot3 legend
    matplot(Ybar, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("bottomleft",  title="RCDPCA class", pch = c(2:(P+1)), col=c("black"), legend = 1:P, ncol = round(P/2))
  }	 # end if



  if (save == 1){
    # identify file path
    if(substring(nb, 2) == 0){
      path <- file.path(getwd(),"RelaxedCDpcaImages", paste(filename, "_noClass.jpeg", sep = ""))
    }else{
      path <- file.path(getwd(),"RelaxedCDpcaImages", paste(filename, ".jpeg", sep=""))
    }
    # create desired directory
    dir.create(dirname(path), showWarnings = FALSE)
    # save dev to jpeg
    dev.copy(jpeg, path)
    dev.off()
  } #end if save


}  # end RelaxedCDpcaplot function
