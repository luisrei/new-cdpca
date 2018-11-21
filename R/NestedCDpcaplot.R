#' NestedCDpcaplot
#'
#' NestedCDpcaplot plots NestedCDpca object.
#'
#' @param data data frame (numeric).
#' @param obj NestedCDpca object.
#' @param class Real classification for comparison (optional).
#'
#' @keywords cluster kmeans pca cdpca
#'
#' @export
#'
#' @examples
#' exampleplot = NestedCDpcaplot(data, obj, class)

NestedCDpcaplot <- function(data, obj, class = 0, save = 0, filename) {


  #verify save/filename relationship
  if (save == 1 && exists("filename")){
    if (!is.character(filename))   stop("Variable filename not of valid type.")
  }else if (save == 1) {
    stop("Filename not provided despite save option enabled.")
  }

  data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))  # normalized data, with variance divided by (I-1)
  I <- nrow(data)  # number of objects I

  #dual cdpca object
  if (class == 0 && length(class) == 1) {
    realclasstrue <- 3
    nb <- 20
    class0 <- numeric(nrow(data))
  }else{
    realclasstrue <- 4
    nb <- 21
    class0 <- data.frame(class)
    maxclass <- max(class0)
  }  # end if

  Ybarorder <- obj[[1]]$Ybar
  Ybar <- obj[[1]]$originalYbar
  Yorder <- obj[[1]]$Y

  dorder <- obj[[1]]$dorder
  classcdpca <- obj[[1]]$tableclass[,ncol(obj[[1]]$tableclass)]
  classcdpca2 <- obj[[2]][[length(obj[[2]])]]
  data.graphic <- data.frame(Yorder, class0, classcdpca, classcdpca2, as.matrix(1:I))

  #
  # PLOTS: CDPCA classification
  #
  displaygraphics <- par(no.readonly = TRUE)
  par(mfrow = c(1, realclasstrue))
  if (realclasstrue == 4) {
    # plot1
    cl <- rainbow(maxclass)
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = "Dim 1", ylab = "Dim 2", type = "n", main = "Real class")
    for (i in 1:maxclass) {
      points(data.graphic[, 1][class == i], data.graphic[, 2][class == i], pch = 1, col = cl[i])
    }
    # plot2
    cl <- rainbow(max(classcdpca))
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "1st CDPCA class")
    for (i in 1:max(classcdpca)) {
      points(data.graphic[, 1][classcdpca == i], data.graphic[, 2][classcdpca == i], pch = 2, col = cl[i])
    }
    points(Ybarorder[, 1], Ybarorder[, 2], pch = 15)  # introduce the centroids into the plot
    # plot3
    cl <- rainbow(max(classcdpca2))
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "2nd CDPCA class")
    for (i in 1:max(classcdpca2)) {
      points(data.graphic[, 1][classcdpca2 == i], data.graphic[, 2][classcdpca2 == i], pch = 3, col = cl[i])
    }
    #points(Ybarorder2[, 1], Ybarorder2[, 2],pch = 15)  # introduce the centroids into the plot
    # plot4 legend
    matplot(Ybar, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("topleft", title="Real class", pch = 1, col=c(cl[1:maxclass]), legend = c(1:maxclass), ncol = round(maxclass/2))
    legend("left", title="1st CDPCA", pch = 2, col=c(cl[1:max(classcdpca)]), legend = 1:max(classcdpca), ncol = round(max(classcdpca)/2))
    legend("bottomleft", title="2nd CDPCA", pch = 3, col=c(cl[1:max(classcdpca2)]), legend = 1:max(classcdpca2), ncol = round(max(classcdpca2)/2))
  }else{
    # plot2
    cl <- rainbow(max(classcdpca))
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "1st CDPCA class")
    for (i in 1:max(classcdpca)) {
      points(data.graphic[, 1][classcdpca == i], data.graphic[, 2][classcdpca == i], pch = 2, col = cl[i])
    }
    points(Ybarorder[, 1], Ybarorder[, 2], pch = 15)  # introduce the centroids into the plot
    # plot3
    cl <- rainbow(max(classcdpca2))
    matplot(data.graphic[, 1], data.graphic[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "2nd CDPCA class")
    for (i in 1:max(classcdpca2)) {
      points(data.graphic[, 1][classcdpca2 == i], data.graphic[, 2][classcdpca2 == i], pch = 3, col = cl[i])
    }
    #points(Ybarorder2[, 1], Ybarorder2[, 2],pch = 15)  # introduce the centroids into the plot
    # plot4 legend
    matplot(Ybar, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("topleft", title="1st CDPCA class", pch = 2, col=c(cl[1:max(classcdpca)]), legend = 1:max(classcdpca), ncol = round(max(classcdpca)/2))
    legend("bottomleft", title="2nd CDPCA class", pch = 3, col=c(cl[1:max(classcdpca2)]), legend = 1:max(classcdpca2), ncol = round(max(classcdpca2)/2))
  }	 # end if

  if (save == 1){
    # identify file path
    if(substring(nb, 2) == 0){
      path <- file.path(getwd(),"NestedCDpcaImages", paste(filename, "_noClass.jpeg", sep = ""))
    }else{
      path <- file.path(getwd(),"NestedCDpcaImages", paste(filename, ".jpeg", sep=""))
    }
    # create desired directory
    dir.create(dirname(path), showWarnings = FALSE)
    # save dev to jpeg
    dev.copy(jpeg, path)
    dev.off()
  } #end if save


}  # end NestedCDpcaplot function
