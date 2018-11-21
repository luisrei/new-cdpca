#' NestedCDpcaSummary
#'
#' Print a variety of useful information regarding a NestedCDpca object.
#'
#' @param obj NestedCDpca object
#'
#' @export

NestedCDpcaSummary <- function(obj) {

  if(length(obj)!=2) stop("Invalid NestedCDpca object.")

  cat(sprintf("CDPCA run #%d\n", 1))
  CDpcaSummary(obj[[1]])

  #Summary of second layer
  for (i in 1:(length(obj[[2]])-1)) {

    cat("\n-----------------------------\n\n")

    cat(sprintf("CDPCA run #%d\n", i+1))
    objs <- (obj[[2]][[i]])

    if (is.null(objs$loop)){
      cat("No clusters were requested in this run.\n")
    }else{
      cat("Summary of Nested Cluster and Disjoint Principal Component Analysis:\n")
      print(rbind("Loop:" = objs$loop, "Number of iterations:" = objs$iter,
                  "F:" = round(objs$Fobj, 5), "Error:" = round(objs$Enorm, 5),
                  "Between cluster deviance (%):" = round(objs$bcdev, 2)))
      print(rbind("Explained variance by CDpca components" = objs$dorder))
    }
  }
  cat("\n")
  invisible(obj)

}  # end NestedCDpcaSummary function
