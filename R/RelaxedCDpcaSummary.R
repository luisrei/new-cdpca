#' RelaxedCDpcaSummary
#'
#' Print a variety of useful information regarding a RelaxedCDpca object.
#'
#' @param obj RelaxedCDpca object
#'
#' @export

RelaxedCDpcaSummary <- function(obj) {

  if(length(obj)!=18) stop("Invalid RelaxedCDpca object.")

  cat("Summary of Relaxed Cluster and Disjoint Principal Component Analysis:\n")
  print(rbind("Loop:" = obj$loop, "Number of iterations:" = obj$iter,
              "F:" = round(obj$Fobj, 5), "Error:" = round(obj$Enorm, 5),
              "Between cluster deviance (%):" = round(obj$bcdev, 2)))
  print(rbind("Explained variance by RelaxedCDpca components" = obj$dorder))
  cat("Pseudo Confusion Matrix:\n")
  print(obj$pseudocm)
  invisible(obj)

}  # end RelaxedCDpcaSummary function
