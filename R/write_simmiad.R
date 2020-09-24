#' Write to disk
#'
#' Write the output of simmiad simulations to disk
#'
#' @param x Output of the function `simmiad()`
#' @param directory Directory where output should be saved without trailing
#' slash.
#'
#' @author Tom Ellis
#' @export
write_simmiad <- function(x, directory){
  dir.create(directory, showWarnings = FALSE)

  for(i in names(x)){
    write.csv(
      x[[i]],
      file = paste(directory, "/", i, '.csv', sep=""),
      row.names = F)
  }
}
