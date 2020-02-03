#' Changes in composition of two transects
#'
#' Calculate the number of sampling points occupied by identical genotypes at
#' two time points.
#'
#' @param x,y Vectors of genotypes along a transect.
#'
#' @return A number between 0 and 1. 0 indicates that no sampling point is
#' occupied by the same genotype, 1 indicates that every sampling point is
#' occupied by the same genotype.
#' @export
#' @author Tom Ellis
transect_stability <- function(x,y){
  mean( x == y , na.rm = T)
}
