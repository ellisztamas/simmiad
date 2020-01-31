#' Estimate spatial clustering
#'
#' Estimate the degree of spatial clustering along a transect
#'
#' `transect_clustering` estimates the degree of non-randomness in the
#' distribution of genotypes along a linear transect by comparing the mean
#' pairwise distances between sampling points with the same genotype to mean
#' pairwise distances between sampling points with non-identical genotypes.
#'
#' @param genotypes Vector of genotype labels
#' @param positions Vector of positions along a linear transect. Should be the
#' same length as `genotypes`.
#'
#' @return Mean distances between sampling points with identical and
#' non-identical genotype.
#'
#' @author Tom Ellis
#' @export
transect_clustering <- function(genotypes, positions){
  if(length(genotypes) != length(positions)){
    stop("Length of vectors for genotypes do not match.")
  }
  # Get the genotypes of all unqiue pairs of sampling points
  geno_pairs <- combn(genotypes, 2)
  # Vector of TRUE/FALSE for whether genotypes in each pair are the same
  geno_ix <- geno_pairs[1,] == geno_pairs[2,]

  # Distances between all unqiue pairs of sampling points.
  dist_pairs <- combn(positions, 2)
  dist_pairs <- abs(dist_pairs[1,] - dist_pairs[2,])

  # Return the distances between identical and non-identical genotypes
  c(
    identical = mean(dist_pairs[ geno_ix], na.rm = T),
    different = mean(dist_pairs[!geno_ix], na.rm = T)
  )
}
