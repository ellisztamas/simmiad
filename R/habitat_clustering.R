#' Do genotypes cluster by microhabitat?
#'
#' Calculates how many pairs of plants of the same genotype are also in the
#' same microhabitat.
#'
#' This calculates Fst as
#'
#' \deqn{
#' \frac{f_0-\bar{f}}{1-\bar{f}}
#' }
#'
#' where \eqn{f_0} is the probability that two individuals in the same habitat
#' are identical and \eqn{f \bar{f}} is the probability that two individuals
#' from the total population are identical.
#'
#' @param genotype Vector of genotype labels.
#' @param habitat Vector of habitat labels of the same length as `genotype`.
#'
#' @return Proportion of pairs of individuals in the same habitat that are
#' identical.
#' @author Tom Ellis
#' @export
habitat_clustering <- function(genotype, habitat){
  if(length(genotype) != length(habitat)){
    stop("`genotype` and `habitat` are different lengths.")
  }
  # For all pairs of plants, check whether genotypes are the same
  genotype_match <- combn(genotype, 2)
  genotype_match <- genotype_match[1,] == genotype_match[2,]
  # For all pairs of plants, check whether habitats are the same
  habitat_match <- combn(habitat, 2)
  habitat_match <- habitat_match[1,] == habitat_match[2,]
  # F statistics
  f_zero <- mean( genotype_match[ habitat_match ], na.rm = TRUE)
  f_bar  <- mean( genotype_match,                 na.rm = TRUE)

  (f_zero - f_bar) / (1 - f_bar)
}
