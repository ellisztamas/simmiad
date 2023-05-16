#' Do genotypes cluster by microhabitat?
#'
#' Calculates how many pairs of plants of the same genotype are also in the
#' same microhabitat.
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
  # Check whether genotype and habitat match.
  both_match <- genotype_match * habitat_match

  sum(both_match, na.rm = TRUE) / sum(habitat_match, na.rm = TRUE)
}
