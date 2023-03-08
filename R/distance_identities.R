#' Decay of identity along a transect
#'
#' Calculates the probability of observing identical genotypes at different
#' distances along a transect.
#'
#' `distance_identities` calculates how often pairs of identical genotypes are
#' observed for all possible distances between sampling points along a transect.
#' Identities are averaged over multiple years for a transect.
#'
#' @param genotypes Vector of genotypes occupying each sampling point.
#' Alternatively a list of vectors can be provided. Each element in the list
#' should be a 1-d vector with an element for every sampling point along the
#' transect, with each element the genotype observed at that sampliong point.
#' @param positions 1-d vector of positions along the transect. This should be
#' the same length of each vector in `genotypes`.
#' @param labels Optional vector of labels for each sampling point. Defaults to
#' integers.
#'
#' @return Returns a data.frame with a row for each pair of sampling points.
#' Columns give the labels for the points being compared, their distance,
#' whether genotypes match, and the number of non-NA observations for
#' each pair. If `genotypes` is a list this will return the average number of
#' matches.
#' Matches will return NaN if sampling points at some distance
#' did not contain a plant.
#'
#' @author Tom Ellis
#' @export
#'
#' @examples
#' positions <- sort(runif(16, -5,5))
#'
#' # Example with vector
#' genotypes <- rep(letters[1:4], each=4)
#' distance_identities(genotypes, positions, labels = LETTERS[1:16])
#'
#' # Example with list
#' genotypes <- list(
#'   rep(letters[1:4], each=4),
#'     rep(letters[1:4], 4)
#'     )
#' distance_identities(genotypes, positions, labels = LETTERS[1:16])
distance_identities <- function(genotypes, positions, labels = 1:length(positions)){
  if(length(positions) != length(labels)){
    stop('labels should be the same length as positions.')
  }
  # Vector of all pairwise distances
  distances <- combn(positions, 2)
  distances <- abs(distances[1,] - distances[2,])

  # Vector of labels for each pair
  labels <- combn(labels, 2)
  labels <- paste(labels[1,], "-", labels[2,], sep="")

  if( is.list(genotypes) ){
    if(all(length(positions) != sapply(genotypes, length))){
      stop("One or more elements in genotypes is not the same length as the vector of positions")
    }
    # Identify whether genotypes match for each pair of sampling points
    # First, do this as a list of logical values
    matches <- vector('list', length(genotypes))
    for(g in 1:length(genotypes)){
      is_match <- combn(genotypes[[g]], 2)
      is_match <- is_match[1,] == is_match[2,]
      matches[[g]] <- is_match
    }
    matches <- do.call('rbind', matches)
    # Number of finite observations
    n <- colSums(!is.na(matches))
    # Average number of matches over rows
    matches <- colMeans(matches, na.rm = T)
    # Tidy output and sort
    output <- data.frame(labels, distances, matches, n)
    output[order(output$distances),]

    return(output)

  }
  if( is.vector(genotypes) ){
    if ( length(genotypes) != length(positions) ) {
      stop("genotypes and positions are different lengths")
    }
    # Identify whether genotypes match for each pair of sampling points
    # First, do this as a list of logical values
    is_match <- combn(genotypes, 2)
    is_match <- is_match[1,] == is_match[2,]
    matches <- as.numeric(is_match)
    # Number of finite observations
    n <- as.numeric(!is.na(matches))
    # Tidy output and sort
    output <- data.frame(labels, distances, matches, n)
    return (output[order(output$distances),])
  }

}
