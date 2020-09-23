#' Decay of identity along a transect
#'
#' Calculates the probability of observing identical genotypes at different
#' distances along a transect.
#'
#' `distance_identities` calculates how often pairs of identical genotypes are
#' observed for all possible distances between sampling points along a transect.
#' Identities are averaged over multiple years for a transect.
#'
#' @param transects List of transects. Each element in the list should be a 1-d
#' vector with an element for every sampling point along the transect, with each
#' element the genotype observed at that sampliong point.
#' @param positions 1-d vector of positions along the transect. This should be
#' the same length of each vector in `transects`.
#' @param labels Optional vector of labels for each sampling point. Defaults to
#' integers.
#'
#' @return Returns a data.frame with a row for each pair of sampling points.
#' Columns give the labels for the points being compared, their distance, the
#' number of matches between them, and the number of non-NA observations for
#' each pair.
#'
#' @author Tom Ellis
#' @export
#'
#' @examples
#' positions <- sort(runif(16, -5,5))
#' transects <- list(
#'   rep(letters[1:4], each=4),
#'     rep(letters[1:4], 4)
#'     )
#' distance_identities(transects, positions, labels = LETTERS[1:16])
distance_identities <- function(transects, positions, labels = 1:length(positions)){
  if(all(length(positions) != sapply(transects, length))){
    stop("One or more elements in transects is not the same length as the vector of positions")
  }
  if(length(positions) != length(labels)){
    stop('labels should be the same length as positions.')
  }
  # Vector of all pairwise distances
  distances <- combn(positions, 2)
  distances <- abs(distances[1,] - distances[2,])

  # Vector of labels for each pair
  labels <- combn(labels, 2)
  labels <- paste(labels[1,], "-", labels[2,], sep="")

  # Identify whether genotypes match for each pair of sampling points
  # First, do this as a list of logical values
  matches <- vector('list', length(transects))
  for(g in 1:length(transects)){
    is_match <- combn(transects[[g]], 2)
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

  output
}
