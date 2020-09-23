#' Take a transect through a Simmiad population
#'
#' Take a transect of equally spaced sampling points through a population
#' of plants with spatial coordinates, and pick the plant closest to each point.
#'
#' @param coords Matrix with a row for each individual and two rows.
#' @param n_sample_points Number of points to sample along the transect
#' @param sample_spacing Distance between sampling points
#'
#' @return A vector of integers giving the row positions of plants to be sampled.
#' If no plants are within 1m of the sampling point, returns NA.
#' @author Tom Ellis
take_transect <- function(coords, n_sample_points, sample_spacing){
  # Total length of the transect
  transect_length <- (n_sample_points-1) * sample_spacing

  # Define sampling points along y=0
  sampling_points <- seq(-transect_length/2, transect_length/2, sample_spacing)
  # Distances of all plants to each sampling point
  distances_to_points <- sqrt(
    outer(coords[,1], rep(0, n_sample_points), "-") ^ 2 +
      outer(coords[,2], sampling_points, '-') ^2
  )
  # Remove plants that are more than 1m away
  distances_to_points[distances_to_points > 1] <- NA

  # If all plants are >1m from a sampling point, return NA
  # Otherwise, get the column index of the plants closest to each sampling point.
  ix <- apply(distances_to_points, 2, function(x) {
    ifelse(all(is.na(x)), yes=NA, no = which.min(x) )
  })

  return(ix)
}
