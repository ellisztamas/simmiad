#' Take a transect through a Simmiad population
#'
#' Take a transect of equally spaced sampling points through a population
#' of plants with spatial coordinates, and pick the plant closest to each point.
#'
#' #' Plants exist within a box whose size is determined by the population size and
#' plant density (more plants at lower density = larger box). `sim_population`
#' will throw and error if the length of the transect
#' (`(n_sample_points-1) * sample_spacing`) is longer than the width of the box.
#' In this case, you can either make the population larger and less dense, or
#' choose fewer, closer sampling points.
#'
#' @param coords Matrix with a row for each individual and two rows.
#' @param n_sample_points Number of points to sample along the transect
#' @param sample_spacing Distance between sampling points
#'
#' @return A vector of integers giving the row positions of plants to be sampled.
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
  # Get the column index of the plants closest to each sampling point.
  ix <- apply(distances_to_points, 2, which.min)

  return(ix)
}
