#' Shift a matrix of x- and y-coordinates in random directions
#'
#' `shift_positions` peturbs the positions of a matrix of positions by
#' distances drawn from an exponential distribution. Dispersal is simulated on
#' a torus, so that if positions are peturbed beyond a pre-defined range limit
#' they are reflected to the other side of the population.
#'
#' @param coords Matrix with a row for each individual and two rows.
#' @param mean_dispersal_distance Float giving the average perturbation
#' distance. The reciprocal of this value is used as a rate for the exponential
#' distribution.
#' @param box_limit Float >1 giving half the circumference of the torus.
#'
#' @return A matrix matching the shape of the input `coords`, but with values
#' shifted by some exponential distance.
#' @export
shift_positions <- function(coords, mean_dispersal_distance, box_limit){
  if(!is.matrix(coords)){
    stop("`coords` should be a matrix.")
  }
  if(ncol(coords) != 2){
    stop("`coords` should have two columns." )
  }
  if(any(coords>box_limit)){
    stop("One or more values in `coords` is beyond the range limit.")
  }
  population_size <- nrow(coords)
  # Draw a total dispersal distance for each individual
  d <- rexp(
    n=population_size,
    rate = 1/mean_dispersal_distance
  )
  # The distance is split along the x- and y-axes.
  # Draw a vector of proporions to determine how much is distributed along each.
  prop <- runif(population_size, 0, 1)

  # Distances to peturb x- and y-axes.
  shift_positions <- matrix(
    c(sqrt(d^2 * prop),
      sqrt(d^2 * (1-prop))
    ),
    ncol = 2
  )
  # Randomly move x- and y-axes in the negative or positive directions
  shift_positions <- shift_positions * sample(c(1,-1), size = 2*population_size, replace = T)

  # Apply peturbations to the current positions
  new_coords <- coords + shift_positions

  # Ensure dispersal is on a torus.
  # When dispersal is beyond the edge of the population, reflect it back to the
  # other side of the range
  ix <- abs(new_coords) > box_limit
  dev <- new_coords[ix] - (sign(new_coords[ix]) * box_limit)
  new_coords[ix] <- -new_coords[ix] + 2*dev

  return(new_coords)
}
