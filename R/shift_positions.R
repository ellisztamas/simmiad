#' Shift a matrix of x- and y-coordinates in random directions
#'
#' `shift_positions` peturbs the positions of a matrix of positions by
#' distances drawn from an exponential distribution. Dispersal is simulated in a
#' square habitat of width `box_limit` centred on x=y=0.
#'
#' Dispersal is in a box.
#' If the simulated dispersal distance takes a seed beyond the boundary of the
#' habitat it will be reflected back into the habitat.
#' Note that if dispersal distances are much larger than the habitat size then
#' it could be that the reflected distance goes beyond the *other* limit of the
#' population.
#'
#' @param coords Matrix with a row for each individual and two rows.
#' @param mean_dispersal_distance Float giving the average perturbation
#' distance. The reciprocal of this value is used as a rate for the exponential
#' distribution.
#' @param box_limit Float >1 giving half the width of the habitat.
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

  # Constrain dispersal to be within the habitat.
  # If distances are beyond the limit of the habitat, reflect them back the same distance.
  ix <- abs(new_coords) > box_limit
  # Plus/negative box_limit, depending on which direction the coordinate goes
  boundary <- sign(new_coords[ix]) * box_limit
  # Deviation of observed coordinates from the boundary
  dev <- new_coords[ix] - boundary
  # Reflect coordinates back across the boundary.
  new_coords[ix] <- boundary - dev

  return(new_coords)
}
