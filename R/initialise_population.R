#' Initialise population
#'
#' Set up an initial population where each genotype starts from a single
#' spot and disperses from there.
#'
#' @param mean_dispersal_distance Float >0. Mean seed dispersal distance. The
#' reciprocal of this is used as the rate parameter to draw from the exponential
#' distribution.
#' @param n_starting_genotypes Int >0. Number of initial genotypes to start with.
#' Defaults to 50
#' @param box_limit Float >1 giving half the circumference of the torus.
#' `sim_population` based on the length of the transect and plant density.
#' @param population_size Integer > 1.
#' @export
initialise_population <- function(
  mean_dispersal_distance,
  n_starting_genotypes,
  population_size,
  box_limit
) {
  stopifnot(mean_dispersal_distance > 0 )
  stopifnot(n_starting_genotypes > 0 )
  stopifnot(population_size > 0 )
  stopifnot(box_limit >= 1 )
  # Initialise the population with randomly dispersed genotypes
  pop <- list(
    # One label for each genotype
    geno = paste("g", 1:n_starting_genotypes, sep=""),
    # Initialise positions for each plant
    coords = matrix(
      runif(n = n_starting_genotypes * 2, min = -box_limit, max = box_limit),
      ncol=2
    )
  )
  # Draw a sample of individuals for generation 1.
  ix <- sample(
    x = 1:n_starting_genotypes,
    size = population_size,
    replace = TRUE
  )
  # Update the genotypes.
  pop <- list(
    geno   = pop$geno[ix],
    coords = pop$coords[ix,]
  )
  # Disperse from the mother
  pop$coords <- shift_positions(
      pop$coords[ix,],
      mean_dispersal_distance = mean_dispersal_distance,
      box_limit = box_limit
  )

  pop
}
