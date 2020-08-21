#' Simulate a population
#'
#' Simulate a single population through time on a torus, and sample a transect
#' for each generation.
#'
#'`sim_population` simulates a population of annuals dispersing in continuous
#' space following exponentially distributed seed dispersal distances.
#'
#' 1. Individuals are drawn from a set of starting genotypes and randomly
#' distributed in a grid. Individuals are labelled by their genotype as 'g'
#' followed by an integer label.
#' 2. Individuals are chosen at random to found the next generation. They
#' move seed dispersal to new location, with distances drawn from an exponential
#' distribution. The population exists on a torus to eliminat edge effects;
#' if dispersal takes a plant of the edge of the range it moves to the other
#' side of the range.
#' 3. A subset of the new generation are chosen to have been germinated from
#' outcrossed seed. If so, they are given a new unique label by appending their
#' current genotype label by the generation number plus a unique integer.
#' 4. At each generation a transect through the population is drawn at x=0 with
#' sampling points at equally spaced intervals. The plant closest to each
#' sampling point is recorded.
#'
#' A seed bank is not currently implemented.
#'
#' Plants exist within a box whose size is determined by the population size and
#' plant density (more plants at lower density = larger box). `sim_population`
#' will throw and error if the length of the transect
#' (`(n_sample_points-1) * sample_spacing`) is longer than the width of the box.
#' In this case, you can either make the population larger and less dense, or
#' choose fewer, closer sampling points.
#'
#' @param population_size Int >0. Number of individuls in the population.
#' @param mean_dispersal_distance Float >0. Mean seed dispersal distance. The
#' reciprocal of this is used as the rate parameter to draw from the exponential
#' distribution.
#' @param outcrossing_rate Float between 0 and 1. Probability that an individual
#' is outcrossed.
#' @param n_generations Int >34. Number of generations to run the simulations.
#' @param n_starting_genotypes Int >0. Number of initial genotypes to start with.
#' Defaults to 50
#' @param density Float >0. Average density of plants per square metre.
#' @param n_sample_points Number of points to sample along the transect.
#' @param sample_spacing Distance between sampling points.
#'
#' @return A list of genotypes recorded at each sampling point in each
#' generation.
#'
#' @author Tom Ellis
#' @export
sim_population <- function(
  population_size,
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes = 50,
  density = 3,
  n_sample_points = 30,
  sample_spacing = 5
){
  # Plants exist in a box centred on zero.
  # Range limit is half the width of the box.
  range_limit <- sqrt(population_size / density) / 2
  # Empty list to store transects in each generation
  samples <- vector(mode = "list", length = n_generations)

  # Initialise the population with randomly dispersed genotypes
  geno <- sample(
    x = paste("g", 1:n_starting_genotypes, sep = ""), # vector of genotype labels
    size = population_size,
    replace = T
  )
  # Initialise positions for each plant
  coords <- matrix(
    runif(
      n = population_size * 2,
      min = -range_limit,
      max = range_limit
    ),
    ncol=2
  )
  # Take a transect through generation 1.
  tx <- take_transect(
    coords,
    n_sample_points = n_sample_points,
    sample_spacing = sample_spacing,
    range_limit = range_limit)
  samples[[1]] <- geno[tx]

  # Loop through subsequent generations
  for(g in 2:n_generations){
    # Sample plants to reproduce at random
    ix <- sample(1:population_size, replace = T)
    # Update the genotypes
    geno <- geno[ix]

    # Choose plants at random to receive outcrossed pollen.
    cross01 <- rbinom(n = population_size, 1, outcrossing_rate)
    cross01 <- as.logical(cross01)
    # Give these plants a new unique genotype by appending generation number and
    # and integer from 1 to the number of outcrossers.
    geno[cross01] <- paste(geno[cross01], "_", g-1, ".", 1:sum(cross01), sep = "")

    # Peturb positions
    coords <- shift_positions(
      coords[ix,],
      mean_dispersal_distance = mean_dispersal_distance,
      range_limit = range_limit
    )
    # Take a transect of the new generation
    tx <- take_transect(
      coords,
      n_sample_points = n_sample_points,
      sample_spacing = sample_spacing,
      range_limit = range_limit
    )
    samples[[g]] <- geno[tx]
  }
  return(samples)
}
