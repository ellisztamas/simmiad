#' Simulate a population
#'
#' Simulate a single population through time on a torus, and sample a transect
#' for each generation.
#'
#'`sim_population` simulates a population of annuals dispersing in continuous
#' space following exponentially distributed seed dispersal distances.
#'
#' 1. Individuals are drawn from a set of starting genotypes and randomly
#' distributed in a square habitat Individuals are labelled by their genotype
#' as 'g' followed by an integer label.
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
#' Plants exist within a box centred around zero with a transect going through
#' zero. The width of the box is defined to be 1.5-fold larger than the length
#' of the transect (but this can be changed with the argument `range_limit`).
#' Total population size is then the density of plants multiplied by the squared
#' width of the box.
#'
#' Genotypes are automatically assigned a phenotype drawn from a standard normal
#' distribution which can be used to simulate directional selection if
#' `selection_gradient` is more or less than zero. The fitness function is a
#' logistic function of minus the product of the phenotype and the selection
#' gradient, and therefore reflects selection on a something like survival.
#'
#' @param mean_dispersal_distance Float >0. Mean seed dispersal distance. The
#' reciprocal of this is used as the rate parameter to draw from the exponential
#' distribution.
#' @param outcrossing_rate Float between 0 and 1. Probability that an individual
#' is outcrossed.
#' @param n_generations Int >12. Number of generations to run the simulations.
#' @param n_starting_genotypes Int >0. Number of initial genotypes to start with.
#' Defaults to 50
#' @param density Float >0. Average density of plants per square metre.
#' @param n_sample_points Number of points to sample along the transect.
#' @param sample_spacing Distance between sampling points.
#' @param dormancy Float between 0 and 1. Probability that a seedling is drawn
#' from the seed bank. Seedlings are drawn from the prior generation with
#' probability 1-dormancy.
#' @param range_limit Float >1 defining how much wider than the transect the
#' width of the habitat should be. This, along with plant density, determines
#' how many plants will be simulated. Defaults to 1.5. Ignored if
#' \code{pop_structure="hardcoded"}
#'
#' @param pop_structure Character string indicating the initial population
#' structure to be simulated. Passing "uniform" simulated a panmictic population.
#' "clusters" simulates clusters of identical individuals that disperse from
#' distinct mothers via exponential dispersal set by \code{mean_dispersal_distance}.
#' This is likely to generate very disperate clumps. Passing "mvnorm" simulates
#' uniformly distributed coordinates for indiduals, as well as centroid
#' positions for genotypes. Individuals are assigned a genotype in proportion to
#' their distance to each genotype centroid based on multivariate-normal
#' probabilities. The variance covariance matrix for this is set as
#' \code{sqrt(habitat_size/n_starting_genotypes) / 3} such that the tail of each
#' genotype just about touches those of its neighbours, on average.
#' If \code{pop_structure='hardcoded'} and a vector of genotypes is passed to
#' \code{n_starting_genotypes}, for example
#' observed genotypes from along all real-world transects, this simulates
#' bands of identical genotypes by copying the vector over an evenly-
#' spaced grid (giving horizontal but not vertical structure). There is one
#' round of dispersal from this initial generation via exponential dispersal
#' (controlled by \code{mixing}) and to get the population to the correct
#' population density.
#' @param pop_structure Character string indicating the initial population
#' structure to be simulated. Passing "uniform" simulated a panmictic population.
#' 'clusters' simulates clusters of identical individuals that disperse from
#' distinct mothers via exponential dispersal set by
#' \code{mean_dispersal_distance}.
#' This is likely to generate very disperate clumps. Passing "mvnorm" simulates
#' uniformly distributed coordinates for indiduals, as well as centroid
#' positions for genotypes. Individuals are assigned a genotype in proportion to
#' their distance to each genotype centroid based on multivariate-normal
#' probabilities. The variance covariance matrix for this is set as
#' \code{sqrt(habitat_size/n_starting_genotypes) / 3} such that the tail of each
#' genotype just about touches those of its neighbours, on average.
#' If \code{pop_structure='hardcoded'} and a vector of genotypes is passed to
#' \code{n_starting_genotypes}, for example
#' observed genotypes from along all real-world transects, this simulates
#' bands of identical genotypes by copying the vector over an evenly-
#' spaced grid (giving horizontal but not vertical structure). There is one
#' round of dispersal from this initial generation via exponential dispersal
#' (controlled by \code{mixing}) and to get the population to the correct
#' population density.
#' @param mixing Float >0. Parameter controlling the degree of spatial
#' clustering of genotypes. Smaller values indicate more structure populations.
#' If \code{pop_structure='mvnorm'} this is a scaler multiplier for the variance
#' of the multivariate normal probability density.
#' If \code{pop_structure="clusters"} or \code{pop_structure='hardcoded'} this is
#' the reciprocal of the rate parameter to draw dispersal distances from the
#' exponential distribution.
#' #' @habitat_labels Optional vector of habitat labels when
#' `pop_structure = 'hard-coded`, with an element for each sample given in
#' `n_starting_genotypes`.
#'
#' @param selection_gradient Float giving strength of directional selection.
#' For no selection, give a zero (the default).
#'
#' @return A list of genotypes recorded at each sampling point in each
#' generation.
#'
#' @author Tom Ellis
#' @export
sim_population <- function(
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes,
  density,
  n_sample_points,
  sample_spacing,
  dormancy,
  range_limit = 1.5,
  pop_structure = 'uniform',
  mixing = mean_dispersal_distance,
  habitat_labels = NULL,
  selection_gradient = 0
){
  stopifnot(range_limit > 1)
  stopifnot(density > 0)
  stopifnot(dormancy >= 0 & dormancy <= 1)
  # Plants exist in a box centred on zero through which the transect runs
  # Range limit is half the width of the box.
  if(pop_structure == "hardcoded"){
    box_limit <- length(n_starting_genotypes) * sample_spacing /2
  } else {
    box_limit <- ((n_sample_points-1) * sample_spacing * range_limit)/2
  }
  # Given a density of plants per sq. metre and a size of the box, calculate
  # how many individuals you need
  population_size <- round(density * (2*box_limit)^2)

  # Empty list to store transects in each generation
  samples <- vector(mode = "list", length = n_generations)

  # Initialise the population with randomly dispersed genotypes
  pop <- initialise_population(
    n_starting_genotypes = n_starting_genotypes,
    population_size = population_size,
    box_limit = box_limit,
    pop_structure=pop_structure,
    mixing = mixing,
    habitat_labels = habitat_labels
  )
  # If they exist, save habitat_labels for later
  if("habitat_labels" %in% names(attributes(pop)) ){
    rotated_habitat_labels <- attributes(pop)$habitat_labels
    # Take the n_sampling_points in the centre of the vector of habitat labels
    # If n_sampling_points is odd, round down
    ix <- floor(
      (length(rotated_habitat_labels) - n_sample_points) / 2 # index where to start counting
    )
    trimmed_habitat_labels <- rotated_habitat_labels[(ix+1) : (ix+n_sample_points)]
  }

  # Initially, the seed bank and current generation are the same, but will change in the loop.
  seed_rain <- seed_bank <- pop
  # Take a transect through generation 1.
  tx <- take_transect(
    pop$coords,
    n_sample_points = n_sample_points,
    sample_spacing = sample_spacing
    )
  # Get genotypes along the transect by indexing pop$geno with tx
  # If tx contains only NAs this returns a huge vector of NAs, causing a crash
  # This hack gets around that.
  if(all(is.na(tx))) {
    samples[[1]] <- tx
  } else {
    samples[[1]] <- pop$geno[tx]
  }

  # Loop through subsequent generations
  for(g in 2:n_generations){
    # Sample plants to reproduce at random
    pop <- update_population(
      seed_rain = seed_rain,
      seed_bank = seed_bank,
      mean_dispersal_distance = mean_dispersal_distance,
      outcrossing_rate = outcrossing_rate,
      dormancy = dormancy,
      generation = g,
      box_limit = box_limit,
      selection_gradient = selection_gradient
    )
    # The previous generation now becomes a seed bank
    seed_bank <- seed_rain
    seed_rain <- pop
    # Take a transect of the new generation
    tx <- take_transect(
      pop$coords,
      n_sample_points = n_sample_points,
      sample_spacing = sample_spacing
    )
    if(all(is.na(tx))) {
      samples[[g]] <- tx
    }  else {
      samples[[g]] <- pop$geno[tx]
    }
  }
  # Include habitat labels, if given
  if( exists('trimmed_habitat_labels') ){
    attr(samples, 'habitat_labels') <- trimmed_habitat_labels
  }

  return(samples)
}
