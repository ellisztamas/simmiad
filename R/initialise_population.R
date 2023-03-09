library(mvtnorm)

#' Initialise population
#'
#' This creates a population of individuals as a 2D matrix of coordinates in a
#' torus and a list of genotypes.
#'
#' @param n_starting_genotypes Int >0. Number of initial genotypes to start with.
#' Defaults to 50
#' @param box_limit Float >1 giving half the circumference of the torus. Ignored
#' if \code{pop_structure} is a vector.
#' @param population_size Integer > 1. Calculated by \code{sim_population} based
#' on the length of the transect and plant density.
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
#'
#' @param mixing Float >0. Parameter controlling the degree of spatial
#' clustering of genotypes. Smaller values indicate more structure populations.
#' If \code{pop_structure='mvnorm'} this is a scaler multiplier for the variance
#' of the multivariate normal probability density.
#' If \code{pop_structure="clusters`} or \code{pop_structure='hardcoded'} this is
#' the reciprocal of the rate parameter to draw dispersal distances from the
#' exponential distribution.
#' @param sample_spacing Positive integer giving the distance between sampling
#' points if \code{pop_structure='hardcoded'}
#'
#' @return A list with two elements: `geno`, a vector of genotype labels;
#' `coords`, a 2D matrix of coordinate positions.
#'
#' @export
#'
initialise_population <- function(
    n_starting_genotypes,
    population_size,
    box_limit,
    pop_structure = "uniform",
    mixing = NULL,
    sample_spacing = 5
) {
  stopifnot(
    population_size > 0,
    box_limit >= 1
  )
  if( !is.null(mixing) ) stopifnot(mixing >= 0)

  if(pop_structure == "uniform"){
    # Start with genotypes well mixed throughout the habitat
    pop <- list(
      # vector of genotype labels
      geno = sample(
        x = paste("g", 1:n_starting_genotypes, sep = ""),
        size = population_size,
        replace = T
      ),
      # Matrix of coordinates.
      coords = matrix(
        runif(n = population_size * 2, min = -box_limit, max = box_limit),
        ncol=2
      )
    )
    return(pop)

  } else if( pop_structure == "clusters" ){
    if( is.null(mixing) ){
      stop("If `pop_structure` is not 'uniform', please supply a value for `mixing`.")
    }
    # Initialise the population with genotypes in normally-distributed clumps
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
      mean_dispersal_distance = mixing,
      box_limit = box_limit
    )

    return(pop)

  } else if (pop_structure == "mvnorm") {

    # Simulate coordinates uniformly, but assign each to different genotypes
    # based on distances to genotype-mean positions
    if( is.null(mixing) ){
      stop("If `pop_structure` is not 'uniform', please supply a value for `mixing`.")
    }
    # Coordinates of individuals
    ind_coords = matrix(
      runif(n = population_size*2, min = -box_limit, max = box_limit),
      ncol=2
    )
    # coordinates of genotype centres
    geno_coords = matrix(
      runif(n = n_starting_genotypes * 2, min = -box_limit, max = box_limit),
      ncol=2
    )

    # Overall area of the torus
    habitat_size <- (2*box_limit)^2
    # Variance-covariance matrix for the MV normal
    var_cov_diagonal <- sqrt(habitat_size/n_starting_genotypes) * mixing
    sigma <- matrix(c(
      var_cov_diagonal,0,
      0,var_cov_diagonal),
      ncol=2)
    # Multivariate-normal probabilities that each individual "belongs" to each genotype
    # Returns a matrix with a row for each individual and a column for each genotype
    # Columns sum to one.
    prob_each_genotype <- apply(geno_coords, 1, function(x){
      probs <- mvtnorm::dmvnorm(
        x = ind_coords,
        mean = c(x[1], x[2]),
        sigma = sigma
      )
    })
    # Normalise so rows sum to 1.
    prob_each_genotype <- prob_each_genotype / rowSums(prob_each_genotype)
    # Choose a genotype for each plant based on those probabilities
    geno <- apply(prob_each_genotype, 1,
                  function(row) sample(1:n_starting_genotypes, size = 1, prob = row)
    )

    pop <- list(
      geno = paste0("g", geno),
      coords = ind_coords
    )
    return(pop)

  } else if( pop_structure == "hardcoded") {
    # Start with a vector genotypes giving a known structure, and simulate
    # vertical bands of the same genotype.
    if( length(n_starting_genotypes) == 1){
      stop("If `pop_structure='hardcoded'` provide a vector of genotypes via `n_starting_genotypes`.")
    }
    real_transect_length <- length(n_starting_genotypes)

    # Holder for population details.
    pop <- list(coords = NULL, geno =NULL)

    # Simulate Generation Zero.
    # Create a vector of genotypes by rotating the input vector by a random amount
    # then replicating genotypes vertically.
    rotate_vector <- function(x, n){
      c(tail(x, n), head(x, -n))
    }
    pop$geno <- rotate_vector(x = n_starting_genotypes, n = sample(1:real_transect_length, 1))
    pop$geno <- rep(pop$geno, real_transect_length)
    # Matrix of x and y positions for each plant.
    # Plants in generation zero are arranged on an evenly spaced grid.
    coords_1d <- (1:real_transect_length) * sample_spacing
    centred_coords <- coords_1d - mean(coords_1d) # centre around zero
    pop$coords <- as.matrix(cbind(
      rep(centred_coords, each = real_transect_length),
      rep(centred_coords, real_transect_length)),
      ncol = 2
    )

    # Simulate Generation One.
    # Draw a sample of individuals for generation 1.
    ix <- sort(
      sample(
        x = 1:length(pop$geno),
        size = population_size,
        replace = TRUE
      )
    )
    # Update the genotypes.
    pop <- list(
      geno   = as.character(pop$geno[ix]),
      coords = pop$coords[ix,]
    )
    # Disperse from the mother
    pop$coords <- shift_positions(
      pop$coords,
      mean_dispersal_distance = mixing,
      box_limit = box_limit
    )
    return(pop)

  } else {
    stop("`pop_structure` should be one of 'uniform', 'clusters' or 'mvnorm'.")
  }
}
