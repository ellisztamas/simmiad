library(mvtnorm)

#' Initialise population
#'
#' This creates a population of individuals as a 2D matrix of coordinates in a
#' torus and a list of genotypes.
#'
#'
#'
#' @param mean_dispersal_distance Float >0. Mean seed dispersal distance. The
#' reciprocal of this is used as the rate parameter to draw from the exponential
#' distribution if `method="clusters`.
#' @param n_starting_genotypes Int >0. Number of initial genotypes to start with.
#' Defaults to 50
#' @param box_limit Float >1 giving half the circumference of the torus.
#' `sim_population` based on the length of the transect and plant density.
#' @param population_size Integer > 1.
#' @param method Character string indicating the initial population structure to
#' be simulated. Passing "uniform" simulated a panmictic population. "clusters"
#' simulates clusters of identical individuals that disperse from distinct
#' mothers via exponential dispersal set by `mean_dispersal_distance`. This is
#' likely to generate very disperate clumps. Passing "mvnorm" simulates
#' uniformly distributed coordinates for indiduals, as well as centroid
#' positions for genotypes. Individuals are assigned a genotype in proportion to
#' their distance to each genotype centroid based on multivariate-normal
#' probabilities. The variance covariance matrix for this is set as
#' \code{sqrt(habitat_size/n_starting_genotypes) / 3} such that the tail of each
#' genotype just about touches those of its neighbours, on average.
#'
#' @return A list with two elements: `geno`, a vector of genotype labels;
#' `coords`, a 2D matrix of coordinate positions.
#'
#' @export
initialise_population <- function(
  mean_dispersal_distance = NULL,
  n_starting_genotypes,
  population_size,
  box_limit,
  method = "uniform"
) {
  stopifnot(mean_dispersal_distance > 0 )
  stopifnot(n_starting_genotypes > 0 )
  stopifnot(population_size > 0 )
  stopifnot(box_limit >= 1 )

  if(method == "uniform"){
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

  } else if( method == "clusters" ){
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
      mean_dispersal_distance = mean_dispersal_distance,
      box_limit = box_limit
    )

    return(pop)

  } else if (method == "mvnorm") {
    #' Simulate coordinates uniformly, but assign each to different genotypes
    #' based on distances to genotype-mean positions

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
    var_cov_diagonal <- sqrt(habitat_size/n_starting_genotypes) / 3 # each genotype overlaps its neighbour by one SD on average
    sigma <- matrix(c(
      var_cov_diagonal,0,
      0,var_cov_diagonal),
      ncol=2)
    # Multivariate-normal probabilities that each individual "belongs" to each genotype
    # Returns a matrix with a row for each individual and a column for each genotype
    # Columns sum to one.
    prob_each_genotype <- apply(geno_coords, 1, function(x){
      probs <- dmvnorm(
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

  } else {
    stop("`method` should be one of 'uniform', 'clusters' or 'mvnorm'.")
  }
}
