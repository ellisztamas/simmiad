#' Simulate replicate populations
#'
#'`simmiad` Simulates multiple replicate populations under a set of input
#'parameters, and outputs a summary of the clustering among genotypes.
#'
#'`simmiad` is a wrapper function to call multiple instances of `sim_population`,
#'estimate the degree of clustering in each, and save results to disk. See the
#'help file `?sim_population` for details of individual population simulations.
#'
#'For each simulation, `simmiad` draws a horizontal transect at each generation
#'across the whole time series. For a subset of those generations it calculates
#'three kinds of summary information about this transect to disk for a subset.
#'First it estimates the degree of clustering estimated using
#'`transect_cluster`, and averages over points between years.
#'This summarises the number and mean distance between pairs of identical
#'genotypes and between pairs of non-identical genotypes.
#'Second, it estimates the temporal stability of the populations by calculating
#'how often a sampling point is occupied by identical genotypes at different
#'at the start of the samplinh period, and at subsequent time points.
#'Third, it calculates the probability that *adjacent* sampling points only are
#'occupied by identical genotypes for each year separately.
#'
#'@inheritParams sim_population
#'@param nsims Int >0. Number of replicate populations to simulate.
#'@param progress If TRUE, a progress bar is printed.
#'@param years_to_sample Vector of integers indexing which generations to sample
#'the transect to calculate spatial structure and temporal stability. Defaults
#'to the last 36 generations of the simulation
#'
#'@return A list of seven dataframes.
#'
#' 1. **parameters** A data.frame giving input parameters.
#' 2. **clustering** The covariance between distance along the transect and the
#' frequency of identical genotypes.
#' 3. **matching_pairs**: The number of pairs of identical genotypes in the
#' transect.
#' 4. **count_NA**: The number of empty sampling points.
#' 5. **n_genotypes**: The number of unique genotypes sampled in the transect
#' (note that this will be different from what you gave as `distance_identity`,
#' because the latter reflects genotypes in *the whole population*, not just in
#'  the transect).
#' 6. **stability**: How often individual sampling points are occupied by the
#' same genotype in the final generations and 1, 2, ..., n generations back.
#' 7. **distance_identity**: Probabilities of finding identical genotypes in
#' pairs of sampling points at all possible distances between transects. For
#' example, if there are five evenly spaced sampling points as in the example
#' above, there are four possible distances between sampling points. Rows
#' indicate replicate simulations.
#' 8. **di_by_year** The probability of finding identical genotypes in pairs
#' of sampling points across generations. Only values for adjacent sampling
#' points are shown.
#'
#' For 2-6 and 8, rows indicate replicate simulations and columns generations.
#'
#'@author Tom Ellis
#'@seealso `sim_population`, `transect_cluster`, `transect_stability`
#'@export
#'
simmiad <- function(
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes = 50,
  density = 3,
  n_sample_points = 30,
  sample_spacing = 5,
  range_limit = 1.5,
  nsims,
  progress = TRUE,
  dormancy,
  years_to_sample = 1 : n_generations,
  pop_structure = "uniform",
  mixing = mean_dispersal_distance,
  habitat_labels = NULL,
  selection_gradient = 0
){
  t0 <- proc.time()[3] # record the starting time.

  if(
    any(years_to_sample > n_generations) | any(years_to_sample < 1)
    ){
    print(years_to_sample)
    stop(strwrap(
      "`years_to_sample` should be a vector of integers indexing over which
      generations spatial stability and atemporal stability should be calculated.
      One or more values are beyond the range of the number of generations
      specified by `n_generations`.")
    )
  }
  if(!is.integer(years_to_sample)){
    stop("years_to_sample should be a vector of integers")
  }
  if(any(table(years_to_sample) > 1)){
    stop("There are replicate entries in years_to_sample")
  }

  # Print message about sims
  cat(
    "\nSimulations of wild Emmer wheat begun on",
    format(Sys.time(), "%a %b %d %X %Y"),
    "using simmiad", as.character(packageVersion('simmiad')),
    "\n"
  )

  # Empty structures to store data
  clustering     <- matrix(NA, nrow = nsims, ncol = n_generations)
  matching_pairs <- matrix(NA, nrow = nsims, ncol = n_generations)
  count_NAs      <- matrix(NA, nrow = nsims, ncol = n_generations)
  n_genotypes    <- matrix(NA, nrow = nsims, ncol = n_generations)
  stability      <- matrix(NA, nrow = nsims, ncol = length(years_to_sample))
  distance_identity <- matrix(NA, nrow=nsims, ncol=n_sample_points-1)
  di_by_year     <-matrix(NA, nrow = nsims, ncol = n_generations)
  if( !is.null(habitat_labels)){
    clustering_by_habitat <- matrix(NA, nrow = nsims, ncol = n_generations)
  } else clustering_by_habitat <- NULL

  # Simulate generations one at a time
  if(progress) pb <- txtProgressBar(min = 2, max = nsims, style = 3)
  for(i in 1:nsims){
    if(progress) setTxtProgressBar(pb, i)
    # Simulate a single replicate population
    sm <- sim_population(
      mean_dispersal_distance = mean_dispersal_distance,
      outcrossing_rate = outcrossing_rate,
      n_generations = n_generations,
      n_starting_genotypes = n_starting_genotypes,
      density = density,
      range_limit = range_limit,
      n_sample_points = n_sample_points,
      sample_spacing = sample_spacing,
      dormancy = dormancy,
      pop_structure = pop_structure,
      mixing = mixing,
      habitat_labels = habitat_labels,
      selection_gradient = selection_gradient
    )

    # Spatial clustering
    sample_positions <- (1:n_sample_points) * sample_spacing
    spatial <- sapply(
      sm,
      function(x) {
        transect_clustering(
          genotypes = x,
          positions = sample_positions,
          deme_size = max(sample_positions)/5
          )
      })
    # Covariance between distance and identity
    clustering[i,] <- spatial['covar',]
    # Number of matching pairs in the transect
    matching_pairs[i,] <- spatial['n_matches',]
    # Number of NA samples in each year
    count_NAs[i,] <- colSums(sapply(sm, is.na))
    # Number of unique genotypes in each year
    n_genotypes[i,] <- sapply(sm, function(x) length(unique(x)))
    # Probabilities of finding identical genotypes in pairs of sampling points
    # at different distances, averaged over years
    di <- distance_identities(
      genotypes = sm[years_to_sample],
      positions = sample_positions
    )
    di <- split(di, di$distances)
    distance_identity[i,] <- sapply(di, function(x){
      sum(x$matches * x$n, na.rm = T) / sum(x$n)
    })
    # Calculate the probability that sampling points within the closest distance
    # classes are occupied by the same DGG.
    di_all_years <- lapply(sm, distance_identities, sample_positions) # distance identity for each year separately
    # Get the average for the smallest distance class
    di_by_year[i,] <- sapply(di_all_years, function(x) {
      j <- x[x$distances == min(x$distances), ]
      sum(j$matches * j$n, na.rm = T) / sum(j$n)
      })
    # Calculate the probability that plants in the same habitat are identical.
    if("habitat_labels" %in% names(attributes(sm))){
      clustering_by_habitat[i,] <- sapply(sm, habitat_clustering, attributes(sm)$habitat_labels)
    }

    # Temporal stability
    temporal <- rep(NA, length(years_to_sample))
    for (g in years_to_sample){
      temporal[g] <- transect_stability(
        x = sm[[1]],
        y = sm[[g]]
      )
    }
    stability[i,] <- temporal

  }
  if(progress) close(pb)

  t1 <- proc.time()[3] # record the end time.
  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes.\n\n\n")

  # Create a table of parameters to export later.
  params <- parameter_table(
    mean_dispersal_distance = mean_dispersal_distance,
    outcrossing_rate = outcrossing_rate,
    n_generations = n_generations,
    n_starting_genotypes = ifelse(pop_structure=="hardcoded", length(n_starting_genotypes), n_starting_genotypes),
    nsims = nsims,
    density = density,
    dormancy = dormancy,
    n_sample_points = n_sample_points,
    sample_spacing = sample_spacing,
    range_limit = range_limit,
    years_to_sample = length(years_to_sample)
  )

  output <- list(
    parameters = params,
    clustering = clustering,
    matching_pairs = matching_pairs,
    count_NAs = count_NAs,
    n_genotypes = n_genotypes,
    stability = stability,
    distance_identity = distance_identity,
    di_by_year = di_by_year
  )
  if( !is.null(clustering_by_habitat) ){
    output$clustering_by_habitat = clustering_by_habitat
  }

  return(output)
}
