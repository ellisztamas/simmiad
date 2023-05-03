#' Simulate replicate populations
#'
#'`simmiad` Simulates multiple replicate populations under a set of input
#'parameters, and outputs a summary of the clustering among genotypes.
#'
#'`simmiad` is a wrapper function to call multiple instances of `sim_population`,
#'estimate the degree of clustering in each, and save results to disk. See the
#'help file `?sim_population` for details of individual population simulations.
#'
#'For each simulation, `simmiad` draws a horizontal transect and saves
#'two kinds of summary information about this transect to disk. First,
#'in the final generation of each population, a single transect is taken through
#'a random row of the grid, and the degree of clustering estimated using
#'`transect_cluster`. This summarises the number and mean distance between pairs of
#'identical genotypes and between pairs of non-identical genotypes.
#'Second, it estimates the spatial stability of the populations by calculating
#'how often a sampling point is occupied by idnetical genotypes at different
#'time points. This is done by comparing the final generation to 1, 2, 4, 6 and
#'twelve years prior to that.
#'
#'@inheritParams sim_population
#'@param nsims Int >0. Number of replicate populations to simulate.
#'@param progress If TRUE, a progress bar is printed.
#'@param stability_years Integer number of generations to sample with a transect
#'@param clustering_years Vector of integers indexing which generations to average
#'over when calculating distance_identities(). Defaults to the last 50% of
#'generations.
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
  stability_years = n_generations,
  clustering_years = n_generations:(n_generations/2),
  pop_structure = "uniform",
  mixing = mean_dispersal_distance
){
  t0 <- proc.time()[3] # record the starting time.

  if(stability_years > n_generations){
    stability_years <- n_generations
    warning(strwrap(
      "The number of generations over which to calculate temporal stability
      given by `stability_years`) is greater than the total number of generations.
      This will be set to the maximum number of generations.")
    )
  }
  if(!is.integer(clustering_years)){
    stop("clustering_years should be a vector of integers")
  }
  if(any(table(clustering_years) > 1)){
    stop("There are replicate entries in clustering_years")
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
  stability      <- matrix(NA, nrow = nsims, ncol = stability_years)
  distance_identity <- matrix(NA, nrow=nsims, ncol=n_sample_points-1)
  di_by_year     <-matrix(NA, nrow = nsims, ncol = n_generations)

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
      mixing = mixing
    )

    # Spatial clustering
    sample_positions <- (1:n_sample_points) * sample_spacing
    spatial <- sapply(
      sm,
      function(x) {
        transect_clustering(
          genotypes = x,
          positions = sample_positions)
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
      genotypes = sm[clustering_years],
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

    # Temporal stability
    temporal <- rep(NA, stability_years)
    for (g in 1:(stability_years)){
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
    stability_years = stability_years
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

  return(output)
}
