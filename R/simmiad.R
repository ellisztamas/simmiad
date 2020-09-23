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
#'@param how_far_back Integer number of generations
#'
#'@return A data.frame giving simulation replicate, number of identical and
#'non-identical genotypes in the transect, mean distances between identical
#'and non-identical plants and the output of transect_stability for each pair
#'of sampling years.
#'
#'@author Tom Ellis
#'@seealso `sim_population`, `transect_cluster`, `transect_stability`
#'@export
#'
simmiad <- function(
  population_size,
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes = 50,
  density = 3,
  n_sample_points = 30,
  sample_spacing = 5,
  nsims,
  progress = TRUE,
  how_far_back = n_generations
){
  t0 <- proc.time()[3] # record the starting time.
  if(n_generations < 13) {
    stop("n_generations must be at least 12.")
  }
  if(how_far_back > n_generations){
    how_far_back <- n_generations
    warning(strwrap(
      "The number of generations over which to calculate temporal stability
      given by `how_far_back`) is greater than the total number of generations.
      This will be set to the maximum number of generations.")
      )
  }
  # Print message about sims
  cat("\nSimulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), "\n")

    # Empty list to store data
  output <- vector(mode = "list", length = nsims)

  if(progress) pb <- txtProgressBar(min = 2, max = nsims, style = 3)
  for(i in 1:nsims){
    if(progress) setTxtProgressBar(pb, i)

    # Simulate a single replicate population
    sm <-sim_population(
      mean_dispersal_distance = mean_dispersal_distance,
      outcrossing_rate = outcrossing_rate,
      n_generations = n_generations,
      n_starting_genotypes = n_starting_genotypes,
      density = density,
      n_sample_points = n_sample_points,
      sample_spacing = sample_spacing
    )

    # Spatial clustering
    sample_positions <- (1:n_sample_points) * sample_spacing
    spatial <- transect_clustering(sm[[n_generations]], sample_positions)

    # Temporal stability
    temporal <- numeric(how_far_back)
    for (g in 1:(how_far_back-1)){
      temporal[[g]] <- transect_stability(
        x = sm[[n_generations]],
        y = sm[[n_generations - g]]
      )
    }

    # send the data to output
    output[[i]] <- c(
      i = i,
      spatial,
      t1  = transect_stability(sm[[n_generations-1]],  sm[[n_generations]]),
      t2  = transect_stability(sm[[n_generations-2]],  sm[[n_generations]]),
      t4  = transect_stability(sm[[n_generations-4]],  sm[[n_generations]]),
      t6  = transect_stability(sm[[n_generations-6]],  sm[[n_generations]]),
      t12 = transect_stability(sm[[n_generations-12]], sm[[n_generations]])
    )
  }
  if(progress) close(pb)

  # Concatenate the list to a data.frame
  output <- do.call('rbind', output)
  output <- as.data.frame(output)
  colnames(output) <- c("i","n_matches", "n_diff", "d_matches", "d_diff", "t1", "t2", "t4", "t6", "t12")
  output$i <- as.integer(output$i)
  output$n <- as.integer(output$n)

  t1 <- proc.time()[3] # record the end time.
  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes.\n\n\n")

  return(output)
}
