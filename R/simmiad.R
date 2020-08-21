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
#'`transect_cluster`. This summarises the mean distance between pairs of
#'identical genotypes and between pairs of non-identical genotypes.
#'Second, it estimates the spatial stability of the populations by calculating
#'how often a sampling point is occupied by idnetical genotypes at different
#'time points. This is done by comparing the final generation to each of the
#'twelve previous generations.
#'
#'@inheritParams sim_population
#'@param nsims Int >0. Number of replicate populations to simulate.
#'@param filename Str. Directory where output and log file should be saved.
#'
#'@return Nothing will be printed on screeen, but two files are saved to disk:
#'1. A CSV file giving genotype distances between identical and
#'non-identical plants, and the output of transect_stability for each pair of
#'sampling years.
#'2. A log file giving simulation details.
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
  filename
){
  t0 <- proc.time()[3] # record the starting time.
  if(n_generations < 12) {
    stop("n_generations must be at least 12.")
  }
  # Filenames to save output and logs
  logfile <- paste(filename, "log", sep=".")
  outfile = paste(filename, "simmiad", sep=".")
  # Print message about sims
  cat("\nSimulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  cat("Output will be save to",outfile,"and progress recorded in",logfile,"\n")
  # Write some messages about simulation parameters to a log file.
  cat("Log file for simulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), ".\n", file = logfile)
  cat("See ", filename,".csv for simulation output.\n\n", sep="", file = logfile, append = TRUE)

  cat("Seed diserpsal distances are exponentially distributed with mean =", mean_dispersal_distance,"metres.\n", file = logfile, append = TRUE)
  cat("Plants outcross with probability ",outcrossing_rate,".\n", sep = "", file = logfile, append=TRUE)
  cat("\nBeginning",nsims,"simulations...\n", file = logfile, append=TRUE)

  # Empty list to store data
  output <- vector(mode = "list", length = nreps)

  pb <- txtProgressBar(min = 2, max = nsims, style = 3)
  for(i in 2:nsims){
    setTxtProgressBar(pb, i)

    # Simulate a single replicate population
    sm <-sim_population(
      population_size = population_size,
      mean_dispersal_distance = mean_dispersal_distance,
      outcrossing_rate = outcrossing_rate,
      n_generations = n_generations,
      n_starting_genotypes = n_starting_genotypes,
      density = density,
      n_sample_points = n_sample_points,
      sample_spacing = sample_spacing
    )
    # send the data to output
    output[[i]] <- c(
      i = i,
      transect_clustering(sm[[n_generations]], (1:n_sample_points) * sample_spacing),
      t1  = transect_stability(sm[[n_generations-1]],  sm[[n_generations]]),
      t2  = transect_stability(sm[[n_generations-2]],  sm[[n_generations]]),
      t4  = transect_stability(sm[[n_generations-4]],  sm[[n_generations]]),
      t6  = transect_stability(sm[[n_generations-6]],  sm[[n_generations]]),
      t12 = transect_stability(sm[[n_generations-12]], sm[[n_generations]])
    )
  }
  close(pb)

  # Concatenate the list to a data.frame
  output <- do.call('rbind', output)
  output <- as.data.frame(output)

  colnames(output) <- c("i", "identical", "different", "t1", "t2", "t4", "t6", "t12")
  output$i <- as.integer(output$i)

  write.csv(output, outfile, row.names = F)

  t1 <- proc.time()[3] # record the end time.

  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes\n",
      file = logfile, append=T)
  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes.\n\n\n")
}
