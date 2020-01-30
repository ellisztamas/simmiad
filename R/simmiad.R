#' Simulate replicate populations
#'
#'`simmiad` Simulates multiple replicate populations under a set of input
#'parameters, and outputs a summary of the clustering among genotypes.
#'
#'`simmiad` is a wrapper function to call multiple instances of `sim_population`,
#'estimate the degree of clustering in each, and save results to disk. See the
#'help file `?sim_population` for details of individual population simulations.
#'
#'In the final generation of each population, a single transect is taken through
#'a random row of the grid, and the degree of clustering estimated using
#'`transect_cluster`. This summarises the mean distance between pairs of
#'identical genotypes and between pairs of non-identical genotypes. Summary data
#'are saved to disk in a CSV file.
#'
#'Need to add sth about parallelising.
#'
#'@inheritParams sim_population
#'@param nsims Int >0. Number of replicate populations to simulate.
#'@param filename Str. Directory where output and log file should be saved.
#'
#'@return Nothing will be printed on screeen, but two files are saved to disk:
#'A CSV file giving genotype distances between identical and
#'non-identical plants, and a log file giving simulation details.
#'
#'@author Tom Ellis
#'@seealso `sim_population`, `transect_cluster`
#'@export
#'
simmiad <- function(
  grid_size,
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes = 124,
  filename='simmiad',
  nsims){

  t0 <- proc.time()[3] # record the starting time.
  # Filenames to save output and logs
  logfile <- paste(filename, "log", sep=".")
  outfile = paste(filename, "simmiad", sep=".")
  # Print message about sims
  cat("Simulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  cat("Output will be save to",outfile,"and progress recorded in",logfile,"\n")
  # Write some messages about simulation properties to a log file.
  cat("Log file for simulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), ".\n", file = logfile)
  cat("See ", filename,".csv for simulation output.\n", sep="", file = logfile, append = TRUE)
  cat("\nSimulating a grid of",grid_size,"by",grid_size,"individuals of", n_starting_genotypes,"intial genotypes for",n_generations,"generations.\n", file = logfile, append = TRUE)
  cat("Seed diserpsal distances are exponentially distributed with mean =", mean_dispersal_distance,"metres.\n", file = logfile, append = TRUE)
  cat("Plants outcross with probability ",outcrossing_rate,".\n", sep = "", file = logfile, append=TRUE)

  cat("\nBeginning",nsims,"simulations...\n", file = logfile, append=TRUE)
  # Empty file to store output
  write.table(matrix(c("i", "identical", "different"), nrow=1),
              file = outfile,
              sep =",",
              col.names = FALSE,
              row.names = FALSE)

  for(i in 1:nsims){
    # Simulate a single population.
    sm <- sim_population(
      grid_size = grid_size,
      n_starting_genotypes = n_starting_genotypes,
      n_generations = n_generations,
      mean_dispersal_distance = mean_dispersal_distance,
      outcrossing_rate = outcrossing_rate,
      verbose = F,
      return_all = F
    )

    # Pick a row to use as a transect.
    transect_row <- sample(1:grid_size, 1)
    # Pull out that transect.
    transect_geno <- sm[transect_row,]

    # Get the mean distances between identical and non-identical genotypes.
    output <- matrix(c(i=i, transect_clustering(transect_geno, 1:grid_size)), nrow = 1)
    output <- round(output, 3)

    # Write to disk.
    write.table(output,
              file = outfile,
              sep =",",
              append    = TRUE,
              col.names = FALSE,
              row.names = FALSE)
  }

  t1 <- proc.time()[3] # record the starting time.
  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes\n",
      file = logfile, append=T)
}



