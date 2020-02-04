#' Simulate replicate populations
#'
#'`simmiad` Simulates multiple replicate populations under a set of input
#'parameters, and outputs a summary of the clustering among genotypes.
#'
#'`simmiad` is a wrapper function to call multiple instances of `sim_population`,
#'estimate the degree of clustering in each, and save results to disk. See the
#'help file `?sim_population` for details of individual population simulations.
#'
#'For each simulation, `simmiad` draws a random horizontal transect and saves
#'two kinds of summary information about this transect to disk. First,
#'in the final generation of each population, a single transect is taken through
#'a random row of the grid, and the degree of clustering estimated using
#'`transect_cluster`. This summarises the mean distance between pairs of
#'identical genotypes and between pairs of non-identical genotypes.
#'Second, it estimates the spatial stability of the populations by calculating
#'how often a sampling point is occupied by idnetical genotypes at different
#'time points. This is done for the equivalent of pairs of successive real
#'sampling years (1984 to 1988, 1988 to 1992 etc).
#'
#'@inheritParams sim_population
#'@param nsims Int >0. Number of replicate populations to simulate.
#'@param transect_length Int between 0 and grid_size. Number of sampling points
#'to draw a transect for. Defaults to `grid_size`.
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
  grid_size,
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes = 124,
  transect_length = NULL,
  filename='simmiad',
  nsims){

  t0 <- proc.time()[3] # record the starting time.
  # Filenames to save output and logs
  logfile <- paste(filename, "log", sep=".")
  outfile = paste(filename, "simmiad", sep=".")
  # Print message about sims
  cat("\nSimulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  cat("Output will be save to",outfile,"and progress recorded in",logfile,"\n")
  # Write some messages about simulation properties to a log file.
  cat("Log file for simulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), ".\n", file = logfile)
  cat("See ", filename,".csv for simulation output.\n", sep="", file = logfile, append = TRUE)
  cat("\nSimulating a grid of",grid_size,"by",grid_size,"individuals of", n_starting_genotypes,"intial genotypes for",n_generations,"generations.\n", file = logfile, append = TRUE)
  cat("Seed diserpsal distances are exponentially distributed with mean =", mean_dispersal_distance,"metres.\n", file = logfile, append = TRUE)
  cat("Plants outcross with probability ",outcrossing_rate,".\n", sep = "", file = logfile, append=TRUE)

  cat("\nBeginning",nsims,"simulations...\n", file = logfile, append=TRUE)
  # Empty file to store output
  cnames <- c("i", "identical", "different", "84-88", "88-92", "92-96", "96-02", "02-14","14-16","16-18")
  write.table(matrix(cnames, nrow=1),
              file = outfile,
              sep =",",
              col.names = FALSE,
              row.names = FALSE)

  pb <- txtProgressBar(min = 0, max = nsims, style = 3)
  for(i in 1:nsims){
    setTxtProgressBar(pb, i)
    # Simulate a single population.
    sm <- sim_population(
      grid_size = grid_size,
      n_starting_genotypes = n_starting_genotypes,
      n_generations = n_generations,
      mean_dispersal_distance = mean_dispersal_distance,
      outcrossing_rate = outcrossing_rate,
      verbose = F,
      return_all = T
    )

    # Get positions for a random transect.
    if(is.null(transect_length)) transect_length <- grid_size
    if(transect_length > grid_size){
      stop("Transect length must be smaller than grid size.")
    }
    # Choose a row.
    tx <- sample(1:grid_size, 1)
    # Choose a transect of length `transect_length` along that row.
    transect_start <- sample(1:(grid_size - transect_length + 1),1)
    ty <- transect_start : (-1 + transect_start + transect_length)

    # Spatial clustering
    # Pull out that transect in the final generation.
    transect_geno <- sm[[nsims]][tx, ty]
    # Get the mean distances between identical and non-identical genotypes.
    spatial <- transect_clustering(transect_geno, 1:transect_length)

    # Temporal stability
    # Positions in generations corresponding to 1984 to 2018 (zero is 2018).
    gx <- length(sm)- c(34, 30, 26,22,16, 4,2,0)
    # Run transect_stability on each pair of sampling years.
    temporal <- sapply(2:8,
           function(g){
             transect_stability(
               sm[[gx[g-1]]],
               sm[[gx[g  ]]]
             )
           }
    )

    # Create a row of output data.
    output <- matrix(c(i=i, spatial, temporal), nrow = 1)
    output <- round(output, 3)

    # Write to disk.
    write.table(output,
              file = outfile,
              sep =",",
              append    = TRUE,
              col.names = FALSE,
              row.names = FALSE)
  }
  close(pb)

  t1 <- proc.time()[3] # record the starting time.
  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes\n",
      file = logfile, append=T)
  cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes.\n\n\n")
}



