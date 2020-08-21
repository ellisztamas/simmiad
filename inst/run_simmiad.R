#!/usr/bin/env Rscript
library("optparse")
library("simmiad")

option_list = list(
  make_option(c("-n", "--population_size"), type="integer",
              help="Number of individuls in the population.", metavar="character"),
  make_option(c("-d", "--dispersal"), type="double",
              help="Mean seed dispersal distance. The reciprocal of this is used as the rate parameter to draw from the exponential distribution.",
              metavar="character"),
  make_option(c("-o", "--outcrossing_rate"), type="double",
              help="Float between 0 and 1. Probability that an individual is outcrossed.",
              metavar="character"),
  make_option(c("-g", "--n_generations"), type="integer",
              help="Int >12. Number of generations to run the simulations.",
              metavar="character"),
  make_option(c("-b", "--n_starting_genotypes"), type="integer",
              help="Int >0. Number of initial genotypes to start with. Defaults to 50",
              metavar="character"),
  make_option(c("-a", "--density"), type="double",
              help="Float >0. Average density of plants per square metre.",
              metavar="character"),
  make_option(c("-p", "--n_sample_points"), type="integer",
              help="Number of points to sample along the transect.",
              metavar="character"),
  make_option(c("-s", "--sample_spacing"), type="double",
              help="Distance between sampling points.",
              metavar="character"),
  make_option(c("-r", "--nsims"), type="integer",
              help="Int >0. Number of replicate populations to simulate.",
              metavar="character"),
  make_option(c("-f", "--filename"), type="character",
              help="Str. Filename for the output data and log file, without a suffix.",
              metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Filenames to save output and logs
logfile <- paste(opt$filename, "log", sep=".")
outfile = paste(opt$filename, "simmiad", sep=".")

# Write some messages about simulation properties to a log file.
cat("Log file for simulations of wild Emmer wheat begun on", format(Sys.time(), "%a %b %d %X %Y"), ".\n",
    file = logfile)
cat("See ", opt$filename,".csv for simulation output.\n", sep="",
    file = logfile, append = TRUE)
cat("\nSimulating a population of", opt$population_size, "individuals of",
    opt$n_starting_genotypes,"intial genotypes for", opt$n_generations,"generations.\n",
    file = logfile, append = TRUE)
cat("Seed diserpsal distances are exponentially distributed with mean =", opt$dispersal,"metres.\n",
    file = logfile, append = TRUE)
cat("Plants outcross with probability ", opt$outcrossing_rate,".\n", sep = "",
    file = logfile, append=TRUE)
cat("\nBeginning", opt$nsims," replicate simulations...\n",
    file = logfile, append=TRUE)

t0 <- proc.time()[3] # record the starting time.

# Run replciate simulations
output <- simmiad(
  population_size = opt$population_size,
  mean_dispersal_distance = opt$dispersal,
  outcrossing_rate = opt$outcrossing_rate,
  n_generations = opt$n_generations,
  n_starting_genotypes = opt$n_starting_genotypes,
  density = opt$density,
  n_sample_points = opt$n_sample_points,
  sample_spacing = opt$sample_spacing,
  nsims = opt$nsims
)
#Save to disk
write.csv(output, opt$filename, row.names = F)

t1 <- proc.time()[3] # record the end time.
cat("\nSimulations completed", format(Sys.time(), "%a %b %d %X %Y"), "after", round((t1-t0)/60, 2), "minutes\n",
    file = logfile, append=T)
