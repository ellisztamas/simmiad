#' Simulate a population
#'
#' Simulate a single population through time.
#'
#' This sets up a population with \eqn{N} squared sampling points arranged at 1m
#' intervals on a grid.
#' Sampling points are then population by randomly assigning one of a set of
#' starting genotypes to each starting point (i.e. there is no initial spatial
#' structure).
#' In the next generation, the genotype at each sampling point is filled by seed
#' dispersal, assuming dispersal distances are exponentially distributed.
#' A single donor plant is chosen from the previous generation in proportion to
#' the probability of seeds travelling that distance under the exponential.
#' To account for outcrossing, each establishing seed is randomly assigned to be
#' outcrossed seed or not in proportion to a given outcrossing rate.
#' Outcrossed seeds are given a new unique genotype.
#'
#' Because the total number of sampling points increases as the square of grid
#' width, the number of individuals is the limiting factor on simulations.
#' To speed things up `sim_population` first draws dispersal distances for
#' individual sampling points for all generations at once, and then fills in
#' the genotypes in a second step.
#' This is valid if we assume no differences in fitness or dispersal between
#' genotypes or though time (at least on average).
#'
#' @param grid_size Int >0. Number of rows in the grid. Population size will be
#' the squre of this number.
#' @param mean_dispersal_distance Float >0. Mean seed dispersal distance.
#' Sampling points are arranged at intervals of 1m, so this should be relative
#' to that.
#' @param outcrossing_rate Float between 0 and 1. Probability that an individual
#' is outcrossed.
#' @param n_generations Int >34. Number of generations to run the simulations.
#' @param n_starting_genotypes Int >0. Number of initial genotypes to start with.
#' Defaults to 124.
#' @param verbose Logical. If TRUE, messages about simulation progress will be
#' printed. Defaults tp
#' @param return_all Logical. If TRUE, the entire population history will be
#' returned as a list. If FALSE, only the final generation is returned as a
#' vector.
#'
#' @return If `return_all` is TRUE, genotypes over the entire population history
#' will be returned as a list. If FALSE, only genotypes in the final generation
#' is returned as a vector.
#'
#' @author Tom Ellis
#' @export
sim_population <- function(
  grid_size,
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes = 124,
  verbose = TRUE,
  return_all = FALSE
  ){

  if(n_generations < 35){
    stop("Simulations must run for at least 35 generations.")
  }

  t0 <- proc.time()[3] # Record the time when the simulation started.
  if(verbose) cat("Initialising population.\n")

  # Initialise a set of starting genotypes.
  iggs <- paste("g", 1:n_starting_genotypes, sep = "")

  # We will keep track of the genotype at each sampling point using a list with
  # one element for each generation
  pop <- vector('list', n_generations)
  # Initialise the first generation with a random draw from the genotypes.
  pop[[1]] <- sample(iggs, size = grid_size^2, replace = T)
  # Data.frame of grid positions for each plant.
  coords <- data.frame(
    x   = rep(1:grid_size, grid_size),
    y   = rep(1:grid_size, each = grid_size)
  )

  # SIMULATE SEED DISPERSAL
  if(verbose) cat("Simulating seed dispersal events for each sampling point...\n")
  # Empty list to store where the seed is coming from in each generation.
  # There's an element for each sampling point.
  # Note that this will get coerced to a matrix later on.
  seed_disp_draws <- vector('list', grid_size^2)
  # For each sampling point draw a position from which seed will disperse in
  # each generation.
  if(verbose) pb <- txtProgressBar(min = 1, max = (grid_size^2), style = 3)
  for(i in 1:(grid_size^2)){
    if(verbose) setTxtProgressBar(pb, i)
    # Vector of distances between this plant and all other plants
    dist_vector <- sqrt((coords$x[i] - coords$x)^2 + (coords$y[i]- coords$y)^2)
    # Vector of probabilities of seed dispersal assuming seed dispersal distances
    # are exponentially distributed. The jth element is the probabilty that the
    # focal sampling point receives a seed from point j.
    pr_dispersal <- dexp(dist_vector, rate = 1/mean_dispersal_distance)
    pr_dispersal <- pr_dispersal / sum(pr_dispersal)
    # Sample where seeds come from in each generation.
    seed_disp_draws[[i]] <- sample(
      1:(grid_size^2),
      size    = n_generations,
      replace = T,
      prob    = pr_dispersal
    )
  }
  if(verbose) close(pb)
  # Coerce seed_disp_draws from a list to a matrix, with a row for each sampling
  # point and a column for each generation.
  seed_disp_draws <- do.call('rbind', seed_disp_draws)

  # OUTCROSSING AND TRACKING GENOTYPES
  # This loop over generations does two things:
  # (1) Uses the seed dispersal events simulated above to put new genotypes into
  # each sampling point.
  # (2) Chooses sampling points at random to receive outcrossed pollen and
  # and create a new offspring (unqiue) genotype at those sampling points.
  if(verbose) cat("Simulating outcrossing and updating genotypes based on dispersal events.\n")
  for(g in 2:n_generations){
    # Update the new generation by drawing seeds from the previous generation.
    pop[[g]] <- pop[[g-1]][seed_disp_draws[,g-1]]

    # Choose plants at random to receive outcrossed pollen.
    cross01 <- rbinom(grid_size^2, 1, outcrossing_rate)
    cross01 <- as.logical(cross01)
    # Give these plants a new unique genotype by appending generation number and
    # and integer from 1 to the number of outcrossers.
    pop[[g]][cross01] <- paste(pop[[g]][cross01],"_", g, ".", 1:sum(cross01), sep="")
  }
  # Each population is currently a vector. Coerce to a grid.
  pop <- lapply(pop, function(x) matrix(x, nrow=grid_size))

  t1 <- proc.time()[3]
  if(verbose) cat("\nSimulation complete after", round(t1-t0, 2), "seconds.\n")
  # output the population history
  if(return_all){
    return(pop)
  } else {
    return(pop[[g]])
  }
}
