#' Simulate a population
#'
#' Simulate a single population through time.
#'
#'`sim_population` simulates a population on an evenly-spaced grid of sampling
#'points through time:
#'
#'\enumerate{
#'    \item{Sampling points are populated at random with a set of unique
#'    starting genotypes.}
#'    \item{In the next generation, each sampling point is filled by one seed
#'    from the same or a different sampling point from a previous generation.
#'    Whichever genotype was at the donor sampling point now occupies the focal
#'    sampling point. I assume seed dispersal distances are exponentially
#'    distributed with some mean value that is used as an input parameter.}
#'    \item{To allow for a seed bank, each seed can be drawn from any previous
#'    generation. The generation is drawn from a poisson distribution.}
#'    \item{Each seed has some probability of having been the product of an
#'    outcrossing event. If so, it is assigned a new unique genotype.}
#' }
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
#' Defaults to 50
#' @param dormancy_rate Float >0. Strength of the seedbank. This is the rate
#' parameter for the poisson distribution. Giving zero will simulate no seed
#' bank.
#' @param verbose Logical. If TRUE, messages about simulation progress will be
#' printed. Defaults to TRUE.
#' @param return_all Logical. If TRUE, the entire population history will be
#' returned as a list. If FALSE, only the final generation is returned as a
#' vector. Defaults to TRUE.
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
  n_starting_genotypes = 50,
  dormancy_rate=0,
  verbose = TRUE,
  return_all = TRUE
  ){

  t0 <- proc.time()[3] # Record the time when the simulation started.
  if(verbose) cat("Initialising population.\n")

  # Initialise a set of starting genotypes.
  iggs <- paste("g", 1:n_starting_genotypes, sep = "")

  # We will keep track of the genotype at each sampling point using a list with
  # one element for each generation
  pop <- matrix(NA, nrow=n_generations, ncol=grid_size^2)
  # pop <- vector('list', n_generations)
  # Initialise the first generation with a random draw from the genotypes.
  pop[1,] <- sample(iggs, size = grid_size^2, replace = T)
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
  seed_disp_draws <- do.call('cbind', seed_disp_draws)

  # OUTCROSSING, SEED BANK AND TRACKING GENOTYPES
  # This loop over generations does three things:
  # (1) Uses the seed dispersal events simulated above to put new genotypes into
  # each sampling point.
  # (2) Draws those genotypes from the seed bank by choosing which of the
  # previous generations each seed should come from.
  # (3) Chooses sampling points at random to receive outcrossed pollen and
  # and create a new (unique) offspring genotype at those sampling points.
  if(verbose) cat("Simulating outcrossing and updating genotypes based on dispersal events.\n")
  for(g in 2:n_generations){
    # positions of seed donors for each sampling point.
    donor_positions <- seed_disp_draws[g-1,]

    # Choose how many generations back to draw seeds from the seed bank.
    seed_bank <- rpois(grid_size^2, lambda = dormancy_rate) + 1
    # If seeds are drawn from a generation before generation 1, set to 1.
    seed_bank[seed_bank >= g] <- g-1

    # Pull out the previous few generations from which to draw seeds.
    # Note that the direction has changed: the previous year is at the top,
    # and subsequent years following that.
    d <- pop[g- (1:max(seed_bank)),]

    # Pull seeds from the seed bank.
    if(all(seed_bank == 1)){
      # If every seed is drawn from the previous generation, just return that.
      # This has to be done separately because suddenly d is a vector, not a matrix.
      drawn_genotypes <- d[donor_positions]
    } else {
      # If there is a seed bank, d is a matrix
      # Pull out the correct genotypes from the rows of seed bank.
      drawn_genotypes <- sapply(1:grid_size^2, function(i) d[seed_bank[i], donor_positions[i]])
    }
    # Update the new generation.
    pop[g,] <- drawn_genotypes

    # Choose plants at random to receive outcrossed pollen.
    cross01 <- rbinom(grid_size^2, 1, outcrossing_rate)
    cross01 <- as.logical(cross01)
    # Give these plants a new unique genotype by appending generation number and
    # and integer from 1 to the number of outcrossers.
    pop[g, cross01] <- paste(pop[g,][cross01],"_", g, ".", 1:sum(cross01), sep="")
  }

  # Pop is currently a matrix of dimensions n_generations x grid_size^2
  # Coerce to a list of n_generations matrices, each of grid_size x grid_size
  pop <- lapply(1:n_generations, function(i) matrix(pop[i,], nrow=grid_size))

  t1 <- proc.time()[3]
  if(verbose) cat("\nSimulation complete after", round(t1-t0, 2), "seconds.\n")
  # output the population history
  if(return_all){
    return(pop)
  } else {
    return(pop[n_generations,])
  }
}
