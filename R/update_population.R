#' Update population
#'
#' Draw individuals to found a new generation. These are drawn from the seed
#' rain from the previous generation or from the seed bank.
#'
#' @param seed_rain List giving (1) genotype labels and (2) coordinates of
#' each individual in the population in the previous generation.
#' @param seed_bank List giving (1) genotype labels and (2) coordinates of
#' each individual in the seed bank.
#' @param mean_dispersal_distance Float >0. Mean seed dispersal distance. The
#' reciprocal of this is used as the rate parameter to draw from the exponential
#' distribution.
#' @param outcrossing_rate Float between 0 and 1. Probability that an individual
#' is outcrossed.
#' @param dormancy Float between 0 and 1. Probability that a seedling is drawn
#' from the seed bank. Seedlings are drawn from the prior generation with
#' probability 1-dormancy.
#' @param box_limit Float >1 giving half the circumference of the torus.
#' @param generation Int indicating the generation in the simulation. This is
#' used to label new genotypes generated by outcrossing.
#'
#' @return List giving (1) genotype labels and (2) coordinates of
#' each individual to found the next generation.
#' @author Tom Ellis
update_population <- function(
  seed_rain,
  seed_bank,
  mean_dispersal_distance,
  outcrossing_rate,
  dormancy,
  generation,
  box_limit
  ){

  stopifnot(is.list(seed_rain) & is.list(seed_bank))
  stopifnot(is.vector(seed_rain$geno) & is.vector(seed_bank$geno))
  stopifnot(ncol(seed_rain$coords) ==2 & ncol(seed_bank$coords) ==2)
  stopifnot(length(seed_rain$geno) == nrow(seed_rain$coords))
  stopifnot(length(seed_bank$geno) == nrow(seed_bank$coords))
  stopifnot(length(seed_rain$geno) == nrow(seed_rain$coords))
  stopifnot(outcrossing_rate >= 0 & outcrossing_rate <= 1)
  stopifnot(dormancy >= 0 & dormancy <= 1)
  stopifnot(generation == round(generation))
  stopifnot(generation > 0)

  # Relative fitness of each genotype
  seed_rain$fitness <- relative_fitness(seed_rain)
  seed_bank$fitness <- relative_fitness(seed_bank)

  # Data.frame containing indexes of plants from seed_rain and seed_bank,
  # with a vector of probabilities of being chosen
  ix <- rbind(
    data.frame(
      i = 1:length(seed_rain$geno),
      p = (1-dormancy ) * seed_rain$fitness
    ),
    data.frame(
      i = (length(seed_bank$geno) + 1) : (2*length(seed_bank$geno)),
      p = dormancy * seed_bank$fitness
    )
  )
  ix$p <- ix$p / sum(ix$p) # normalise probabilties to sum to one.

  # Draw indices for genotypes to found the next generation.
  ix <- sample(
    x = ix$i,
    size = length(seed_rain$geno),
    prob = ix$p,
    replace = TRUE
  )

  # Create a new population
  pop <- list(
    geno   = c(seed_rain$geno, seed_bank$geno)[ix],
    coords = do.call(rbind, list(seed_rain$coords, seed_bank$coords))[ix,],
    phenotype = c(seed_rain$phenotype, seed_bank$phenotype)[ix]
  )

  # Choose plants at random to receive outcrossed pollen.
  cross01 <- rbinom(n = length(pop$geno), 1, outcrossing_rate)
  cross01 <- as.logical(cross01)
  # Give these plants a new unique genotype by appending generation number and
  # and integer from 1 to the number of outcrossers.
  pop$geno[cross01] <- paste(pop$geno[cross01], "_", generation-1, ".", 1:sum(cross01), sep = "")
  pop$phenotype[cross01] <- rnorm(
    n    = sum(cross01),
    mean = mean(pop$phenotype, na.rm=TRUE),
    sd   = sd(pop$phenotype, na.rm=TRUE) / 2
    )

  # Peturb positions
  pop$coords <- shift_positions(
    pop$coords,
    mean_dispersal_distance = mean_dispersal_distance,
    box_limit = box_limit
  )

  pop
}
