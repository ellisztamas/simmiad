update_population <- function(
  pop,
  mean_dispersal_distance,
  range_limit,
  size = nrow(pop$coords)
  ){
  stopifnot(is.list(pop))
  stopifnot(is.vector(pop$geno))
  stopifnot(ncol(coords) ==2)
  stopifnot(length(pop$geno) == nrow(pop$coords))
  stopifnot(mean_dispersal_distance > 0)
  # Draw a sample of individuals
  ix <- sample(
    x = 1:size,
    size = population_size,
    replace = TRUE
  )
  #
  list(
    geno   = pop$geno[ix],
    coords = shift_positions(
      pop$coords[ix,],
      mean_dispersal_distance = mean_dispersal_distance,
      range_limit = range_limit
    )
  )
}
