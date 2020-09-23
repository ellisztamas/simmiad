mean_dispersal_distance = 0.5
outcrossing_rate = 0.01
n_starting_genotypes = 40
density = 1
n_generations <- 3
range_limit <- 0.8
n_sample_points = 30
sample_spacing = 5

sm <- sim_population(
  mean_dispersal_distance = mean_dispersal_distance,
  outcrossing_rate = outcrossing_rate,
  n_generations = n_generations,
  n_starting_genotypes = n_starting_genotypes,
  density = density,
  n_sample_points = n_sample_points,
  sample_spacing = sample_spacing,
  range_limit = range_limit
)

test_that("sim_population returns a list of the right length", {
  expect_true(class(sm) == "list")
  expect_true(length(sm) == n_generations)
  expect_true(all(sapply(sm, length) > 0))
})

test_that("sim_population throws an error for range limit < 1")
expect_error(
  sim_population(
    mean_dispersal_distance = mean_dispersal_distance,
    outcrossing_rate = outcrossing_rate,
    n_generations = n_generations,
    n_starting_genotypes = n_starting_genotypes,
    density = density,
    n_sample_points = n_sample_points,
    sample_spacing = sample_spacing,
    range_limit = 0.8
  ),
  "range_limit must be greater than one"
)
