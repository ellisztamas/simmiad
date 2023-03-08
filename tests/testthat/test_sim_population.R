mean_dispersal_distance = 1
outcrossing_rate = 0.01
n_starting_genotypes = 10
density <- 3
n_generations <- 10
range_limit <- 1.5
n_sample_points = 5
sample_spacing = 5

test_that("sim_population returns a list of the right length", {
  sm <- sim_population(
    mean_dispersal_distance = mean_dispersal_distance,
    outcrossing_rate = outcrossing_rate,
    n_generations = n_generations,
    n_starting_genotypes = n_starting_genotypes,
    density = density,
    n_sample_points = n_sample_points,
    sample_spacing = sample_spacing,
    range_limit = range_limit,
    dormancy = 0.3
  )

  expect_true(class(sm) == "list")
  expect_true(length(sm) == n_generations)
  expect_true(all(sapply(sm, length) > 0))
})

test_that("sim_population throws an error for range limit < 1", {
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
    "range_limit > 1 is not TRUE")
})

test_that("sim_population returns sensible data for `method='mvnorm'.", {
  sm <- sim_population(
    mean_dispersal_distance = mean_dispersal_distance,
    outcrossing_rate = outcrossing_rate,
    n_generations = n_generations,
    n_starting_genotypes = n_starting_genotypes,
    density = density,
    n_sample_points = n_sample_points,
    sample_spacing = sample_spacing,
    range_limit = range_limit,
    dormancy = 0.3,
    method = "mvnorm"
  )

  expect_true(class(sm) == "list")
  expect_true(length(sm) == n_generations)
  expect_true(all(sapply(sm, length) > 0))
})
