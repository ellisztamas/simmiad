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

test_that("sim_population returns sensible data for `pop_structure='mvnorm'.", {
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
    pop_structure = "mvnorm"
  )

  expect_true(class(sm) == "list")
  expect_true(length(sm) == n_generations)
  expect_true(all(sapply(sm, length) > 0))
})

test_that(
  "sim_population returns sensible data for `pop_structure='hardcoded'.",
  {
    set.seed(3)
    pop_structure <- c( "G18", "G18", "G1", "G1", "G2", "G3", "G4", "G5", "G3", "G3", "G3", "G3",
                        "G6", "G7", "G8", "G3", "G3", "G9", "G3", "G10", "G11", "G12","G13", "G13",
                        "G13", "G13", "G11", "G11", "G14", "G15", "G13", "G16", NA, "G17", "G11",
                        "G17", "G17", "G19", "G19", "G20", "G20", "G20", "G20", "G17", "G17", "G20",
                        "G20", NA, "G13", "G13", "G21", "G17", "G22", "G11", "G17", "G23", "G22",
                        "G22", "G22", "G17", "G22", "G19", "G24", "G19", "G25", "G19", "G20", "G26",
                        "G22", "G25", "G27", "G28", "G28", "G28", "G28", "G29", "G30", "G30", "G30",
                        "G31", "G25", "G25", "G20", "G25", "G25", NA, "G20", "G25", NA, "G32", "G32",
                        NA, "G20", "G33", "G25", NA, "G25", "G19", "G34", "G25", NA, "G19")
    mixing = 1
    n_sample_points=30
    sample_spacing=5
    real_transect_length <- length(pop_structure)
    range_limit = 1 + 1/(real_transect_length-1) # Make sure this is one. Add a warning
    density<- 4

    sm <- sim_population(
      mean_dispersal_distance = mixing,
      outcrossing_rate = 0.04,
      n_generations = 2,
      n_starting_genotypes = pop_structure,
      density = density,
      n_sample_points = n_sample_points,
      sample_spacing = sample_spacing,
      range_limit = range_limit,
      dormancy = 0.3,
      pop_structure = "hardcoded"
    )
    # Check the output formatting
    expect_true(class(sm) == "list")
    expect_true(length(sm) == 2)
    expect_true(all(sapply(sm, length) > 0))
    # Check that no single genotype makes up the whole sample
    # That's what would happen if the transect was going along the wrong axis
    expect_gt(
      length(  unique(na.exclude(sm[[1]])) ),
      1)
  })
