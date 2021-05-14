ip <- initialise_population(
  mean_dispersal_distance = 1,
  n_starting_genotypes = 10,
  population_size = 30,
  box_limit = 20
)

test_that(
  "update_population gives expected output", {
  pop <- update_population(
    seed_rain = ip,
    seed_bank = ip,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.1,
    dormancy = 0.3,
    generation = 2,
    box_limit = 20
  )
  expect_true(is.list(pop))
  expect_true(length(pop) == 2)
  expect_true(length(pop$geno) == 30)
  expect_true(nrow(pop$coords) == 30)
  })
