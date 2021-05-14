test_that("initialise_population returns correct output", {
  ip <- initialise_population(
    mean_dispersal_distance = 1,
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 20
  )
  expect_true(is.list(ip))
  expect_true(length(ip$geno) == 30)
  expect_true(nrow(ip$coords) == 30)
  expect_true(ncol(ip$coords) == 2)
})
