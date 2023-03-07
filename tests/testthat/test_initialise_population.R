test_that("initialise_population returns correct output for uniform genotypes", {
  uniform <- initialise_population(
    mean_dispersal_distance = 1,
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 200,
    method = "uniform"
  )
  expect_true(is.list(uniform))
  expect_true(length(uniform$geno) == 30)
  expect_true(nrow(uniform$coords) == 30)
  expect_true(ncol(uniform$coords) == 2)
})

test_that("initialise_population returns correct output for founding genotypes", {
  founders <- initialise_population(
    mean_dispersal_distance = 1,
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 200,
    method = "clusters"
  )
  expect_true(is.list(founders))
  expect_true(length(founders$geno) == 30)
  expect_true(nrow(founders$coords) == 30)
  expect_true(ncol(founders$coords) == 2)
})

test_that("initialise_population returns correct output for 'mvnorm'.", {
  multiv <- initialise_population(
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 200,
    method = "mvnorm"
  )
  expect_true(is.list(multiv))
  expect_true(length(multiv$geno) == 30)
  expect_true(nrow(multiv$coords) == 30)
  expect_true(ncol(multiv$coords) == 2)
})

test_that("initialise_population throws an error if the method is garbage.", {
  expect_error(
    initialise_population(
      mean_dispersal_distance = 1,
      n_starting_genotypes = 10,
      population_size = 30,
      box_limit = 200,
      method = "whoops"
    ),
    regexp = "`method` should be one of 'uniform', 'clusters' or 'mvnorm'."
  )
})
