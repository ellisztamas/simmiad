test_that("relative_fitness returns 1 is selection is 0",{
  ip <- initialise_population(
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 20,
    var_w = 0
  )
  w <- relative_fitness(ip)
  expect_true(all(
    w ==1
    ))
})

test_that("relative_fitness returns means relative fitness of 1",{
  ip <- initialise_population(
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 20,
    var_w = 1
  )
  w <- relative_fitness(ip)
  expect_equal( mean(w), 1 )
})

test_that("When phenotypes are NA, relative_fitness returns 0", {
  ip <- initialise_population(
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 20,
    var_w = 1
  )

  ip$phenotype[c(1,5,9)] <- NA
  w <- relative_fitness(ip)
  expect_true( all(w[c(1,5,9)] == 0) )
})
