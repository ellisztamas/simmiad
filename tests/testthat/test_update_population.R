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

test_that("Check dormancy setting works", {
  pop1 <- pop2 <- ip
  pop1$geno <- rep("1", 30)
  pop2$geno <- rep("2", 30)
  # Setting dormancy to 0 should return only genotypes of 1
  pop3 <- update_population(
    pop1,
    pop2,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0,
    dormancy = 0,
    generation = 2,
    box_limit = 20
    )
  expect_true(all(pop3$geno == 1))
  # Setting dormancy to 2 should return only genotypes of 2
  pop4 <- update_population(
    pop1,
    pop2,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0,
    dormancy = 1,
    generation = 2,
    box_limit = 20
  )
  expect_true(all(pop4$geno == 2))
})


