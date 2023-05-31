ip <- initialise_population(
  n_starting_genotypes = 10,
  population_size = 30,
  box_limit = 20
)

test_that("update_population gives expected output", {
  pop <- update_population(
    seed_rain = ip,
    seed_bank = ip,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.1,
    dormancy = 0.3,
    generation = 2,
    box_limit = 20,
    selection_gradient = 0
  )
  expect_true(is.list(pop))
  expect_true(length(pop) == 3)
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
    box_limit = 20,
    selection_gradient = 0
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
    box_limit = 20,
    selection_gradient = 0
  )
  expect_true(all(pop4$geno == 2))
})

test_that("Phenotypes are inherited correctly.", {
  pop <- update_population(
    seed_rain = ip,
    seed_bank = ip,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.1,
    dormancy = 0.3,
    generation = 2,
    box_limit = 20,
    selection_gradient = 0
  )

  non_crossed_genotypes <- grep("_", pop$geno, invert = TRUE)
  original_z <- ip$phenotype[ match( pop$geno[non_crossed_genotypes], ip$geno) ]
  new_z <- pop$phenotype[non_crossed_genotypes]
  expect_true(
    all(original_z == new_z)
  )
})

test_that("There is a change in phenotype under directional selection in update_population.", {
  # Check that phenotypes are higher in generation 2 when selection is positive
  pop <- update_population(
    seed_rain = ip,
    seed_bank = ip,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.1,
    dormancy = 0.3,
    generation = 2,
    box_limit = 20,
    selection_gradient = 5
  )
  gen1 <- mean(ip$phenotype, na.rm = TRUE)
  gen2 <- mean(pop$phenotype, na.rm = TRUE)
  expect_gt(gen2, gen1)

  # Check that phenotypes stay roughly the same when selection is zero
  pop <- update_population(
    seed_rain = ip,
    seed_bank = ip,
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.1,
    dormancy = 0.3,
    generation = 2,
    box_limit = 20,
    selection_gradient = 0
  )
  gen1 <- mean(ip$phenotype, na.rm = TRUE)
  gen2 <- mean(pop$phenotype, na.rm = TRUE)
  expect_equal(gen2, gen1, tolerance = 0.2)

})
