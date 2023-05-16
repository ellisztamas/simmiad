test_that("initialise_population returns correct output for uniform genotypes", {
  uniform <- initialise_population(
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 200,
    pop_structure = "uniform"
  )
  expect_true(is.list(uniform))
  expect_true(length(uniform$geno) == 30)
  expect_true(nrow(uniform$coords) == 30)
  expect_true(ncol(uniform$coords) == 2)
})

test_that("initialise_population returns correct output for founding genotypes", {
  founders <- initialise_population(
    n_starting_genotypes = 10,
    population_size = 30,
    box_limit = 200,
    pop_structure = "clusters",
    mixing = 1
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
    pop_structure = "mvnorm",
    mixing = 1
  )
  expect_true(is.list(multiv))
  expect_true(length(multiv$geno) == 30)
  expect_true(nrow(multiv$coords) == 30)
  expect_true(ncol(multiv$coords) == 2)
})

test_that("initialise_population returns correct output for a specfic structure.", {
  pop_structure <- c( "G18", "G18", "G1", "G1", "G2", "G3", "G4", "G5", "G3", "G3", "G3", "G3",
                      "G6", "G7", "G8", "G3", "G3", "G9", "G3", "G10", "G11", "G12","G13", "G13",
                      "G13", "G13", "G11", "G11", "G14", "G15", "G13", "G16", NA, "G17", "G11",
                      "G17", "G17", "G19", "G19", "G20", "G20", "G20", "G20", "G17", "G17", "G20",
                      "G20", NA, "G13", "G13", "G21", "G17", "G22", "G11", "G17", "G23", "G22",
                      "G22", "G22", "G17", "G22", "G19", "G24", "G19", "G25", "G19", "G20", "G26",
                      "G22", "G25", "G27", "G28", "G28", "G28", "G28", "G29", "G30", "G30", "G30",
                      "G31", "G25", "G25", "G20", "G25", "G25", NA, "G20", "G25", NA, "G32", "G32",
                      NA, "G20", "G33", "G25", NA, "G25", "G19", "G34", "G25", NA, "G19")

  habitat_labels <- rep(LETTERS[1:10], rmultinom(1, 102, rep(1/10, 10)))

  mixing = 1
  n_sample_points=30
  sample_spacing=5
  real_transect_length <- length(pop_structure)
  # range_limit = 1 + 1/(real_transect_length-1) # Make sure this is one. Add a warning
  density<- 4
  box_limit <- ( real_transect_length * sample_spacing ) / 2
  # Given a density of plants per sq. metre and a size of the box, calculate
  # how many individuals you need
  population_size <- density * (2*box_limit)^2

  ip <- initialise_population(
    n_starting_genotypes = pop_structure,
    population_size = population_size,
    box_limit = box_limit,
    pop_structure = "hardcoded",
    mixing = 1,
    habitat_labels = habitat_labels
  )

  # Visualise the population. I don't know how to make this a formal test.
  # ix <- sample(1:nrow(ip$coords), 5000, replace = FALSE)
  # plot(ip$coords[ix,1], ip$coords[ix,2], col = substr(ip$geno[ix],2,5))

  expect_true(is.list(ip))
  expect_true(length(ip$geno) == population_size)
  expect_true(nrow(ip$coords) == population_size)
  expect_true(ncol(ip$coords) == 2)
  expect_true(
    length(attributes(ip)$habitat_labels) == length(habitat_labels)
  )
})


test_that(
  "Throws and error if `pop_structure='hardcoded'` but `n_starting_genotypes is not a vector.",
  {
    n_starting_genotypes = rep(1:5, each = 20)
    density <- 4
    sample_spacing <- 1
    box_limit <- length(n_starting_genotypes) * sample_spacing /2
    population_size <- round(density * (2*box_limit)^2)

    expect_error(
      initialise_population(
        n_starting_genotypes = "n_starting_genotypes",
        population_size = population_size,
        box_limit = box_limit,
        pop_structure = "hardcoded",
        mixing = 1
      ))
  })

test_that("initialise_population throws an error is `mixing` is required but not given.", {
  expect_error(
    initialise_population(
      n_starting_genotypes = 10, population_size = 30,box_limit = 200,
      pop_structure = "clusters",
      mixing = NULL
    )
  )
  expect_error(
    initialise_population(
      n_starting_genotypes = 10, population_size = 30,box_limit = 200,
      pop_structure = "mvnorm",
      mixing = NULL
    )
  )

})

test_that("initialise_population throws an error if the pop_structure is garbage.", {
  expect_error(
    initialise_population(
      n_starting_genotypes = 10,
      population_size = 30,
      box_limit = 200,
      pop_structure = "whoops"
    ),
    regexp = "`pop_structure` should be one of 'uniform', 'clusters' or 'mvnorm'."
  )
})
