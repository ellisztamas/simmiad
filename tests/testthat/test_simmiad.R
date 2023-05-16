test_that("simmiad() gives the expected output", {
  rs <- simmiad(
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.01,
    n_generations = 10,
    n_starting_genotypes = 10,
    density = 3,
    n_sample_points = 5,
    sample_spacing = 5,
    nsims = 3,
    dormancy = 0.3,
    progress = FALSE
  )

  expect_true(class(rs) == "list")
  expect_length(rs, 8)
  expect_true(ncol(rs$clustering) == 10)
  expect_true(nrow(rs$clustering) == 3)
  expect_true(
    all( rs$stability[,1] == 1 ) # The first column should all be one, because there is no change
  )
})

test_that("simmiad() gives the expected output when `pop_structure='mvnorm'`.", {
  rs <- simmiad(
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.01,
    n_generations = 10,
    n_starting_genotypes = 10,
    density = 3,
    n_sample_points = 5,
    sample_spacing = 5,
    nsims = 3,
    dormancy = 0.3,
    progress = FALSE,
    pop_structure = 'mvnorm',
    mixing = 3
  )

  expect_true(class(rs) == "list")
  expect_length(rs, 8)
  expect_true(ncol(rs$clustering) == 10)
  expect_true(nrow(rs$clustering) == 3)
})

test_that("simmiad() gives the expected output when `pop_structure='mvnorm'`.", {
  rs <- simmiad(
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.01,
    n_generations = 10,
    n_starting_genotypes = 10,
    density = 3,
    n_sample_points = 5,
    sample_spacing = 5,
    nsims = 3,
    dormancy = 0.3,
    progress = FALSE,
    pop_structure = 'mvnorm',
    mixing = 3
  )

  expect_true(class(rs) == "list")
  expect_length(rs, 8)
  expect_true(ncol(rs$clustering) == 10)
  expect_true(nrow(rs$clustering) == 3)
})

test_that("simmiad() gives the expected output when `pop_structure='hardcoded'`.", {
  set.seed(11)
  pop_structure <- c(
    "G18", "G18", "G1", "G1", "G2", "G3", "G4", "G5", "G3", "G3", "G3", "G3",
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
  range_limit = 1 + 1/(real_transect_length-1) # Make sure this is one. Add a warning
  density<- 4
  nsims<- 3

  rs <- simmiad(
    mean_dispersal_distance = mixing,
    outcrossing_rate = 0.04,
    n_generations = 2,
    n_starting_genotypes = pop_structure,
    density = density,
    n_sample_points = n_sample_points,
    sample_spacing = 5,
    nsims = 3,
    dormancy = 0.3,
    progress = FALSE,
    pop_structure = 'hardcoded'
  )

  expect_true(class(rs) == "list")
  expect_length(rs, 8)
  expect_true(ncol(rs$clustering) == 2)
  expect_true(nrow(rs$clustering) == 3)
  expect_equal(
    rs$parameters$value[rs$parameters$parameter=="n_starting_genotypes"], length(pop_structure)
    )
  expect_equal( nrow(rs$distance_identity), 3 )
  expect_equal( ncol(rs$distance_identity), n_sample_points-1 )
})


test_that("simmiad() gives the expected output with a vector of habitat labels.", {
  set.seed(11)
  pop_structure <- c(
    "G18", "G18", "G1", "G1", "G2", "G3", "G4", "G5", "G3", "G3", "G3", "G3",
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
  range_limit = 1 + 1/(real_transect_length-1) # Make sure this is one. Add a warning
  density<- 4
  nsims<- 3

  rs <- simmiad(
    mean_dispersal_distance = mixing,
    outcrossing_rate = 0.04,
    n_generations = 2,
    n_starting_genotypes = pop_structure,
    density = density,
    n_sample_points = n_sample_points,
    sample_spacing = 5,
    nsims = 3,
    dormancy = 0.3,
    progress = FALSE,
    pop_structure = 'hardcoded',
    habitat_labels = habitat_labels
  )

  expect_true( is.matrix(rs$clustering_by_habitat) )
})
