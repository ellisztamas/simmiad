  rs <- simmiad(
    mean_dispersal_distance = 1,
    outcrossing_rate = 0.01,
    n_generations = 10,
    n_starting_genotypes = 10,
    density = 3,
    n_sample_points = 5,
    sample_spacing = 5,
    nsims = 3,
    dormancy = 0.3
  )

test_that("simmiad() gives the expected output", {
  expect_true(class(rs) == "list")
  expect_length(rs, 7)
  expect_true(ncol(rs$clustering) == 10)
  expect_true(nrow(rs$clustering) == 3)
})
