n_starting_genotypes <- 10
box_limit <- 100

coords = matrix(
  runif(n = n_starting_genotypes * 2, min = -box_limit, max = box_limit),
  ncol=2
)

test_that("shift_positions gives expected output", {
  sf <- shift_positions(
    coords,
    mean_dispersal_distance = 1,
    box_limit = box_limit
  )
  expect_true(is.matrix(sf))
  expect_true(nrow(sf) == n_starting_genotypes)
  expect_true(ncol(sf) == 2)
})
