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

pos_diff <- function(a,b){
  sqrt(
    ((a[,1] - b[,1])^2) + ((a[,2] - b[,2])^2)
  )
}



test_that(
  "shift_positions moves things by the amount I expect",
  {
    # Baseline population
    coords = matrix(
      runif(n = 100000 * 2, min = -box_limit, max = box_limit),
      ncol=2
    )
    # Offset by 0.1 metre on average
    offset <- 0.1
    coords_new <- shift_positions(
      coords,
      mean_dispersal_distance = offset,
      box_limit = box_limit * 10
    )
    expect_equal(
      mean(pos_diff(coords, coords_new)), offset, tolerance = 0.01
    )
    # Offset by 1 metre on average
    offset <- 1
    coords_new <- shift_positions(
      coords,
      mean_dispersal_distance = offset,
      box_limit = box_limit * 10
    )
    expect_equal(
      mean(pos_diff(coords, coords_new)), offset, tolerance = 0.01
    )
    # Offset by 10 metres on average
    offset <- 10
    coords_new <- shift_positions(
      coords,
      mean_dispersal_distance = offset,
      box_limit = box_limit * 10
    )
    expect_equal(
      mean(pos_diff(coords, coords_new)), offset, tolerance = 0.01
    )

    # Offset by 0 metres on average
    offset <- 0
    coords_new <- shift_positions(
      coords,
      mean_dispersal_distance = offset,
      box_limit = box_limit * 10
    )
    expect_equal(
      mean(pos_diff(coords, coords_new)), offset, tolerance = 0.01
    )
  }
  )
