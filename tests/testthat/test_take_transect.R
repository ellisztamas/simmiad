# Define a squre habitat with width=10 centred on zero
range_limit <- 5
# This gives a transect at
# x = -4 -3 -2 -1  0  1  2  3  4
# y = 0
n_sample_points = 9
sample_spacing = 1

# Generate a dense population that should run without error
set.seed(24)
population_size <- 1000
coords <- matrix(
  runif(
    n = population_size,
    min = -range_limit,
    max = range_limit
  ),
  ncol=2
)
# Check the function returns what I expect for the dense population
test_that("take_transect returns a vector of integers for a dense population.",{
  tr <- take_transect(coords, n_sample_points, sample_spacing)
  expect_true(class(tr) == "integer")
  expect_equal(length(tr), n_sample_points)
  expect_true(all(tr < population_size))
})

# A population parallel to the transect, but 1m away
coords<- cbind(1, -4:4)
test_that("plants 1m from sampling points are sampled", {
  tr <- take_transect(coords, n_sample_points, sample_spacing)
  expect_equal(tr, 1:9)
})
# A population parallel to the transect, but 2m away
coords<- cbind(2, -4:4)
test_that("plants 2m from sampling points are not sampled", {
  tr <- take_transect(coords, n_sample_points, sample_spacing)
  expect_true(all(is.na(tr)))
})


