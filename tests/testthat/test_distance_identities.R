positions <- sort(runif(16, -5,5))

transects <- list(
  rep(letters[1:4], each=4),
  rep(letters[1:4], 4)
)

test_that("distance_identities returns a data.frame of the right dimensions", {
  d <- distance_identities(transects, positions, labels = LETTERS[1:16])
  expect_true(class(d) == "data.frame")
  expect_true(ncol(d) == 4)
  expect_true(nrow(d) == (length(positions) * (length(positions)-1)) / 2)
})

test_that("distance_identities returns expected output if genotypes is vector not list", {
  # Check `matches` can only be zero or one.
  expect_true(all(
    distance_identities(transects[[1]], positions)$matches %in% c(0,1)
  ))
  # Check that only one comparson is made.
  expect_true(all(
    distance_identities(transects[[1]], positions)$n == 1
  ))
})

test_that("distance_identities throws an error if arguments are different lengths", {
  expect_error(distance_identities(transects, positions[-1]))
  expect_error(distance_identities(transects, positions, labels = 1:4))
})

