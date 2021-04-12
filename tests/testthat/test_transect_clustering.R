positions <- 1:16

geno_clustered <- rep(letters[1:4], each=4)
geno_dispersed <- rep(letters[1:4], 4)

test_that("Clustered genotypes give a negative covariance with distance", {
  expect_true(transect_clustering(geno_clustered, positions)[3] < 0)
})

test_that("Dispersed genotypes give a positive covariance with distance", {
  expect_true(transect_clustering(geno_dispersed, positions)[3] > 0)
})

