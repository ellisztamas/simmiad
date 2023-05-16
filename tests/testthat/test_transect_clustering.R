positions <- 1:16

geno_clustered <- rep(letters[1:4], each=4)
geno_dispersed <- rep(letters[1:4], 4)

test_that("Across the whole transect, 20% of pairs match", {
  expect_true(  transect_clustering(geno_clustered, positions, 16)[3] == 24/120  )
  expect_true(  transect_clustering(geno_dispersed, positions, 16)[3] == 24/120  )
})


test_that("Within four metres there is more clustering in geno_clustered than geno_dispersed", {
  expect_true(  transect_clustering(geno_clustered, positions, 4)[3] == 24/120  )
  expect_true(  transect_clustering(geno_dispersed, positions, 4)[3] == 0  )
})



