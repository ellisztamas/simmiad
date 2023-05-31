test_that(
  "habitat_clustering gives expected values", {
  # No concordance between habitat and genotype
  t1 <- habitat_clustering(
    genotype = paste0("g", 1:12),
    habitat = rep(paste0('h', LETTERS[1:3]), each = 4)
  )
  expect_true(is.na(t1))
  # Perfect concordance between habitat and genotype
  t2 <- habitat_clustering(
    genotype = rep(paste0('g', LETTERS[1:3]), each =  4),
    habitat = rep(paste0('h', LETTERS[1:3]), each = 4)
  )
  expect_equal(t2, 1)
  })
