test_that(
  "habitat_clustering gives expected values", {
    # No concordance between habitat and genotype
    t1 <- habitat_clustering(
      genotype = paste0("g", 1:12),
      habitat = rep(paste0('h', LETTERS[1:3]), each = 4)
    )
    expect_equal(t1, 0)
    # Perfect concordance between habitat and genotype
    t2 <- habitat_clustering(
      genotype = rep(paste0('g', LETTERS[1:3]), each =  4),
      habitat = rep(paste0('h', LETTERS[1:3]), each = 4)
    )
    expect_equal(t2, 1)
  })

test_that(
  "habitat_clustering is not affected by NA",{
    expect_equal(
      habitat_clustering(
        genotype = c(rep(paste0('g', LETTERS[1:3]), each =  4), NA),
        habitat = c(rep(paste0('h', LETTERS[1:3]), each = 4), "hD")
      ),
      1
    )
  }
)
