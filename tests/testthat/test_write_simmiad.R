test_that(
  "write_simmiad writes the correct files.",
  {
    outdir <- tempdir()
    # double check the folder doesn't exist before running this test
    if( dir.exists(outdir) ) unlink(outdir, recursive = TRUE)

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
      var_w = 1
    )
    write_simmiad(rs, outdir)

    # Check the correct files are created.
    expected_files <- c(
      "clustering.csv", "count_NAs.csv", "di_by_year.csv", "distance_identity.csv",
      "matching_pairs.csv", "n_genotypes.csv", "parameters.csv", "stability.csv")
    expect_true(
      all( list.files(outdir) == expected_files )
    )
    # Check the files are non-zero size
    expect_true( all(
      sapply( Sys.glob( paste0(outdir,"/*csv") ), function(f) file.info(f)$size > 0)
    ) )
    # Remove the test folder so as not to interfere with future tests.
    unlink(outdir, recursive = TRUE)
  }
)
