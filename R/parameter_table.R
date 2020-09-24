#' Data.frame of parameter values
#'
#' @inheritParams simmiad
#' @param habitat_width Length of the edge of the habitat square
#' @param population_size Number of individual plants in the population,
#' derived from the habitat size and population density
#' @author Tom Ellis
#' @return A data.frame giving parameter names and values
parameter_table <- function(
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes,
  density,
  nsims,
  n_sample_points,
  sample_spacing,
  range_limit,
  how_far_back,
  habitat_width = (n_sample_points * sample_spacing * range_limit),
  population_size = density * habitat_width^2
){
  data.frame(
    parameter = c(
      'mean_dispersal_distance',
      'outcrossing_rate',
      'n_generations',
      'n_starting_genotypes',
      'density',
      'nsims',
      'n_sample_points',
      'sample_spacing',
      'range_limit',
      'habitat_width',
      'population_size',
      'how_far_back',
      'npairs'
    ),
    value = c(
      mean_dispersal_distance,
      outcrossing_rate,
      n_generations,
      n_starting_genotypes,
      density,
      nsims,
      n_sample_points,
      sample_spacing,
      range_limit,
      habitat_width,
      population_size,
      how_far_back,
      (n_sample_points * (n_sample_points-1))/2
    )
  )
}
