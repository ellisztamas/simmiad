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
  years_to_sample,
  dormancy,
  habitat_width = (n_sample_points * sample_spacing * range_limit),
  population_size = density * habitat_width^2,
  selection_gradient
){
  data.frame(
    parameter = c(
      'mean_dispersal_distance',
      'outcrossing_rate',
      'n_generations',
      'n_starting_genotypes',
      'density',
      'dormancy',
      'nsims',
      'n_sample_points',
      'sample_spacing',
      'range_limit',
      'habitat_width',
      'population_size',
      'years_to_sample',
      'npairs',
      'selection'
    ),
    value = c(
      mean_dispersal_distance,
      outcrossing_rate,
      n_generations,
      n_starting_genotypes,
      density,
      dormancy,
      nsims,
      n_sample_points,
      sample_spacing,
      range_limit,
      habitat_width,
      population_size,
      years_to_sample,
      (n_sample_points * (n_sample_points-1))/2,
      selection_gradient
    )
  )
}
