#' Summarise output from `simmiad`
#'
#' Function to loop through folders and summarise simulation results over many
#' parameter combinations.
#'
#' @param base_folder String giving the folder containing simulation results.
#' This folder should contain subfolders for each combination of parameters run,
#' each with a set of CSV files for that analysis.
#' @param parameter String giving the output parameter to be summarised (i.e.
#' the output CSV file to be opened).
#'
#' @return This returns a list giving means and CIs for the final generation, plus
#' means through time across generations.
#' @author Tom Ellis
#' @export
summarise_simmiad <- function(base_folder, parameter){
  folders <- list.dirs(base_folder)

  over_reps    <- vector('list', length(folders))
  through_time <- vector('list', length(folders))

  for(f in 2:length(folders)){ # start from 2 because element 1 is the parent folder
    # Open the parameter file
    p <-  read.csv(
      paste(folders[f], "/parameters.csv", sep=""),
      header = FALSE,
      col.names = c('parameter', 'value')
    ) %>%
      filter(parameter %in% c(
        'mean_dispersal_distance', 'outcrossing_rate', 'n_generations', 'n_starting_genotypes', 'density', 'dormancy'
      ))
    # Open the data file
    df <- read.csv(paste(folders[f], "/", parameter, sep=""), header=FALSE)

    # Get the mean and CIs for the final generation.
    final_generation <- df[, ncol(df)]
    over_reps[[f]] <- c(
      p$value,
      mean = mean(final_generation, na.rm = T),
      quantile(final_generation, c(0.025, 0.975), na.rm = T)
    )
    # Get the mean of each generation through time
    through_time[[f]] <- cbind(
      sapply(p$value, rep, ncol(df)), # parameter values
      generation = 1:ncol(df),
      mean = colMeans(df, na.rm = T),
      lower = apply(df, 2, quantile, 0.025, na.rm=T),
      upper = apply(df, 2, quantile, 0.975, na.rm=T)
    )

  }
  # Concatenate those lists
  over_reps <- as.data.frame(do.call("rbind", over_reps))
  colnames(over_reps) <- c(
    as.character(p$parameter),
    'mean', 'lower', 'upper'
  )
  through_time <- as.data.frame(do.call('rbind', through_time))
  colnames(through_time) <- c(
    as.character(p$parameter),
    'generation',
    'mean', 'lower', 'upper'
  )

  list(
    over_reps   = over_reps,
    through_time = through_time
  )
}
