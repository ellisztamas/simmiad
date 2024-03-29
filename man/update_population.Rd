% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_population.R
\name{update_population}
\alias{update_population}
\title{Update population}
\usage{
update_population(
  seed_rain,
  seed_bank,
  mean_dispersal_distance,
  outcrossing_rate,
  dormancy,
  generation,
  box_limit
)
}
\arguments{
\item{seed_rain}{List giving (1) genotype labels and (2) coordinates of
each individual in the population in the previous generation.}

\item{seed_bank}{List giving (1) genotype labels and (2) coordinates of
each individual in the seed bank.}

\item{mean_dispersal_distance}{Float >0. Mean seed dispersal distance. The
reciprocal of this is used as the rate parameter to draw from the exponential
distribution.}

\item{outcrossing_rate}{Float between 0 and 1. Probability that an individual
is outcrossed.}

\item{dormancy}{Float between 0 and 1. Probability that a seedling is drawn
from the seed bank. Seedlings are drawn from the prior generation with
probability 1-dormancy.}

\item{generation}{Int indicating the generation in the simulation. This is
used to label new genotypes generated by outcrossing.}

\item{box_limit}{Float >1 giving half the circumference of the torus.}
}
\value{
List giving (1) genotype labels and (2) coordinates of
each individual to found the next generation.
}
\description{
Draw individuals to found a new generation. These are drawn from the seed
rain from the previous generation or from the seed bank.
}
\author{
Tom Ellis
}
