% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_population.R
\name{sim_population}
\alias{sim_population}
\title{Simulate a population}
\usage{
sim_population(
  population_size,
  mean_dispersal_distance,
  outcrossing_rate,
  n_generations,
  n_starting_genotypes,
  density,
  n_sample_points,
  sample_spacing
)
}
\arguments{
\item{population_size}{Int >0. Number of individuls in the population.}

\item{mean_dispersal_distance}{Float >0. Mean seed dispersal distance. The
reciprocal of this is used as the rate parameter to draw from the exponential
distribution.}

\item{outcrossing_rate}{Float between 0 and 1. Probability that an individual
is outcrossed.}

\item{n_generations}{Int >12. Number of generations to run the simulations.}

\item{n_starting_genotypes}{Int >0. Number of initial genotypes to start with.
Defaults to 50}

\item{density}{Float >0. Average density of plants per square metre.}

\item{n_sample_points}{Number of points to sample along the transect.}

\item{sample_spacing}{Distance between sampling points.}
}
\value{
A list of genotypes recorded at each sampling point in each
generation.
}
\description{
Simulate a single population through time on a torus, and sample a transect
for each generation.
}
\details{
\code{sim_population} simulates a population of annuals dispersing in continuous
space following exponentially distributed seed dispersal distances.
\enumerate{
\item Individuals are drawn from a set of starting genotypes and randomly
distributed in a grid. Individuals are labelled by their genotype as 'g'
followed by an integer label.
\item Individuals are chosen at random to found the next generation. They
move seed dispersal to new location, with distances drawn from an exponential
distribution. The population exists on a torus to eliminat edge effects;
if dispersal takes a plant of the edge of the range it moves to the other
side of the range.
\item A subset of the new generation are chosen to have been germinated from
outcrossed seed. If so, they are given a new unique label by appending their
current genotype label by the generation number plus a unique integer.
\item At each generation a transect through the population is drawn at x=0 with
sampling points at equally spaced intervals. The plant closest to each
sampling point is recorded.
}

A seed bank is not currently implemented.

Plants exist within a box whose size is determined by the population size and
plant density (more plants at lower density = larger box). \code{sim_population}
will throw and error if the length of the transect
(\code{(n_sample_points-1) * sample_spacing}) is longer than the width of the box.
In this case, you can either make the population larger and less dense, or
choose fewer, closer sampling points.
}
\author{
Tom Ellis
}
