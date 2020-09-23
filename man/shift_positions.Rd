% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift_positions.R
\name{shift_positions}
\alias{shift_positions}
\title{Shift a matrix of x- and y-coordinates in random directions}
\usage{
shift_positions(coords, mean_dispersal_distance, range_limit)
}
\arguments{
\item{coords}{Matrix with a row for each individual and two rows.}

\item{mean_dispersal_distance}{Float giving the average perturbation
distance. The reciprocal of this value is used as a rate for the exponential
distribution.}

\item{range_limit}{Float giving the diameter of the torus.}
}
\value{
A matrix matching the shape of the input \code{coords}, but with values
shifted by some exponential distance.
}
\description{
\code{shift_positions} peturbs the positions of a matrix of positions by
distances drawn from an exponential distribution. Dispersal is simulated on
a torus, so that if positions are peturbed beyond a pre-defined range limit
they are reflected to the other side of the population.
}