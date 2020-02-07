# simmiad
R package to simulate populations of wild Emmer wheat from the Kibbutz Ammiad

## Table of contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Dependencies](#dependencies)
4. [Usage](#usage)
    1. [Simulate a single population](#simulate-a-single-population)
    2. [Describing spatial structure](#describing-spatial-structure)
    3. [Replicate simulations](#replicate-simulations)
5. [Author and license information](#author-and-license-information)

## Introduction

An R package for simulating a population of wild Emmer wheat to ask whether the amount of spatial clustering of unique genotypes can be explained by purely neutral forces.

This simulates a population on an evenly-spaced grid of sampling points through time:

* Sampling points are populated at random with a set of unqiue starting genotypes.
* In the next generation, each sampling point is filled by one seed from the same or a different sampling point from a previous generation. Whichever genotype was at the donor sampling point now occupies the focal sampling point. I assume seed dispersal distances are exponentially distributed with some mean value that is used as an input parameter.
* To allow for a seed bank, each seed can be drawn from any previous generation. The generation is drawn from a poisson distribution.
* Each seed has some probability of having been the product of an outcrossing event. If so, it is assigned a new unique genotype.

This makes certain assumptions that it is good to be explicit about:

* There are no differences in fitness between genotypes, or genotype-by-environment interactions for fitness with a heterogeneous landscape.
* Population density is even across the landscape.
* Seed dispersal distances are exponentially distributed. I am hoping dispersal is primarily through gravity and is fairly short scale. If there is something more complicated happening, for example additional longer-distance dispersal by rodents, this could be modelled with some kind of mixture of distributions, This would complicate things.
* Outcrossing is random. In reality there will be some kind of pollen dispersal kernel shape, but I have no idea how that should look here.
* Seed dispersal and outcrossing rates/distances do not change through time.
* Plants live in a grid. They don't...

## Installation

Installation is easiest straight from GitHub using the package `devtools` from within R.
If necessary, install this with

```
install.packages("devtools")
```

Then you can install with

```
devtools::install_github("ellisztamas/simmiad")
```

## Dependencies

`simmiad` uses base R functions only.

## Usage
### Simulate a single population

Functions in `simmiad` simulate populations given a set of input parameters:

* **grid_size**: Number of rows in the grid. Population size will be the square of this number.
* **mean_dispersal_distance** Mean seed dispersal distance in metres.
* **outcrossing_rate** Probability that an individual is outcrossed.
* **n_generations** Number of generations to run the simulations.
* **n_starting_genotypes** Number of initial genotypes to start with.

To simulate a single population you can use `sim_population`. For example, this runs a single simulation of a population with 3x3=9 plants of 124 genotypes for 100 generations, with mean dispersal distance of 3m and an outcrossing rate of 1%.

```
library('simmiad')
set.seed(124) # so you get the same answer as me

sm <- sim_population(
  grid_size = 3,
  mean_dispersal_distance = 3,
  outcrossing_rate = 0.01,
  n_generations = 100,
  n_starting_genotypes = 124
  )
```

This returns a matrix of the genotypes of the population in the final generation. You can set it to return genotypes from every generation by setting `return_all=TRUE`, but I don't recommend printing that to the console for larger simulations. The final generation looks like this:

```
[,1]      [,2]      [,3]
[1,] "g5"      "g45_8.1" "g5"
[2,] "g45_8.1" "g5"      "g5"
[3,] "g5"      "g45"     "g5"
```
'g' stands for genotype, and is followed by a number between 1 and 124 indicating the id of the initial genotype. If an outcrossing event has occurred the genotype label is appended by the generation outcrossing occured (8) and a unique integer within that generation. That ensures every outcrossed genotype is a new unique label. As this is the last of 10 generations, and the outcrossing event occured in generation 8, we can see that the outcrossed individual has left two offspring two generations later. Note that if outcrossed genotypes outcross again the names will keep getting longer (and messier!).

### Describing spatial structure
I am not sure this is the best way to acheive this, but one simple measure of clustering of unique genotypes is to compare the average distance between pairs of identical genotypes to pairs of non-identical genotypes. `transect_clustering` will do this for samples along a single transect.

```
# Make a larger population
sm2 <- sim_population(
  grid_size = 100,
  mean_dispersal_distance = 3,
  outcrossing_rate = 0.01,
  n_generations = 30,
  n_starting_genotypes = 126
)

transect_clustering(
  genotypes = sm2[4,], # Genotypes along the 4th row in the population
  positions = 1:100 # Spaial positions along that transect
  )
```

### Replicate simulations
Most of the time you will want to simulate multiple replicate populations with a set of input parameters. This can be done with the function `simmiad` using the same input parameters as [before](#Simulate-a-single-population), plus a file path to output the results and a value for the number of simulations to run:

```
simmiad(
  grid_size = 100,
  mean_dispersal_distance = 3,
  outcrossing_rate = 0.01,
  n_generations = 100,
  n_starting_genotypes = 124,
  filename='test_run',
  nsims = 3
  )
```
This function simulates individual populations through time, then takes a single horizontal transect through the final generation. It then estimates [spatial clustering](#describing-spatial-structure) through that transect.

This will not create an object in R, but will output the results directly to a CSV file called 'test_run.csv', and also save a log file called 'test_run.log'. You can change the path and file names by changing the `filename` argument. Open the results with
```
sm <- read.csv('test_run.csv')
```

Ther are three columns: the first column gives the simulation number; the second and third give the output of `transect_clustering`.

## Author and license information

Tom Ellis (thomas.ellis@gmi.oeaw.ac.at)

`simmiad` is available under the MIT license. See LICENSE for more information.
