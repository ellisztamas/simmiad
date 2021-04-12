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

An R package for simulating a population of wild Emmer wheat to ask whether the amount of spatial clustering of unique genotypes and the stability of that clustering through time can be explained by purely neutral forces. The idea is to simulate a population of plants evolving under seed dispersal and limited, random outcrossing only, then to sample plants along a transect in the same way that the real population is sampled.

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

## How simulations work

### Initial generation:

* This simulates a population of plants at a given density on a torus (i.e. there are no edges). The radius of the torus is 50% longer than the length of the transect
* Population size is determined as the number of plants needed to fill the torus given its area and population density.
* Plants are initially distributed at random throughout the habitat and assigned one of a set of unique genotypes.

### Simulating through time:

* In the next generation, each plant generates an average of one seed. Offspring numbers are drawn from a multinomial distribution of size N, with probablity 1/N for each mother, where N is population size.
* Each seed disperses in a random direction at a distance drawn from an exponential distribution.
* Each seed has some probability of having been the product of an outcrossing event. If so, it is assigned a new unique genotype. If not, it is assumed to have been selfed and shares the genotype of its mother.
* This is repeated for many generations (I used 500 generations)

### Transect samples

* At the end of the simulation, plants are sampled along a transect over multiple generations.
* A transect is drawn through the middle of the population with evenly spaced sampling points.
* At each sampling point, we sample the plant closest to the sampling point. If no plant is within 1m of the sampling point, then no plant is recorded.
* This is repeated for some number of years back into the past from the final generation.

### Assumptions
This makes certain assumptions that it is good to be explicit about:

* There are no differences in fitness between genotypes, or genotype-by-environment interactions for fitness with a heterogeneous landscape.
* Population density is even across the landscape.
* Seed dispersal distances are exponentially distributed. I am hoping dispersal is primarily through gravity and is fairly short scale. If there is something more complicated happening, for example additional longer-distance dispersal by rodents, this could be modelled with some kind of mixture of distributions, This would complicate things.
* Outcrossing is random. In reality there will be some kind of pollen dispersal kernel shape, but I have no idea how that should look here.
* Seed dispersal and outcrossing rates/distances do not change through time.

## Usage
### Simulate a single population

Functions in `simmiad` simulate populations given a set of input parameters:

* **mean_dispersal_distance** Mean seed dispersal distance in metres.
* **outcrossing_rate** Probability that an individual is outcrossed.
* **n_generations** Number of generations to run the simulations.
* **n_starting_genotypes** Number of initial genotypes to start with.
* **density** Average density of plants per square metre.

To simulate a single population you can use `sim_population`. For example, this runs a single simulation of a population with 3x3=9 plants of 124 genotypes for 100 generations, with mean dispersal distance of 3m and an outcrossing rate of 1%.

```
library('simmiad')
set.seed(124) # so you get the same answer as me

# Set input parameters
mean_dispersal_distance = 0.5
outcrossing_rate = 0.01
n_generations = 10
n_starting_genotypes = 10
density = 1
how_far_back <- n_generations
n_sample_points = 30
sample_spacing = 5

sm <- sim_population(
  mean_dispersal_distance = mean_dispersal_distance,
  outcrossing_rate = outcrossing_rate,
  n_generations = n_generations,
  n_starting_genotypes = n_starting_genotypes,
  density = density,
  n_sample_points = n_sample_points,
  sample_spacing = sample_spacing,
  )
```

This returns a list of genotypes in each generation. The final generation looks like this:

```
sm[[9]]
 [1] "g5"       "g4"       "g7"       "g6"       NA         "g3"       "g9"      
 [8] "g6"       "g2"       NA         "g2"       "g5"       "g2"       "g3"      
[15] "g2"       "g5"       "g4"       "g4"       "g7"       "g4"       "g9"      
[22] "g10"      NA         "g5"       "g1"       "g1"       "g10"      "g9"      
[29] "g8_6.377" "g10"          
```
'g' stands for genotype, and is followed by a number between 1 and 10 indicating the id of the initial genotype. Individual 29 shows what happens if outcrossing occurs: the genotype label is appended by the generation outcrossing occured and a unique integer within that generation. That ensures every outcrossed genotype is a new unique label. Note that if outcrossed genotypes outcross again the names will keep getting longer (and messier!).

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
Most of the time you will want to simulate multiple replicate populations with a set of input parameters. This can be done with the function `simmiad` using similar input parameters as [before](#Simulate-a-single-population)

```
rs <- simmiad(
  mean_dispersal_distance = 0.5,
  outcrossing_rate = 0.001,
  n_generations = 12,
  n_starting_genotypes = 10,
  density = 3,
  n_sample_points = 5,
  sample_spacing = 2,
  nsims = 3
)
```
This function simulates multiple individual populations through time, and returns a list of different data:

1. **parameters** A data.frame giving input parameters.
2. **clustering** The covariance between distance along the transect and the frequency of identical genotypes.
3. 

## Author and license information

Tom Ellis (thomas.ellis@gmi.oeaw.ac.at)

`simmiad` is available under the MIT license. See LICENSE for more information.
