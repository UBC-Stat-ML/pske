# pske: Probabilistic Solutions to Kolmogorov's Equations

## Description

`pske` (pronounced "peh·skee") is a collection of probabilistic solvers of Kolmogorov's equations for Markov jump processes.

## Installation

``` r
if (!require(remotes)) {
    install.packages('remotes')
}
remotes::install_github("UBC-Stat-ML/pske")
```

## TODO

- Basic parameter checking to avoid crashing R session by passing wrong things to C code
- More solvers
