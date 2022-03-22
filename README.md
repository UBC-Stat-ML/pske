# pske: Probabilistic Solutions to Kolmogorov's Equations

## Description

`pske` (pronounced "peh·skee") is a collection of probabilistic solvers of Kolmogorov's equations for Markov jump processes.

## Installation

### Step 1: system dependencies

Building the underlying C/C++ code requires a working C++11 toolchain. Linux systems usually ship with such tools, but macOS and Windows users might require additional setup:

- Windows: install the [Rtools toolchain](https://cran.r-project.org/bin/windows/Rtools/).
- macOS: Follow [these instructions](https://github.com/rmacoslib/r-macos-rtools#how-do-i-use-the-installer) to configure the C toolchain.


### Step 2: install `pske`

``` r
if (!require(remotes)) {
    install.packages('remotes')
}
remotes::install_github("UBC-Stat-ML/pske")
```

## TODO

- Basic parameter checking to avoid crashing R session by passing wrong things to C code
- More solvers
