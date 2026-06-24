# Principal Curves of Oriented Points

[![R-CMD-check](https://github.com/kmfrick/Rpcop/actions/workflows/r.yml/badge.svg)](https://github.com/kmfrick/Rpcop/actions/workflows/r.yml)
[![pkgcheck](https://github.com/kmfrick/Rpcop/actions/workflows/pkgcheck.yaml/badge.svg)](https://github.com/kmfrick/Rpcop/actions/workflows/pkgcheck.yaml)
[![Codecov test coverage](https://codecov.io/gh/kmfrick/Rpcop/graph/badge.svg)](https://app.codecov.io/gh/kmfrick/Rpcop)

This repository provides an implementation of principal curves on oriented points that can be compiled on modern systems, using Rcpp.

The original implementation, of which this repository is an adaptation, is provided [here](https://www-eio.upc.es/~delicado/PCOP/index.html) by the authors.

Principal curves on oriented points are introduced in Delicado and Huerta ([2003](https://link.springer.com/article/10.1007/s001800300145)).

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("kmfrick/Rpcop")
```

## Usage

```r
library(Rpcop)

set.seed(1)
n <- 120
t <- runif(n, -1, 1)
x <- cbind(t, t^2 + rnorm(n, sd = 0.08))

fit <- pcop(x, Ch = 1.5, Cd = 0.3, plot.true = FALSE)
names(fit)
summary(fit)
```

`pcop()` accepts a finite numeric matrix or numeric data frame. `Ch` must be a
single number between 0.5 and 1.5, and `Cd` must be a single number between 0.25
and 0.5. Missing, infinite, and non-numeric values are rejected before native
code is called.

PCOP is distance-based, so coordinate scaling can affect the fitted curve. The
wrapper first runs the native backend on the submitted scale. If that backend
fails, it retries on internally standardized coordinates and maps the returned
curve coordinates back to the submitted scale.

The returned object keeps the original `pcop.f1` and `pcop.f2` elements and also
supports `print()`, `summary()`, and `plot()` methods.

## Lifecycle and Prior Art

`Rpcop` is a maintained Rcpp port of the original PCOP C++ implementation
published by Delicado and Huerta. Development is focused on keeping the original
algorithm buildable on modern R toolchains, adding targeted validation, and
documenting package behavior for review rather than extending the method.

## Tests

Use these commands from repository root:

```bash
# Build the source tarball, then run a CRAN-style check outside the repo tree
rv sync
pkgdir=$PWD
libdir=$(rv library)
mkdir -p /private/tmp/rpcop-check
cd /private/tmp/rpcop-check
R_LIBS_USER="$libdir" R CMD build --no-manual "$pkgdir"
R_LIBS_USER="$libdir" R CMD check --as-cran --no-manual --output=/private/tmp/rpcop-check Rpcop_1.2.2.tar.gz
```

```bash
# Run Linux ASAN/UBSAN diagnostics
vagrant up --provision
```
