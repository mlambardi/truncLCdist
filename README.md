# truncLCdist
A non-uniform random number generation package in `R` for simulating extremely truncated log-concave distributions

This package implements an acceptance-rejection algorithm for sampling from a log-concave distribution and all its truncations to intervals in the form `]a,b]`.

There exist several `R` packages aimed at simulating truncated distributions, but often they are targeted at only few distributions. The package `truncLCdist` implements support for Tweedie distributions, like normal, binomial, Poisson, negative binomial, geometric, exponential, gamma, inverse-Gaussian. The user can add other distributions, provided they are **log-concave**. In fact, one can simulate any such distribution via Devroye's method (1986, 1987), just by providing:

* `f(x)` probability mass/density function, with its logarithmic version `f(x, log=T)=log f(x)`
* `F(x)=P(X≤x)` cumulative distribution function, with its logarithmic version `F(x, log.p=T)=log F(x)`
* `m=arg max f(x)` mode of the distribution

## Installation

It is possible to install `truncLCdist` via `devtools`, but there are few prerequisites depending on the operating system.

* **GNU/Linux**: The current version of `truncLCdist` requires two system libraries:
  * `gmp`, GNU Multiple Precision arithmetic library
  * `mpfr`, GNU Multiple Precision Floating-point Reliable library
  
  In future versions of this package, perhaps `gmp` and `mpfr` won't be required anymore. Fow now, if not installed, you can install them via the system console. On Ubuntu:
  ```
  sudo apt-get install libgmp-dev libmpfr-dev
  ```
  Tested on Ubuntu 22.04, `R 4.2`.

* **Windows**: The `R` package `devtools` works best if the `Rtools` utility is installed first. You can download `Rtools` from: https://cran.r-project.org/bin/windows/Rtools/, then install it manually. Alternatively, `RStudio` automatically proposes and carries out the installation of `Rtools` when a `C++`/`python`/`R markdown` file is opened.

  Tested on a Windows 11 virtual machine (`R 4.2`), see https://developer.microsoft.com/it-it/windows/downloads/virtual-machines/.

* **Mac OS X**: No additional library seems needed for installation on Mac OS X.

  Tested on Mac OS X Sierra, `R 4.2`.

After fixing these dependencies, you can install `devtools`, if missing, via the `R` console:
```
install.packages("devtools")
```
and then install `truncLCdist` via
```
devtools::install_github("mlambardi/truncLCdist")
```

Other `R`-only dependencies are installed automatically on some operating systems. If it fails, they can be installed manually via `install.packages()`, since they are available on `CRAN`.

## Usage

Try it out with:
```
truncLCdist::rtruncLC(n, spec, ..., a=-Inf, b=+Inf)
```
It will attempt to generate `n` variates distributed as `spec`, with optional parameters `...`, but truncated to the interval `[a,b]`. Try for instance
```
truncLCdist::rtruncLC(1000, "norm", mean=3, a=20)
```

Built-in support for few distributions and parameter definitions:

* `"binom"`: binomial, indexed by `size` (1, 2, ...) and `prob` (between 0 and 1)
* `"pois"`: Poisson, indexed by `lambda` (positive)
* `"nbinom"`: negative-binomial, indexed by `size` (1, 2, ...) and `prob` (between 0 and 1)
* `"geom"`: geometric, indexed by `prob` (between 0 and 1)
* `"gamma"`: gamma, indexed by `shape` (positive) and `rate` or `scale` (both positive)
* `"norm"`: normal, indexed by `mean` (real) and `sd` (positive)
* `"invgauss"`: inverse-Gaussian, indexed by `mean` (positive) and `shape` or `dispersion` (positive)

The binomial distribution was implemented by means of the command
```
truncLCdist::truncLCdistoptions(
  nbinom.continuous=F, ## is it a continuous distribution?
  nbinom.m=function(size, prob) pmax(floor((size - 1)*(1 - prob)/prob), 0), ## the mode of the distribution
  nbinom.d=Rmpfr::dnbinom, ## the probability mass/density function
  nbinom.p=pnbinom ## the cumulative distribution function
  )
```
Custom distributions can be added by providing the required `.continuous`, `.m`, `.d` and `.p` functions.

* The `.d` function requires a logarithmic version that activates with the option `log=T`.
* The `.p` function requires a logarithmic version that activates with the option `log.p=T`.
* All the functions must be parameterized consistently.

## References

* Devroye, L. (1986) *Non-Uniform Random Variate Generation*. Springer-Verlag, New York. https://doi.org/10.1007/978-1-4613-8643-8
* Devroye, L. (1987) *A simple generator for discrete log-concave distributions*. Computing **39**, 87--91. https://doi.org/10.1007/BF02307716

Also see our working paper for the rationale:

* Lambardi di San Miniato, M., Kenne-Pagui, E.C. (2022) *Scalable random number generation for truncated log-concave distributions*. arXiv preprints. https://doi.org/10.48550/arXiv.2204.01364
