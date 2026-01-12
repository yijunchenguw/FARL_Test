# FARL: Use Factor-augmented Regularized Latent Regression for Large Scale Assessment ![](reference/figures/logo.png)

The goal of `FARL` is to provide computationally efficient tools for
large scale assessment, to generate Plausible Values, with high
accuracy. The package contains several example datasets and functions
for

- FARLR: Factor-augmented Regularized Latent Regression
- FARLR-debias: Factor-augmented Regularized Latent Regression - Debias

## Installation

To install this package from source:

1.  **Windows** users may need to install the
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) and include
    the checkbox option of installing Rtools to their path for easier
    command line usage. **Mac** users may need to install Xcode command
    line tools by `sudo xcode-select --install` in the terminal, and
    then install [GNU Fortran
    compiler](https://mac.r-project.org/tools/). Most **Linux**
    distributions should already have up to date compilers (or if not
    they can be installed/updated easily).

2.  Install the `devtools` package (if necessary), and install the
    package from [GitHub](https://github.com/) with

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("MAP-LAB-UW/FARL", build_vignettes = T)
torch::install_torch()
```

## Tutorial

After installing the `FARL` package, open its tutorial by running

``` r
vignette("FARL")
```
