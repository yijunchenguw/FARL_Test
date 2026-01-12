# Fit a Latent Regression Model via the EM Algorithm

Fits a latent regression model using an expectation–maximization (EM)
algorithm, where latent abilities are linked to observed covariates and
item response data through an IRT measurement model.

## Usage

``` r
mml(X, Y, parTab, p, n_sam = 30, verbose = TRUE)
```

## Arguments

- X:

  Numeric matrix. Design matrix of covariates (\\N \times p\\).

- Y:

  Matrix. Item response matrix (\\N \times J\\).

- p:

  Integer. Number of covariates included in the latent regression.

- n_sam:

  Integer. Number of Monte Carlo or quadrature samples used to
  approximate posterior expectations in the E-step. Default is `30`.

- verbose:

  Logical. If `TRUE`, prints progress information during estimation.
  Default is `TRUE`.

- a:

  Numeric vector or matrix. Item discrimination parameters.

- d:

  Numeric vector or matrix. Item difficulty (intercept) parameters.

## Value

A list containing:

- `coef`: Estimated coefficients from the latent regression model.

- `Sigma`: Estimated covariance matrix of the latent variables.

- `theta`: Posterior mean estimates of latent abilities.

- `converged`: Logical indicator of whether the EM algorithm converged.

- `iter`: Number of EM iterations performed.

## Details

The EM algorithm alternates between computing the conditional
expectations of latent abilities given the observed responses and
current parameter estimates (E-step) and maximizing the expected
complete-data log-likelihood with respect to the model parameters
(M-step). Numerical optimization is employed for item parameters, while
closed-form updates are available for latent regression and covariance
parameters.
