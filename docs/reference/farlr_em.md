# FARLR Approach II via EM

Implements FARLR Approach II using an EM algorithm. The method estimates
latent regression parameters from covariates with structured covariance
and supports iterative updating with convergence control.

## Usage

``` r
farlr_em(
  n,
  resp,
  parTab,
  K_hat,
  p,
  lambda_all,
  delta.criteria = 0.001,
  iter.max = 500,
  n_sam = 50,
  window.size = 50,
  theta_est_irt.mean,
  theta_est_irt.se,
  resp_rep,
  Z.em,
  main,
  verbose = TRUE
)
```

## Arguments

- n:

  Integer. Sample size.

- resp:

  Matrix. Item response matrix (individuals \\\times\\ items).

- K_hat:

  Integer. Estimated number of latent dimensions.

- p:

  Integer. Number of covariates.

- lambda_all:

  Numeric vector or matrix. Regularization parameters.

- delta.criteria:

  Numeric. Convergence threshold for parameter updates. Default is
  `1e-3`.

- iter.max:

  Integer. Maximum number of EM iterations. Default is `500`.

- n_sam:

  Integer. Number of Monte Carlo or posterior samples. Default is `50`.

- window.size:

  Integer. Window size used for convergence diagnostics. Default is
  `50`.

- theta_est_irt.mean:

  Numeric matrix. Initial latent ability means obtained from an IRT
  model.

- theta_est_irt.se:

  Numeric matrix. Standard errors of IRT-based latent ability estimates.

- resp_rep:

  Matrix. Replicated or augmented response matrix used in EM.

- Z.em:

  Matrix. Design matrix for EM-based latent regression.

- main:

  Logical. If `TRUE`, uses main effects only in the regression.

- a:

  Numeric vector or matrix. Item discrimination parameters.

- d:

  Numeric vector or matrix. Item difficulty (or intercept) parameters.

- c:

  Numeric vector. Guessing parameters (if applicable).

## Value

A list containing:

- `coef`: Estimated coefficients from the latent regression model.

- `sigma`: Estimated covariance matrix of the latent variables.

- `minBIC`: Minimum Bayesian Information Criterion (BIC) value achieved
  during model fitting.

- `converged`: Logical indicator of whether the EM algorithm converged.

- `cov`: Number of EM iterations required until convergence.

## Details

FARLR Approach II models covariate-dependent heteroskedasticity through
a low-rank random-effects formulation of the latent covariance
structure. Estimation is carried out via a Gaussian Variational EM
algorithm, where closed-form updates are available for latent regression
and covariance parameters, while convergence is monitored using a
sliding window criterion.
