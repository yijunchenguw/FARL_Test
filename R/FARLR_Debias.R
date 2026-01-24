library(torch)
#' FARLR Debiased Estimation for Regularized Latent Regression
#'
#' Fits a Factor-Augmented Regularized Latent Regression (FARLR) model with a
#' post-selection debiasing procedure. Latent traits are approximated using a
#' Monte Carlo scheme based on a normal approximation to IRT posterior estimates
#' (\code{theta_est_irt.mean}, \code{theta_est_irt.se}). Regression coefficients
#' are first obtained via weighted LASSO regularization using \code{glmnet}, after
#' which a debiased refit is performed on the selected active set using weighted
#' least squares. To improve numerical stability, coefficient updates are smoothed
#' across iterations using a sliding-window averaging scheme. The regularization
#' parameter \code{lambda} is selected by minimizing a BIC-type criterion over
#' \code{lambda_all}.
#'
#' @param n Integer. Number of individuals (sample size).
#' @param resp Matrix. Observed item response matrix of dimension \code{n x J}.
#' @param parTab Data frame. Item parameter table containing at least the columns
#'   \code{slope}, \code{difficulty}, and \code{guessin}.
#' @param K_hat Integer. Number of latent components included in the regression
#'   design matrix \code{Z.em}.
#' @param p Integer. Number of observed covariates included in \code{Z.em}.
#' @param lambda_all Numeric vector. Candidate regularization parameters supplied
#'   to \code{glmnet}.
#' @param delta.criteria Numeric. Convergence tolerance for iterative updates.
#'   Defaults to \code{1e-3}.
#' @param iter.max Integer. Maximum number of iterations for each value of
#'   \code{lambda}. Defaults to \code{500}.
#' @param n_sam Integer. Number of Monte Carlo samples per individual used to
#'   approximate latent trait uncertainty. Defaults to \code{5}.
#' @param window.size Integer. Window size for sliding-window averaging of
#'   regression coefficient updates. Defaults to \code{50}.
#' @param theta_est_irt.mean Numeric vector of length \code{n}. Posterior mean
#'   estimates of the latent trait obtained from an IRT model.
#' @param theta_est_irt.se Numeric vector of length \code{n}. Posterior standard
#'   error estimates of the latent trait obtained from an IRT model.
#' @param resp_rep Matrix. Replicated or expanded response matrix used for Monte
#'   Carlo integration (see \code{q_num_NA}).
#' @param Z.em Matrix. Regression design matrix, typically of dimension
#'   \code{n x (K_hat + p)}.
#' @param Uupdate Internal object. Passed to internal update routines controlling
#'   coefficient smoothing.
#' @param hatU Internal object. Passed to internal routines used in the debiasing
#'   refit.
#' @param Fan Internal object. Passed to internal routines for adaptive weighting
#'   or regularization.
#' @param main Integer vector. Indices of predictors to be treated as unpenalized
#'   in \code{glmnet} via \code{penalty.factor = 0}. Remaining predictors are
#'   penalized unless they correspond to the first \code{K_hat} latent components.
#' @param verbose Logical. If \code{TRUE}, progress messages and a progress
#'   indicator are displayed. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{coefficients}}{Debiased regression coefficient estimates
#'   corresponding to the selected \code{lambda} (length \code{K_hat + p}).}
#'   \item{\code{sigma}}{Estimated residual standard deviation.}
#'   \item{\code{LogLik}}{Value of the BIC-type objective function at the selected
#'   \code{lambda}. Returned as \code{LogLik} for compatibility.}
#'   \item{\code{minBIC}}{Index of \code{lambda_all} that minimizes the BIC-type
#'   criterion.}
#'   \item{\code{Convergence}}{Character string indicating convergence status of
#'   the iterative procedure.}
#' }
#'
#' @details
#' For each candidate value in \code{lambda_all}, coefficient estimates are updated
#' iteratively until convergence or until \code{iter.max} iterations are reached.
#' Convergence is assessed using the maximum absolute change in regression
#' coefficients and the residual scale parameter:
#' \deqn{
#' \delta = \max\left(
#' \lvert \sigma^{(t)} - \sigma^{(t-1)} \rvert,
#' \max_j \lvert \beta_j^{(t)} - \beta_j^{(t-1)} \rvert
#' \right).
#' }
#'
#' The debiasing step refits the regression model on the active set selected by the
#' penalized estimator using weighted least squares, reducing shrinkage bias and
#' improving finite-sample interpretability of the coefficient estimates.
#'
#' This function depends on auxiliary routines (not shown here), including
#' \code{q_num_NA()} for Monte Carlo integration and \code{add_to_window()} for
#' sliding-window averaging.
#'
#' @seealso \code{\link[glmnet]{glmnet}}, \code{\link[mirt]{simdata}}
#'
#' @examples
#' \dontrun{
#' fit <- farlr_debias(
#'   n = nrow(resp),
#'   resp = resp,
#'   parTab = parTab,
#'   K_hat = K_hat,
#'   p = p,
#'   lambda_all = seq(0.001, 0.1, length.out = 10),
#'   theta_est_irt.mean = theta_mean,
#'   theta_est_irt.se = theta_se,
#'   resp_rep = resp_rep,
#'   Z.em = Z.em,
#'   main = 1,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
farlr_debias <- function(n,resp, parTab, K_hat, p, lambda_all, delta.criteria = 1e-3, iter.max = 500, n_sam = 5, window.size = 50, theta_est_irt.mean, theta_est_irt.se,resp_rep = NA, Z.em = NA, Uupdate,hatU, Fan, main, verbose = TRUE)
{
  # method debias
  results <- list()
  a <- parTab$slope
  d <- -parTab$slope * parTab$difficulty
  c <- parTab$guessing
  output <- if (verbose) stderr() else nullfile()
  cat(file = output, 'Fitting the latent regression model...\n')
  Uupdate.em <- rep(1, n_sam) %x% Uupdate
  Uupdate.em.torch <- torch_tensor(Uupdate.em)
  resp_rep <- rep(1, n_sam) %x% resp
  hatU.em <- rep(1, n_sam) %x% hatU
  hatU.em.torch <- torch_tensor(hatU.em)
  Fan.em <- rep(1, n_sam) %x% Fan
  temp <- Uupdate.em.torch$matmul(Uupdate.em.torch$t())
  P <- (1 / (n * n_sam)) * temp
  #P <- 1/(n*n_sam)*torch_matmul(Uupdate.em.torch,Uupdate.em.torch$t())
  Z.em <- rep(1, n_sam) %x% Fan
  if (verbose) pb <- txtProgressBar(file = output, 0, length(lambda_all), style = 3)
  for (ll in 1:length(lambda_all)) {
    # initial value
    #if (verbose) cat("\n", file = output)
    beta_gamma_old <- rep(0,K_hat+p)
    sigma_old <- 1
    beta_gamma_t <- list()
    sigma_t <- list()
    window_size = window.size

    delta <- 1
    iter <- 1
    theta_t <- 1
    cov_t <- c()
    cv.fit.em <- list()
    cv.fit.em.debias <- list()
    lambda <- lambda_all[ll]

    while ((delta > delta.criteria) && (iter < iter.max)) {
      # if (verbose && iter %% 10 == 0) {
      #   cat(".", file = output)
      #   if (iter %% 30 == 0) cat("\n", file = output)
      # }
      # E-step : sample the theta
      theta_sample <- rnorm(n*n_sam, theta_est_irt.mean, (theta_est_irt.se + 0.2))
      theta_sample.torch <- torch_tensor(theta_sample)
      q_num_sample <- q_num_NA(a,d, c, theta_sample,resp_rep,Z.em,beta_gamma_old,sigma_old)

      h_sample <- dnorm(theta_sample, theta_est_irt.mean, (theta_est_irt.se + 0.2))

      den_all <- q_num_sample/h_sample

      den_i <- (1/n_sam)*as.numeric(tapply(den_all, (seq_along(den_all) - 1) %% nrow(resp) + 1, sum))

      w_ik <- (1/den_i)*q_num_sample/h_sample
      w_ik.torch <- torch_tensor(w_ik)
      theta_1 <-theta_sample.torch-torch_matmul(P,theta_sample.torch)
      factors <- rep(1,p)
      factors[main] <- 0
      cv.fit.em[[iter]] <- glmnet(hatU.em, as.matrix(theta_1),
                                  weights = w_ik, penalty.factor=factors,intercept = FALSE, lambda = lambda)
      coef_hat_em <- coef(cv.fit.em[[iter]],complete=TRUE)[-1]
      coef_hat_em.torch <- torch_tensor(coef_hat_em)
      weights <- torch_diag(w_ik.torch)
      C<-torch_diag(rep(1,p))
      # T<-c()
      # for( j in 1:ncol(hatU.em)){
      #
      #   C[j,-j]<-0
      #   T<-c(T,1/(n*n_sam)*sum(w_ik*(hatU.em[,j])^2))
      #   # if (j%%100==0){
      #   #   print(j)
      #   # }
      # }
      # Ensure w_ik is [n, 1] for broadcasting
      w_ik_exp <- w_ik.torch$unsqueeze(2)  # [n, 1]

      # Weighted square of hatU_em: [n, p]
      weighted_sq <- w_ik_exp * hatU.em.torch$pow(2)

      # Sum across rows (i.e., for each column j): result [p]
      col_sums <- weighted_sq$sum(dim = 1)

      # Compute T vector
      T <- (1 / (n * n_sam)) * col_sums
      # Create diagonal matrix C: [p, p]
      T1<-torch_diag(1/T)
      Theta1<-torch_matmul(T1,C)
      temp <- (1/(n*n_sam))* torch_matmul(Theta1,hatU.em.torch$t())
      residuals <- theta_1 - torch_matmul(hatU.em.torch, coef_hat_em.torch)
      coef_hat_debias <- coef_hat_em + torch_matmul(torch_matmul(temp, weights), residuals)
      coef_hat_debias <- as.array(coef_hat_debias)
      #debias ----
      Ut <- Uupdate.em.torch$transpose(1, 2)

      # Full left-hand matrix: (Uᵗ W U)
      lhs <- Ut$matmul(weights)$matmul(Uupdate.em.torch)

      # Right-hand side: Uᵗ W θ
      rhs <- Ut$matmul(weights)$matmul(theta_sample)
      #if (verbose) cat(".", file = output)
      # Solve the system (lhs)^(-1) * rhs
      phi_hat <- torch_inverse(lhs)$matmul(rhs)
      #phi_hat <-(solve(t(Uupdate.em)%*%Uupdate.em))%*%t(Uupdate.em) %*% theta_sample
      phi_hat <- as.array(phi_hat)
      coef_hat_em_debias <- ifelse(coef_hat_em !=0, (coef_hat_debias), 0)
      coef_hat_em_debias <- append(coef_hat_em_debias, (phi_hat), after = 0)

      #add into sliding window
      beta_gamma_t <- add_to_window(coef_hat_em_debias, beta_gamma_t, window_size)
      #calculate mean
      matrix_beta_gamma <- do.call(cbind, beta_gamma_t)

      beta_gamma_means <- rowMeans(matrix_beta_gamma)

      fitted_values.em <- matrix(beta_gamma_means, nrow = 1)%*%t(Z.em)
      residuals_2 <- theta_sample - fitted_values.em
      WSSR <- sum(w_ik*(residuals_2^2))#
      sigma_temp <- sqrt(WSSR/(nrow(Z.em)-p))

      matrix_sigma <- do.call(cbind, sigma_t)
      sigma_means <- sigma_temp #rowMeans(matrix_sigma)

      #joint response ???hg
      bic <- -2*(-n/2*log(2*pi*sigma_temp^2)- 1/(2*n_sam*sigma_temp^2) * sum(w_ik*(theta_sample - Z.em%*%beta_gamma_means)^2)) + sum(beta_gamma_means!=0)*log(n)

      #calculate delta
      delta_s <- abs(sigma_old-sigma_means)
      delta_b <- max(abs(beta_gamma_old - beta_gamma_means))
      delta <- max(delta_s, delta_b)

      sigma_old <- sigma_means
      beta_gamma_old <- beta_gamma_means
      iter <- iter + 1
      #print(delta)
    }
    #print(ll)
    #print(iter)
    if (verbose) {setTxtProgressBar(pb, pb$getVal() + 1)}
    results[[ll]] <- list(beta_gamma_old,sigma_old,bic)
  }

  minBIC <- which.min(unlist(lapply(results, function(x) x[[3]])))
  beta_hat <- results[[minBIC]][[1]]; sigma <- results[[minBIC]][[2]]
  n_sam <- 100
  Uupdate.em <- rep(1, n_sam) %x% Uupdate
  resp_rep <- rep(1, n_sam) %x% resp
  hatU.em <- rep(1, n_sam) %x% hatU
  Z.em <- rep(1, n_sam) %x% Fan

  theta_sample <- rnorm(n*n_sam, theta_est_irt.mean, (theta_est_irt.se + 0.2))

  q_num_sample <- q_num_NA(a,d,c, theta_sample,resp_rep,Z.em,beta_gamma_old,sigma_old)

  h_sample <- dnorm(theta_sample, theta_est_irt.mean, (theta_est_irt.se + 0.2))

  den_all <- q_num_sample/h_sample

  den_i <- (1/n_sam)*as.numeric(tapply(den_all, (seq_along(den_all) - 1) %% nrow(resp) + 1, sum))

  w_ik <- (1/den_i)*q_num_sample/h_sample

  fitted_values.em <- matrix(beta_hat, nrow = 1)%*%t(Z.em)
  residuals_2 <- theta_sample - fitted_values.em
  WSSR <- sum(w_ik*(residuals_2^2))#
  sigma_temp <- sqrt(WSSR/(nrow(Z.em)-p))

  return(list(
    coefficients  = beta_hat,
    sigma = sigma_temp,
    LogLik = results[[minBIC]][[3]],
    minBIC = minBIC,
    Convergence = ifelse(iter<iter.max,"Converged","Did not converge")
  ))
}
