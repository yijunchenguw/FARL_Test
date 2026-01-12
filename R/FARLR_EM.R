q_num_NA <- function(a, d, c, theta, resp, Z, beta_gamma_t, sigma_t) {
  # pij
  pij <- c + (1 - c) * (1 / (1 + exp(-(as.matrix(theta) %*% a + rep(d, each = length(theta))))))

  # mask non NA
  mask <- !is.na(resp)

  lpij_irt <- matrix(0, nrow = nrow(pij), ncol = ncol(pij))  # 0
  lpij_irt[mask] <- resp[mask] * log(pij[mask]) + (1 - resp[mask]) * log(1 - pij[mask])

  pi_irt <- exp(rowSums(lpij_irt))  # summing over items

  p_reg <- as.vector(dnorm(theta, beta_gamma_t %*% t(Z), sigma_t))
  num <- p_reg * pi_irt
  return(num)
}
# Function to add new results to the sliding window
add_to_window <- function(new_result, results, window_size) {
  # Add the new result
  results <- c(results, list(new_result))

  # If the window exceeds the size, remove the oldest result
  if (length(results) > window_size) {
    results <- results[-1]  # Remove the first element
  }

  return(results)
}
#' FARLR Approach II via EM
#'
#' Implements FARLR Approach II using an EM algorithm. The method estimates latent
#' regression parameters from covariates with structured
#' covariance and supports iterative updating with convergence control.
#'
#' @param n Integer. Sample size.
#' @param resp Matrix. Item response matrix (individuals \eqn{\times} items).
#' @param a Numeric vector or matrix. Item discrimination parameters.
#' @param d Numeric vector or matrix. Item difficulty (or intercept) parameters.
#' @param c Numeric vector. Guessing parameters (if applicable).
#' @param K_hat Integer. Estimated number of latent dimensions.
#' @param p Integer. Number of covariates.
#' @param lambda_all Numeric vector or matrix. Regularization parameters.
#' @param delta.criteria Numeric. Convergence threshold for parameter updates.
#'   Default is \code{1e-3}.
#' @param iter.max Integer. Maximum number of EM iterations.
#'   Default is \code{500}.
#' @param n_sam Integer. Number of Monte Carlo or posterior samples.
#'   Default is \code{50}.
#' @param window.size Integer. Window size used for convergence diagnostics.
#'   Default is \code{50}.
#' @param theta_est_irt.mean Numeric matrix. Initial latent ability means
#'   obtained from an IRT model.
#' @param theta_est_irt.se Numeric matrix. Standard errors of IRT-based latent
#'   ability estimates.
#' @param resp_rep Matrix. Replicated or augmented response matrix used in EM.
#' @param Z.em Matrix. Design matrix for EM-based latent regression.
#' @param main Logical. If \code{TRUE}, uses main effects only in the regression.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coef}: Estimated coefficients from the latent regression model.
#'   \item \code{sigma}: Estimated covariance matrix of the latent variables.
#'   \item \code{minBIC}: Minimum Bayesian Information Criterion (BIC) value achieved during model fitting.
#'   \item \code{converged}: Logical indicator of whether the EM algorithm converged.
#'   \item \code{cov}: Number of EM iterations required until convergence.
#' }
#'
#' @details
#' FARLR Approach II models covariate-dependent heteroskedasticity through
#' a low-rank random-effects formulation of the latent covariance structure.
#' Estimation is carried out via a Gaussian Variational EM algorithm, where
#' closed-form updates are available for latent regression and covariance
#' parameters, while convergence is monitored using a sliding window criterion.
#'
#'
#' @export
farlr_em <- function(n,resp, parTab, K_hat, p, lambda_all, delta.criteria = 1e-3, iter.max = 500, n_sam = 50,  window.size = 50, theta_est_irt.mean, theta_est_irt.se, resp_rep, Z.em, main, verbose = TRUE){
  set.seed(234)
  results <- list()
  a <- parTab$slope
  d <- -parTab$slope * parTab$difficulty
  c <- parTab$guessin
  output <- if (verbose) stderr() else nullfile()
  cat(file = output, 'Fitting the latent regression model...\n')
  if (verbose) pb <- txtProgressBar(file = output, 0, length(lambda_all), style = 3)
  for (ll in 1:length(lambda_all)) {
    # initial value
    if (verbose) setTxtProgressBar(pb, pb$getVal() + 1)
    beta_gamma_old <- rep(0,K_hat+p)
    beta_gamma_t <- c()
    sigma_old <- 1

    delta <- 1
    iter <- 1
    theta_t <- 1
    cov_t <- c()
    cv.fit.em <- list()
    cv.fit.em.debias <- list()
    lambda <- lambda_all[ll]
    window_size <- window.size

    while ((delta > delta.criteria) && (iter < iter.max)) {
      # E-step : sample the theta

      theta_sample <- rnorm(n*n_sam, theta_est_irt.mean, (theta_est_irt.se + 0.2))

      q_num_sample <- q_num_NA(a,d,c, theta_sample,resp_rep,Z.em,beta_gamma_old,sigma_old)

      h_sample <- dnorm(theta_sample, theta_est_irt.mean, (theta_est_irt.se + 0.2))

      den_all <- q_num_sample/h_sample

      den_i <- (1/n_sam)*as.numeric(tapply(den_all, (seq_along(den_all) - 1) %% nrow(resp) + 1, sum))
      w_ik <- (1/den_i)*q_num_sample/h_sample

      factors <- rep(1,p)
      factors[main] <- 0
      cv.fit.em[[iter]] <- glmnet(Z.em, theta_sample, standardize=FALSE, family="gaussian",
                                  penalty.factor=c(rep(0,K_hat),factors), weights = w_ik, intercept = FALSE, lambda = lambda)
      #plot(cv.fit.em)
      coef_hat_em <- coef(cv.fit.em[[iter]],complete=TRUE)[-1]

      #lm ----
      cv.fit.em.debias[[iter]] <- lm(theta_sample~0+Z.em[,which(coef_hat_em!=0)], weights = w_ik)
      coef_hat_em_debias <- rep(0,p+K_hat)
      coef_hat_em_debias[which(coef_hat_em!=0)] <- coef(cv.fit.em.debias[[iter]])
      coef_hat_em_debias <- ifelse(is.na(coef_hat_em_debias), 0 ,coef_hat_em_debias )
      # cross validation
      # cv.fit
      beta_gamma_t <- add_to_window(coef_hat_em_debias, beta_gamma_t, window_size)

      #calculate mean
      matrix_beta_gamma <- do.call(cbind, beta_gamma_t)

      beta_gamma_means <- rowMeans(matrix_beta_gamma)

      cov <- vcov(cv.fit.em.debias[[iter]])
      fitted_values.em <- matrix(beta_gamma_means, nrow = 1, ncol = p+K_hat)%*%t(Z.em)
      residuals <- theta_sample - fitted_values.em
      WSSR <- sum(w_ik*(residuals^2))#
      #var <- posterior_var
      sigma_temp <- sqrt((WSSR)/(sum(w_ik)-p))

      #joint response ???hg
      bic <- -2*(-n/2*log(2*pi*sigma_temp^2)- 1/(2*n_sam*sigma_temp^2) * sum(w_ik*(theta_sample - Z.em%*%coef_hat_em_debias)^2)) + sum(coef_hat_em_debias!=0)*log(n)

      delta_s <- abs(sigma_old-sigma_temp)
      delta_b <- max(abs(beta_gamma_old - beta_gamma_means))
      #delta_t <- abs(theta_t-theta_debiased)
      delta <- max(delta_s, delta_b)
      sigma_old <- sigma_temp
      cov_t <- cov
      beta_gamma_old <- beta_gamma_means
      #theta_t <- theta_debiased
      iter <- iter+1
      #print(c(delta))
    }
    #print(ll)
    #print(iter)
    results[[ll]] <- list(beta_gamma_old,sigma_old,bic,cov_t)
  }
  minBIC <- which.min(unlist(lapply(results, function(x) x[[3]])))
  return(list(
    coefficients  = results[[minBIC]][[1]],
    sigma = results[[minBIC]][[2]],
    LogLik = results[[minBIC]][[3]],
    minBIC = minBIC,
    Convergence = ifelse(iter<iter.max,"Converged","Did not converge")
  ))
}

