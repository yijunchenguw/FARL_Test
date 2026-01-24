library(paran)
library(lavaan)
library(mirt)
library(glmnet)
library(torch)
#' Fit a Latent Regression Model via the EM Algorithm
#'
#' Fits a latent regression model using an expectation--maximization (EM)
#' algorithm, where latent abilities are linked to observed covariates and
#' item response data through an IRT measurement model.
#'
#' @param X Numeric matrix. Design matrix of covariates (\eqn{N \times p}).
#' @param Y Matrix. Item response matrix (\eqn{N \times J}).
#' @param a Numeric vector or matrix. Item discrimination parameters.
#' @param d Numeric vector or matrix. Item difficulty (intercept) parameters.
#' @param p Integer. Number of covariates included in the latent regression.
#' @param n_sam Integer. Number of Monte Carlo or quadrature samples used to
#'   approximate posterior expectations in the E-step. Default is \code{30}.
#' @param verbose Logical. If \code{TRUE}, prints progress information during
#'   estimation. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coef}: Estimated coefficients from the latent regression model.
#'   \item \code{Sigma}: Estimated covariance matrix of the latent variables.
#'   \item \code{theta}: Posterior mean estimates of latent abilities.
#'   \item \code{converged}: Logical indicator of whether the EM algorithm
#'     converged.
#'   \item \code{iter}: Number of EM iterations performed.
#' }
#'
#' @details
#' The EM algorithm alternates between computing the conditional expectations
#' of latent abilities given the observed responses and current parameter
#' estimates (E-step) and maximizing the expected complete-data log-likelihood
#' with respect to the model parameters (M-step). Numerical optimization is
#' employed for item parameters, while closed-form updates are available for
#' latent regression and covariance parameters.
#'
#' @export
mml <- function(X, Y, parTab, n_sam = 5, method = "FARLR_EMM", lambda = seq(0.1, 0.5, by = 0.1),delta.criteria = 1e-3,iter.max = 200, window.size = 50, verbose = TRUE) {
  PA <- paran(X, iterations = 500, centile = 0, quiet = TRUE)
  K_hat <- PA$Retained
  fa <- factor.analysis(X, K_hat, method = "ml")
  Wupdate.t <- fa$Gamma
  Sgm_inv <- solve(diag(fa$Sigma))
  p = ncol(X)
  n = nrow(X)
  orthg <- t(Wupdate.t) %*% Sgm_inv %*% Wupdate.t/p
  V <- eigen(orthg)$vectors
  Wupdate <- t(Wupdate.t %*% V)
  Uupdate <- t(solve(Wupdate %*% Sgm_inv %*% t(Wupdate)) %*% Wupdate %*% Sgm_inv %*% t(X))
  Z <- cbind(Uupdate,X)
  colnames(Y) <- paste0("item", c(1:ncol(Y)), sep = "")
  storage.mode(Y) <- "integer"
  est_mirt <- suppressMessages(
    suppressWarnings(
      mirt(Y, 1, verbose = FALSE)
    )
  )
  theta_est_irt <- fscores(est_mirt,full.scores.SE = TRUE)
  theta_est_irt.mean <- theta_est_irt[,1]
  theta_est_irt.se <- theta_est_irt[,2]
  bin <- c(1, 2, 29, 15, 45)
  itemNames <- if (is.null(colnames(Y))) {
    paste0("i", sprintf("%03d", seq_len(ncol(Y))))
  } else {
    colnames(Y)
  }
  colnames(Y) <- itemNames
  subject <- factor(c(1:nrow(X)))
  stuItems <- reshape(data=data.frame(cbind(Y,subject)), varying=itemNames, idvar="subject",
                      direction="long", v.names="score",
                      times=itemNames, timevar="key")
  Y_back <- reshape(
    stuItems,
    idvar = "subject",
    timevar = "key",
    direction = "wide"
  )
  rownames(Y_back) <- subject
  Y_back <- Y_back[,-1]
  if (method == "FARLR_Debias") {

    hatB <- t( (1/n) * t(Uupdate) %*% X )
    hatU <- X - Uupdate %*% t(hatB)
    Fan  <- cbind(Uupdate, hatU)
    Z.em <- NA
    resp_rep <- NA
    fn <- farlr_debias

  } else if (method == "FARLR_EMM") {
    hatB <- NA
    hatU <- NA
    Fan  <- NA
    resp_rep <- rep(1, n_sam) %x% Y
    Z.em <- rep(1, n_sam) %x% Z
    fn <- farlr_emm

  } else {
    stop("Unknown method: ", method)
  }
  resultII <-fn(n, resp, parTab, K_hat, ncol(X), lambda_all = lambda, delta.criteria = 1e-3,iter.max = 200, n_sam = n_sam, window.size = 50,theta_est_irt.mean, theta_est_irt.se, resp_rep, Z.em, Uupdate,hatU, Fan,bin, verbose = TRUE)
  #resultII <- fn(nrow(Y), Y_back, parTab, K_hat, ncol(X), lambda, delta.criteria = 1e-3,iter.max = 200, n_sam = 30, window.size = 50,theta_est_irt.mean, theta_est_irt.se, resp_rep, Z.em, bin, verbose = TRUE)
  resultII$stuDat <- cbind(subject,Z)
  resultII$stuItems <- stuItems
  invisible(resultII)

  # hatB<-t(1/n*t(Uupdate)%*%X) #Estimated Factor Loading
  # hatU<-X-Uupdate%*%t(hatB)
  # Fan <- cbind(Uupdate,hatU)
  # resultDebias <- farlr_debias(nrow(Y_back), resp, parTab, K_hat, ncol(X), lambda, delta.criteria = 1e-3,iter.max = 200, n_sam = 10, window.size = 50,theta_est_irt.mean, theta_est_irt.se, resp_rep, Z.em, Uupdate,hatU, Fan,bin, verbose = TRUE)
  # resultII$stuDat <- cbind(subject,Z)
  # resultII$stuItems <- stuItems
  # invisible(resultII)

}
mml_test <- function(){
  mmlcomp <-  with(sim_a1, mml(X, Y, parTab, method = "FARLR_Debias")) #mml(X = sim_a1$X, Y = sim_a1$Y, parTab = sim_a1$parTab,  method = "FARLR_Debias")
  colnames(sim_a1$X) <- paste0("X", c(1:ncol(sim_a1$X)))
  mmlcomp$X <- sim_a1$X
  mmlcomp$item_params <- sim_a1$parTab
  invisible(mmlcomp)
  PVs <- drawPVs(mmlcomp, 10L)
  return(PVs)
}
mml_test2 <- function(){
  with(sim_a1, mml(X, Y, parTab, method = "FARLR_Debias"))
  }
