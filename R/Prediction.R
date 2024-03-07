#' Prediction of new observations
#'
#' @description
#' Computing a prediction for new observations
#'
#' @details
#' It provides a summary of the output of the ZIHR function, including parameter estimations.
#'
#'
#' @param object 	an object inheriting from class ZIHR
#' @param data dataset of observed variables with the same format as the data in the object
#'
#'
#' @return
#' Estimation, standard errors and 95% credible intervals for predictions
#'
#'
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZIHR.R
#'
#' @md
#' @seealso \code{\link{ZIHR}}
#' @export

Prediction <- function(object, data) {
  family <- object$family
  modelY <- object$modelY
  modelZ <- object$modelZ
  n <- length(data[, 1])
  z <- rep(0, n)
  y <- data[all.vars(modelY)][, 1]
  z[y == 0] <- 1

  mfX1 <- stats::model.frame(modelY, data = data)
  X1 <- stats::model.matrix(modelY, mfX1)

  ## Number of patients and number of longitudinal observations per patient
  Nbeta <- dim(X1)[2]


  data1 <- data.frame(z, data)

  mfX2 <- stats::model.frame(modelZ, data = data1)
  X2 <- stats::model.matrix(modelZ, mfX2)
  Nalpha <- dim(X2)[2]

  beta <- object$MCMC$beta
  alpha <- object$MCMC$alpha
  m <- dim(alpha)[1]


  if (family == "Gaussian") {
    Out <- "Ongoing"
  }
  ## $$
  if (family == "Gamma") {
    Out <- "Ongoing"
  }
  ## $$
  if (family == "Exponential") {
    Out <- "Ongoing"
  }

  if (family == "Weibull") {
    Out <- "Ongoing"
  }

  if (family == "Beta") {
    Out <- "Ongoing"
  }
  ## $$
  if (family == "inverse.gaussian") {
    Out <- "Ongoing"
  }
  ## $$
  if (family == "Poisson") {
    mu <- a <- b <- E <- matrix(0, n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        mu[i, j] <- exp(X1[i, ] %*% beta[j, ])
        a[i, j] <- (1 + exp(X2[i, ] %*% alpha[j, ]))^-1
        b[i, j] <- 1 - exp(-mu[i, j])
        E[i, j] <- a[i, j] * mu[i, j] / b[i, j]
      }
    }
    E_est <- apply(E, 1, mean)
    E_sd <- apply(E, 1, sd)
    E_CI <- apply(E, 1, quantile, c(0.025, 0.975))
    Out <- cbind(E_est, E_sd, t(E_CI), data)
    colnames(Out) <- c("est_pred", "sd_pred", "l_pred", "u_pred", names(data))
  }

  ## $$
  if (family == "NB") {
    r <- object$MCMC$Dispersion

    mu <- a <- b <- E <- matrix(0, n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        mu[i, j] <- exp(X1[i, ] %*% beta[j, ])
        a[i, j] <- (1 + exp(X2[i, ] %*% alpha[j, ]))^-1
        b[i, j] <- 1 - (r[j] / (r[j] + mu[i, j])^r[j])
        E[i, j] <- a[i, j] * mu[i, j] / b[i, j]
      }
    }
    E_est <- apply(E, 1, mean)
    E_sd <- apply(E, 1, sd)
    E_CI <- apply(E, 1, quantile, c(0.025, 0.975))
    Out <- cbind(E_est, E_sd, t(E_CI), data)
    colnames(Out) <- c("est_pred", "sd_pred", "l_pred", "u_pred", names(data))
  }

  ## $$
  if (family == "Logarithmic") {
    mu <- a <- E <- matrix(0, n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        mu[i, j] <- -(log(1 - exp(X1[i, ] %*% beta[j, ])))^-1 *
          exp(X1[i, ] %*% beta[j, ]) / (1 - exp(X1[i, ] %*% beta[j, ]))

        a[i, j] <- (1 + exp(X2[i, ] %*% alpha[j, ]))^-1
        E[i, j] <- a[i, j] * mu[i, j]
      }
    }
    E_est <- apply(E, 1, mean)
    E_sd <- apply(E, 1, sd)
    E_CI <- apply(E, 1, quantile, c(0.025, 0.975))
    Out <- cbind(E_est, E_sd, t(E_CI), data)
    colnames(Out) <- c("est_pred", "sd_pred", "l_pred", "u_pred", names(data))
  }

  ## $$
  if (family == "Binomial") {
    Out <- "Ongoing"
  }


  ## $$
  if (family == "Bell") {
    mu <- a <- b <- E <- matrix(0, n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        mu[i, j] <- exp(X1[i, ] %*% beta[j, ]) * exp(exp(X1[i, ] %*% beta[j, ]))
        a[i, j] <- (1 + exp(X2[i, ] %*% alpha[j, ]))^-1
        b[i, j] <- 1 - exp(-exp(mu[i, j]) + 1)
        E[i, j] <- a[i, j] * mu[i, j] / b[i, j]
      }
    }
    E_est <- apply(E, 1, mean)
    E_sd <- apply(E, 1, sd)
    E_CI <- apply(E, 1, quantile, c(0.025, 0.975))
    Out <- cbind(E_est, E_sd, t(E_CI), data)
    colnames(Out) <- c("est_pred", "sd_pred", "l_pred", "u_pred", names(data))
  }


  ## $$
  if (family == "GP") {
    phi <- object$MCMC$Dispersion

    mu <- a <- b <- E <- matrix(0, n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        mu[i, j] <- exp(X1[i, ] %*% beta[j, ])
        a[i, j] <- (1 + exp(X2[i, ] %*% alpha[j, ]))^-1
        b[i, j] <- 1 - exp(-mu[i, j] / (1 + phi[j] * mu[i, j]))
        E[i, j] <- a[i, j] * mu[i, j] / b[i, j]
      }
    }
    E_est <- apply(E, 1, mean)
    E_sd <- apply(E, 1, sd)
    E_CI <- apply(E, 1, quantile, c(0.025, 0.975))
    Out <- cbind(E_est, E_sd, t(E_CI), data)
    colnames(Out) <- c("est_pred", "sd_pred", "l_pred", "u_pred", names(data))
  }



  list(Prediction = Out)
}
