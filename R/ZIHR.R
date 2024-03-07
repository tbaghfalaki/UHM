#' Zero-inflation hurdle regression models
#'
#' @description
#' Fits zero-inflated hurdle regression models
#'
#' @details
#' A function utilizing the 'JAGS' software to estimate the linear hurdle regression model.
#'
#' @param modelY a formula for the mean of the count response. This argument is identical to the one in the "glm" function.
#' @param modelZ a formula for the probability of zero. This argument is identical to the one in the "glm" function.
#' @param data data set of observed variables.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param family Family objects streamline the specification of model details for functions like glm. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson models, respectively, while the others are utilized for the corresponding models based on their names.
#'
#' @importFrom
#'
#' @return
#'
#' - MCMC chains for the unknown parameters
#' - Est list of posterior mean for each parameter
#' - SD list of standard error for each parameter
#' - L_CI list of 2.5th percentiles of the posterior distribution serves as the lower bound of the Bayesian credible interval
#' - U_CI list of 97.5th percentiles of the posterior distribution serves as the lower bound of the Bayesian credible interval
#' - Rhat Gelman and Rubin diagnostic for all parameter
#' - beta the regression coefficients of mean of the hurdle model
#' - alpha the regression coefficients of probability of the hurdle model
#' - The variance, over-dispersion, dispersion, or scale parameters of models depend on the family used
#' - DIC deviance information criterion
#' - LPML Log Pseudo Marginal Likelihood (LPML) criterion
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZIHR.R
#'
#' @md
#' @export

ZIHR <- function(modelY, modelZ, data, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, family = "Gaussian") {
  n <- length(data[, 1])
  z <- rep(0, n)
  y <- data[all.vars(modelY)][, 1]
  z[y == 0] <- 1
  mfX1 <- stats::model.frame(modelY, data = data)
  X1 <- stats::model.matrix(modelY, mfX1)

  ## Number of patients and number of longitudinal observations per patient
  Nbeta <- dim(X1)[2]

  data <- data.frame(z, data)

  mfX2 <- stats::model.frame(modelZ, data = data)
  X2 <- stats::model.matrix(modelZ, mfX2)
  Nalpha <- dim(X2)[2]

  K <- sum(y) + 10000
  ###########
  GP <- "model{

  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

 ll[i]<-(1-z[i])*(y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(lambda[i])<-inprod(beta[],X[i,])
     logit(muz[i])<-inprod(alpha[],X2[i,])
  }

  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  phiz~dnorm(0,.001)


  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  #####


  Bell <- "model{
for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

    log(theta[i])<- inprod(beta[],X[i,])
    logit(muz[i])<-inprod(alpha[],X2[i,])
}

  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  ###########
  binomial <- "model{

  m<-max(y)
  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(p[i])+(m-y[i])*log(1-p[i])-
      log(1-pow(1-p[i],m)))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])
           logit(p[i])<-inprod(beta[],X[i,])
    logit(muz[i])<-inprod(alpha[],X2[i,])

  }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  ###########
  Gamma <- "model{

  for(i in 1:n){
 cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

     zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


    mu1[i]<-sigma1/mu[i]
    log(mu[i])<-inprod(beta[],X[i,])
        logit(muz[i])<-inprod(alpha[],X2[i,])

  }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  sigma1~dgamma(.1,.1)
  sigma<-1/sigma1
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"


  ###########
  Exp <- "model{

  for(i in 1:n){
 cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

     zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


    log(mu[i])<-inprod(beta[],X[i,])
        logit(muz[i])<-inprod(alpha[],X2[i,])

  }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"


  ###########
  Beta <- "model{

  for(i in 1:n){
 cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

     zeros[i]~dpois(phiz[i])
    phiz[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.beta(y[i],a[i], b[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    a[i]<-mu[i]*phi1
    b[i]<-(1-mu[i])*phi1
    logit(mu[i])<-inprod(beta[],X[i,])
        logit(muz[i])<-inprod(alpha[],X2[i,])

  }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  phi1~dgamma(.1,.1)
  phi<-1/phi1
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"


  ###########
  Weibull <- "model{

  for(i in 1:n){
 cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

     zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


    log(mu[i])<-inprod(beta[],X[i,])
        logit(muz[i])<-inprod(alpha[],X2[i,])

  }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  kappa~dgamma(.1,.1)
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"


  ###########
  Gaussian <- "model{

  for(i in 1:n){
   cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K
 ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


    mu[i]<-inprod(beta[],X[i,])
        logit(muz[i])<-inprod(alpha[],X2[i,])

 }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  tau~dgamma(.01,.01)
  sigma<-1/tau
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  ###########
  IGauss <- "model{

  for(i in 1:n){
        cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])



    log(mu[i])<-inprod(beta[],X[i,])
    logit(muz[i])<-inprod(alpha[],X2[i,])
 }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  lambda~dgamma(.1,.1)
  sigma<-1/lambda
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  ###########
  logar <- "model{


  for(i in 1:n){
        cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-(1-z[i])*(log(-1/log(1-p[i]))+y[i]*log(p[i])-log(y[i]))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


 logit(p[i])<-inprod(beta[],X[i,])
     logit(muz[i])<-inprod(alpha[],X2[i,])

 }

  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  ###########
  NB <- "model{


  for(i in 1:n){
        cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-(1-z[i])*(loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+y[i]*log(lambda[i]/(lambda[i]+r))-
    log(1-pow(r/(r+lambda[i]),r)))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(lambda[i])<-inprod(beta[],X[i,])
         logit(muz[i])<-inprod(alpha[],X2[i,])

 }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  r~ dgamma(.1,.1)

  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"

  ###########
  Poisson <- "model{


  for(i in 1:n){
        cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-(1-z[i])*(y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(lambda[i])<-inprod(beta[],X[i,])
         logit(muz[i])<-inprod(alpha[],X2[i,])

 }


  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
}"
  ##############


  if (family == "Gaussian") {
    model.file <- textConnection(Gaussian)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta), tau = 1)
    }
    parameters <- c("alpha", "beta", "sigma", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, zeros = rep(0, n), K = K, Nalpha = Nalpha)

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, sigma = sim1$sims.list$sigma)
    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Variance"

      Variance <- cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)

      colnames(Variance) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Variance = Variance)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Variance"

      Variance <- cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)

      colnames(Variance) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Variance = Variance)
    }
  }

  ## $$
  if (family == "Weibull") {
    model.file <- textConnection(Weibull)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta), kappa = 1)
    }
    parameters <- c("alpha", "beta", "kappa", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, zeros = rep(0, n), K = K, Nalpha = Nalpha)

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, sigma = sim1$sims.list$sigma)
    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <-
        names(sim1$Rhat$kappa) <- "Scale"

      Scale <- cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa, sim1$Rhat$kappa)

      colnames(Scale) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Scale = Scale)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <- "Scale"

      Scale <- cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa)

      colnames(Scale) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Scale = Scale)
    }
  }

  ## $$
  if (family == "Gamma") {
    model.file <- textConnection(Gamma)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta), sigma = 1)
    }
    parameters <- c("alpha", "beta", "sigma", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, zeros = rep(0, n), K = K, Nalpha = Nalpha)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, sigma = sim1$sims.list$sigma)


    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Scale"

      Scale <- cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)

      colnames(Scale) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Scale = Scale)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Scale"

      Scale <- cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)

      colnames(Scale) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Scale = Scale)
    }
  }

  ## $$
  if (family == "Beta") {
    model.file <- textConnection(Beta)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta), phi1 = 1)
    }
    parameters <- c("alpha", "beta", "phi", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, zeros = rep(0, n), K = K, Nalpha = Nalpha)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, phi = sim1$sims.list$phi)


    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$phi) <-
        names(sim1$sd$phi) <-
        names(sim1$q2.5$phi) <-
        names(sim1$q97.5$phi) <-
        names(sim1$Rhat$phi) <- "Dispersion"

      Dispersion <- cbind(sim1$mean$phi, sim1$sd$phi, sim1$q2.5$phi, sim1$q97.5$phi, sim1$Rhat$phi)

      colnames(Dispersion) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Dispersion = Dispersion)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$phi) <-
        names(sim1$sd$phi) <-
        names(sim1$q2.5$phi) <-
        names(sim1$q97.5$phi) <- "Dispersion"

      Dispersion <- cbind(sim1$mean$phi, sim1$sd$phi, sim1$q2.5$phi, sim1$q97.5$phi)

      colnames(Dispersion) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Dispersion = Dispersion)
    }
  }

  ## $$
  if (family == "Exponential") {
    model.file <- textConnection(Exp)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta))
    }
    parameters <- c("alpha", "beta", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, zeros = rep(0, n), K = K, Nalpha = Nalpha)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha)


    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)


      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    }
  }

  ## $$
  if (family == "inverse.gaussian") {
    model.file <- textConnection(IGauss)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta), lambda = 1)
    }
    parameters <- c("alpha", "beta", "sigma", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha, zeros = rep(0, n), K = K)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, sigma = sim1$sims.list$sigma)

    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Shape parameter"

      Shape_parameter <- cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)

      colnames(Shape_parameter) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Shape_parameter = Shape_parameter)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Shape parameter"

      Shape_parameter <- cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)

      colnames(Shape_parameter) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Shape_parameter = Shape_parameter)
    }
  }


  ## $$
  if (family == "Poisson") {
    model.file <- textConnection(Poisson)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta))
    }
    parameters <- c("alpha", "beta", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha, zeros = rep(0, n), K = K)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )
    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha)

    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    }
  }

  ## $$
  if (family == "NB") {
    model.file <- textConnection(NB)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta), r = 1)
    }
    parameters <- c("alpha", "beta", "r", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha, zeros = rep(0, n), K = K)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )
    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, Dispersion = sim1$sims.list$r)

    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"

      Dispersion <- cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r, sim1$Rhat$r)

      colnames(Dispersion) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Dispersion = Dispersion)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <- "Dispersion"

      Dispersion <- cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r)

      colnames(Dispersion) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Dispersion = Dispersion)
    }
  }

  ## $$
  if (family == "Logarithmic") {
    model.file <- textConnection(logar)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta))
    }
    parameters <- c("alpha", "beta", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha, zeros = rep(0, n), K = K)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )
    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha)


    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    }
  }

  ## $$
  if (family == "Binomial") {
    model.file <- textConnection(binomial)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta))
    }

    parameters <- c("alpha", "beta", "cpoinvy")

    d.jags <- list(n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha, zeros = rep(0, n), K = K)
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )
    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha)

    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    }
  }


  ## $$
  if (family == "Bell") {
    model.file <- textConnection(Bell)

    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta))
    }

    parameters <- c("alpha", "beta", "cpoinvy")



    C <- c()
    for (i in 1:n) {
      C[i] <- log(numbers::bell(y[i])) - lfactorial(y[i])
    }


    d.jags <- list(
      n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha,
      zeros = rep(0, n), K = K, C = C
    )

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )
    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha)

    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)


      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)


      results <- list(Count_model = MM, Zero_inflated_model = ZPM)
    }
  }
  ####



  ## $$
  if (family == "GP") {
    model.file <- textConnection(GP)
    i.jags <- function() {
      list(alpha = stats::rnorm(Nalpha), beta = stats::rnorm(Nbeta))
    }

    parameters <- c("alpha", "beta", "phiz", "cpoinvy")

    d.jags <- list(
      n = n, y = y, X = X1, Nbeta = Nbeta, z = z, X2 = X2, Nalpha = Nalpha,
      zeros = rep(0, n), K = K
    )
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = NULL,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )
    MCMC <- list(beta = sim1$sims.list$beta, alpha = sim1$sims.list$alpha, Dispersion = sim1$sims.list$phiz)

    if (n.chains > 1) {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <-
        names(sim1$Rhat$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <-
        names(sim1$Rhat$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha, sim1$Rhat$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta, sim1$Rhat$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"

      Dispersion <- cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz, sim1$Rhat$phiz)
      colnames(Dispersion) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Dispersion = Dispersion)
    } else {
      names(sim1$mean$beta) <-
        names(sim1$sd$beta) <-
        names(sim1$q2.5$beta) <-
        names(sim1$q97.5$beta) <- colnames(X1)

      names(sim1$mean$alpha) <-
        names(sim1$sd$alpha) <-
        names(sim1$q2.5$alpha) <-
        names(sim1$q97.5$alpha) <- colnames(X2)

      ZPM <- cbind(sim1$mean$alpha, sim1$sd$alpha, sim1$q2.5$alpha, sim1$q97.5$alpha)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$beta, sim1$sd$beta, sim1$q2.5$beta, sim1$q97.5$beta)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM)

      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <- "Dispersion"

      Dispersion <- cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz)

      colnames(Dispersion) <- c("Est", "SD", "L_CI", "U_CI")
      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Dispersion = Dispersion)
    }
  }


  ##################

  DIC <- sim1$DIC - 2 * n * K

  LPML <- -sum(log(sim1$mean$cpoinvy))

  list(
    modelY = modelY, modelZ = modelZ,
    family = family, MCMC = MCMC, Estimation = results,
    DIC = DIC, LPML = LPML
  )
}
