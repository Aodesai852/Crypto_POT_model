options(warn=-1) #suppress the warnings for numerical optimization
####################################################################################################
# CDF, PDF and quantiles of truncated observations GPD (with the flexibility of specifying an extra k)
# PDF of GPD
dGPD <- function(y, tau, k, sigma, xi){
  indicator <- y>0
  result <- (indicator>0)*(k*xi*sigma^xi*(tau+y)^(-xi-1)) + (indicator==0)*(1-k*(tau/sigma)^(-xi))
  return(result)
}
# quantile of GPD
qGPD <- function(q, tau, k, sigma, xi){
  indicator <- q<1-k*(tau/sigma)^(-xi)
  (indicator>0)*tau + (indicator==0)*((k/(1-q))^(1/xi)*sigma)
} 
# CDF of GPD
pGPD <- function(y, tau, k, sigma, xi){
  1-k*(tau/sigma)^(-xi)*(1+y/tau)^(-xi)
}
#######################################
# Truncated t or normal distribution
trunc_normal <- function(tau, sigma){
  u <- runif(1)
  while(u>pnorm(tau, sd=sigma, mean=0)){
    u <- runif(1)
  }
  return(qnorm(u, sd=sigma, mean=0))
}

trunc_t <- function(tau, sigma, xi){
  u <- runif(1)
  while(u>pt(tau/sigma, df=xi)){
    u <- runif(1)
  }
  return(sigma*qt(u, df=xi))
}
#######################################
### Generate the underlying dynamic univariate POT model###
#######################################
# Simulate the dynamic univariate POT-GPD time series with t distribution for e_t
sim_univPOT_t <-function(true_para, total, burnin=1000,Sigma_type){
  dyn_para_xi <- true_para[[1]]; dyn_para_sigma <- true_para[[2]]
  phi1 <- dyn_para_xi[1]; phi2 <- dyn_para_xi[2]; phi3 <- dyn_para_xi[3];
  psi1 <- dyn_para_sigma[1]; psi2 <- dyn_para_sigma[2]; psi3 <- dyn_para_sigma[3]
  
  N <- total + burnin
  raw_obs <- xi_obs <- sigma_obs <- rep(0, N)
  xi_obs[1] <- 4; sigma_obs[1] <- 0.5; raw_obs[1] <- sigma_obs[1]*rt(n=1, df=xi_obs[1])
  # cur_para <- c(tau, k, sigma_obs[1], xi_obs[1])
  if (Sigma_type == "Log_Sigma") {
    for(t in 2:N){
      update1 <- psi3*exp(-abs(raw_obs[t-1]))
      update2 <- phi3*exp(-abs(raw_obs[t-1]))
      # print(update)
      sigma_obs[t] <- exp(psi1 + psi2*log(sigma_obs[t-1]) + update1)
      xi_obs[t] <- exp(phi1 + phi2*log(xi_obs[t-1]) + update2)
      # cur_para <- c(tau, k, sigma_obs[t], xi_obs[t])
      raw_obs[t] <- sigma_obs[t]*rt(n=1, df=xi_obs[t])
      raw_obs[t] <- min(raw_obs[t],10); raw_obs[t] <- max(raw_obs[t],-10)
    }
  } else if (Sigma_type == "Square_Sigma") {
    for(t in 2:N){
      update1 <- psi3*(raw_obs[t-1])^2
      update2 <- phi3*exp(-abs(raw_obs[t-1]))
      # print(update)
      sigma_obs[t] <- sqrt(psi1 + psi2*(sigma_obs[t-1])^2 + update1)
      xi_obs[t] <- exp(phi1 + phi2*log(xi_obs[t-1]) + update2)
      # cur_para <- c(tau, k, sigma_obs[t], xi_obs[t])
      raw_obs[t] <- sigma_obs[t]*rt(n=1, df=xi_obs[t])
      raw_obs[t] <- min(raw_obs[t],10); raw_obs[t] <- max(raw_obs[t],-10)
    }  
  } 
  else {
    stop("Invalid `Sigma_type` value")
  }
  return(list(raw_data=raw_obs[-c(1:burnin)], xi=xi_obs[-c(1:burnin)], sigma=sigma_obs[-c(1:burnin)]))
}

# Simulate the dynamic univariate POT-GPD time series with emprical distribution for e_t
sim_univPOT_garch <-function(para, total, empirical, burnin=1000, Sigma_type){
  GPD_para <- para[[1]]; dyn_para_xi <- para[[2]]; dyn_para_sigma <- para[[3]]
  tau <- GPD_para[1]; k <- GPD_para[2]
  phi1 <- dyn_para_xi[1]; phi2 <- dyn_para_xi[2]; phi3 <- dyn_para_xi[3]
  psi1 <- dyn_para_sigma[1]; psi2 <- dyn_para_sigma[2]; psi3 <- dyn_para_sigma[3]
  
  N <- total + burnin
  raw_obs <- xi_obs <- sigma_obs <- rep(0, N)
  xi_obs[1] <- 4;  sigma_obs[1] <- 0.5;  raw_obs[1] <- sample(empirical, 1)
  cur_para <- c(tau, k, sigma_obs[1], xi_obs[1])
  if (Sigma_type == "Log_Sigma") {
    for(t in 2:N){
      update1 <- psi3*exp(-abs(raw_obs[t-1]))
      update2 <- phi3*exp(-abs(raw_obs[t-1]))
      # print(update)
      sigma_obs[t] <- exp(psi1 + psi2*log(sigma_obs[t-1]) + update1)
      xi_obs[t] <- exp(phi1 + phi2*log(xi_obs[t-1]) + update2)
      cur_para <- c(tau, k, sigma_obs[t], xi_obs[t])
      randomU <- runif(1)
      exceed_prob <- 1-k*(tau/sigma_obs[t])^(-xi_obs[t])
      if(randomU < exceed_prob){
        raw_obs[t] <- quantile(empirical, randomU/exceed_prob)
      } else {
        raw_obs[t] <- min(qGPD(randomU, tau, k, sigma_obs[t], xi_obs[t]), 8.2) # trim the GPD
      }
    }
  } else if (Sigma_type == "Square_Sigma") {
    for(t in 2:N){
      update1 <- psi3*(raw_obs[t-1])^2
      update2 <- phi3*exp(-abs(raw_obs[t-1]))
      # print(update)
      sigma_obs[t] <- sqrt(psi1 + psi2*(sigma_obs[t-1])^2 + update1)
      xi_obs[t] <- exp(phi1 + phi2*log(xi_obs[t-1]) + update2)
      cur_para <- c(tau, k, sigma_obs[t], xi_obs[t])
      randomU <- runif(1)
      exceed_prob <- 1-k*(tau/sigma_obs[t])^(-xi_obs[t])
      if(randomU < exceed_prob){
        raw_obs[t] <- quantile(empirical, randomU/exceed_prob)
      } else {
        raw_obs[t] <- min(qGPD(randomU, tau, k, sigma_obs[t], xi_obs[t]), 8.2) # trim the GPD
      }
    }
  }
  else {
    stop("Invalid `Sigma_type` value")
  }  
  return(list(raw_data=raw_obs[-c(1:burnin)], xi=xi_obs[-c(1:burnin)], sigma=sigma_obs[-c(1:burnin)]))
}

####################################################################################################
# Recover dynamics of sigma, tail index with raw observations from the dynamic univariate POT-GPD time series
recoverXiSigma_univPOT_garch<- function(para, raw_obs, init=NULL, Sigma_type){
  GPD_para <- para[[1]]; dyn_para_xi <- para[[2]]; dyn_para_sigma <- para[[3]]
  tau <- GPD_para[1]; k <- GPD_para[2]
  phi1 <- dyn_para_xi[1]; phi2 <- dyn_para_xi[2]; phi3 <- dyn_para_xi[3]
  psi1 <- dyn_para_sigma[1]; psi2 <- dyn_para_sigma[2]; psi3 <- dyn_para_sigma[3]
  N <- length(raw_obs); xi_obs <- sigma_obs <- rep(0, N)
  if(is.null(init)){ # initial values for xi and sigma
    xi_obs[1] <- 4
    sigma_obs[1] <- 0.5
  } else {
    xi_obs[1] <- init[1]; sigma_obs[1] <- init[2]
  }
  # constrict xi and sigma in [lower_bound, upper_bound] for numerical stability
  upper_bound_sigma <- tau; lower_bound_sigma <- 1e-2 # model resrriction: the threshold tau needs to be larger than the scale parameter sigma
  upper_bound_xi <- 10; lower_bound_xi <- 1e-2
  if (Sigma_type == "Log_Sigma") {
    for(t in 2:N){
      update1 <- psi3*exp(-abs(raw_obs[t-1]))
      update2 <- phi3*exp(-abs(raw_obs[t-1]))
      # print(update)
      sigma_obs[t] <- exp(psi1 + psi2*log(sigma_obs[t-1]) + update1)
      xi_obs[t] <- exp(phi1 + phi2*log(xi_obs[t-1]) + update2)
      if(xi_obs[t]>upper_bound_xi|is.na(xi_obs[t])){xi_obs[t] <- upper_bound_xi}
      if(xi_obs[t]<lower_bound_xi|is.na(xi_obs[t])){xi_obs[t] <- lower_bound_xi}
      if(sigma_obs[t]>upper_bound_sigma|is.na(sigma_obs[t])){sigma_obs[t] <- upper_bound_sigma}
      if(sigma_obs[t]<lower_bound_sigma|is.na(sigma_obs[t])){sigma_obs[t] <- lower_bound_sigma}
    }
  }else if (Sigma_type == "Square_Sigma") {
    for(t in 2:N){
      update1 <- psi3*(raw_obs[t-1])^2
      update2 <- phi3*exp(-abs(raw_obs[t-1]))
      # print(update)
      sigma_obs[t] <- sqrt(psi1 + psi2*(sigma_obs[t-1])^2 + update1)
      xi_obs[t] <- exp(phi1 + phi2*log(xi_obs[t-1]) + update2)
      if(xi_obs[t]>upper_bound_xi|is.na(xi_obs[t])){xi_obs[t] <- upper_bound_xi}
      if(xi_obs[t]<lower_bound_xi|is.na(xi_obs[t])){xi_obs[t] <- lower_bound_xi}
      if(sigma_obs[t]>upper_bound_sigma|is.na(sigma_obs[t])){sigma_obs[t] <- upper_bound_sigma}
      if(sigma_obs[t]<lower_bound_sigma|is.na(sigma_obs[t])){sigma_obs[t] <- lower_bound_sigma}
    }    
  } else {
    stop("Invalid `Sigma_type` value")
  }
  return(list(xi=xi_obs, sigma=sigma_obs))
}
####################################################################################################
# Recover the CDF of Y (GPD) based on the marginal parameters
recoverUniFy <- function(y, k, tau, xiSigma_marg){
  y <- y; 
  xi <- xiSigma_marg$xi; sigma <- xiSigma_marg$sigma
  Fy <- pGPD(y, tau, k, sigma, xi)
  return(Fy)
}
####################################################################################################
# Likelihood function for the dynamic univariate POT-GPD time series with emprical distribution for e_t 
loglikUnivariate_univPOT_garch <- function(para, tau, y_obs, raw_obs, init=NULL, burnin=10, ktype='nonfixed', Sigma_type){ # init: initial values for xi and sigma
  if(ktype=='nonfixed'){ # k！=1
    k <- para[1]
    GPD_para <- c(tau, k); dyn_para_xi <- para[2:4]; dyn_para_sigma <- para[-c(1:4)]
    para <- list(GPD_para, dyn_para_xi, dyn_para_sigma)
    tmp <- recoverXiSigma_univPOT_garch(para, raw_obs, init, Sigma_type)
    xi <- tmp$xi; sigma <- tmp$sigma
  }else{
    para <- c(1,para) # k=1
    k <- para[1]
    GPD_para <- c(tau, k); dyn_para_xi <- para[2:4]; dyn_para_sigma <- para[-c(1:4)]
    para <- list(GPD_para, dyn_para_xi, dyn_para_sigma)
    tmp <- recoverXiSigma_univPOT_garch(para, raw_obs, init, Sigma_type)
    xi <- tmp$xi; sigma <- tmp$sigma
  }
  loglik <- -sum(log(dGPD(y_obs[-c(1:burnin)], tau, k, sigma[-c(1:burnin)], xi[-c(1:burnin)])))
  return(loglik)
}

####################################################################################################
# initial estimation of sigma and alpha
initEstUniGPDMarg <- function(para, y, tau){
  xi <- para[1]
  sigma <- para[2]
  k <- 1
  result <- -sum(log(dGPD(y, tau, k, sigma, xi)))
  return(result)
}

####################################################################################################
# initial estimation of sigma and alpha
initEstUniGPDMarg_ktype <- function(para, y, tau, ktype){
  if(ktype=='nonfixed'){ # k！=1
    k <- para[1]
    xi <- para[2]
    sigma <- para[3]
    result <- -sum(log(dGPD(y, tau, k, sigma, xi)))
  }else{
    xi <- para[1]
    sigma <- para[2]
    k <- 1
    result <- -sum(log(dGPD(y, tau, k, sigma, xi)))
  }
  return(result)
}