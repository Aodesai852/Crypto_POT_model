#' @title Connectedness Approach
#' @description This function provides a modular framework combining various models and connectedness frameworks.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param nfore H-step ahead forecast horizon
#' @param window.size Rolling-window size or Bayes Prior sample size
#' @param model Estimation model
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param connectedness Type of connectedness approach
#' @param VAR_config Config for VAR model
#' @param DCC_config Config for DCC-GARCH model
#' @param Connectedness_config Config for connectedness approach
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("acg2020")
#' dca = ConnectednessApproach(acg2020, 
#'                             nlag=1, 
#'                             nfore=12,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.96,
#'                                             prior="MinnesotaPrior", gamma=0.1)))
#' dca$TABLE
#' }
#' @import progress
#' @importFrom rugarch ugarchspec
#' @importFrom rugarch multispec
#' @importFrom rmgarch dccspec
#' @importFrom rmgarch dccfit
#' @references
#' Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal, 119(534), 158-171.\\
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting, 28(1), 57-66.\\
#' Barunik, J., & Krehlik, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics, 16(2), 271-296.\\
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management, 13(4), 84.\\
#' Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.\\
#' Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy, 73, 102219.\\
#' Chatziantoniou, I., & Gabauer, D. (2021). EMU risk-synchronisation and financial fragility through the prism of dynamic connectedness. The Quarterly Review of Economics and Finance, 79, 1-14.\\
#' Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters, 204, 109891.\\
#' Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management, 60, 100680.\\
#' Gabauer, D., Gupta, R., Marfatia, H., & Miller, S. (2020). Estimating US Housing Price Network Connectedness: Evidence from Dynamic Elastic Net, Lasso, and Ridge Vector Autoregressive Models (No. 202065). University of Pretoria, Department of Economics.\\
#' Cunado, J, Chatziantoniou, I., Gabauer, D., Hardik, M., & de Garcia, F.P. (2022). Dynamic spillovers across precious metals and energy realized volatilities: Evidence from quantile extended joint connectedness measures.\\
#' Chatziantoniou, I., Gabauer, D., & Gupta, R. (2021). Integration and Risk Transmission in the Market for Crude Oil: A Time-Varying Parameter Frequency Connectedness Approach (No. 202147).\\
#' Chatziantoniou, I., Aikins Abakah, E. J., Gabauer, D., & Tiwari, A. K. (2021). Quantile time-frequency price connectedness between green bond, green equity, sustainable investments and clean energy markets: Implications for eco-friendly investors. Available at SSRN 3970746.\\
#' Gabauer, D. (2020). Volatility impulse response analysis for DCC‐GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.
#' @author David Gabauer
#' @export
ConnectednessApproach = function(x,
                                 nlag=1, 
                                 nfore=10, 
                                 window.size=NULL, 
                                 corrected=FALSE,
                                 model=c("VAR", "QVAR", "LASSO", "Ridge", "Elastic", "TVP-VAR", "DCC-GARCH"),
                                 connectedness=c("Time","Frequency", "Joint", "Extended Joint"),
                                 VAR_config=list(
                                   QVAR=list(tau=0.5),
                                   ElasticNet=list(nfolds=10, alpha=NULL, loss="mae", delta_alpha=0.1),
                                   TVPVAR=list(kappa1=0.99, kappa2=0.99, prior="BayesPrior", gamma=0.01)),
                                 DCC_config=list(standardize=FALSE),
                                 Connectedness_config = list(
                                   TimeConnectedness=list(generalized=TRUE),
                                   FrequencyConnectedness=list(partition=c(pi,pi/2,0), generalized=TRUE, scenario="ABS")
                                 )) {
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  model = match.arg(model)
  connectedness = match.arg(connectedness)
  
  NAMES = colnames(x)
  k = ncol(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  t = nrow(x)
  if (is.null(window.size)) {
    window.size = nrow(x)
    t0 = 1
  } else {
    window.size = window.size - nlag
    t0 = t - window.size + 1
  }

  if (model=="VAR") {
    var_model = VAR
    configuration = list(nlag=nlag)
  } else if (model=="QVAR") {
    var_model = QVAR
    configuration = list(nlag=nlag, tau=VAR_config$QVAR$tau)
  } else if (model=="LASSO") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=1, nfolds=VAR_config$ElasticNet$nfolds, loss=VAR_config$ElasticNet$loss)
  } else if (model=="Ridge") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=0, nfolds=VAR_config$ElasticNet$nfolds, loss=VAR_config$ElasticNet$loss)
  } else if (model=="Elastic") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=VAR_config$ElasticNet$alpha, nfolds=VAR_config$ElasticNet$nfolds,
                         loss=VAR_config$ElasticNet$loss, delta_alpha=VAR_config$ElasticNet$delta_alpha)
  } else if (model=="TVP-VAR") {
    prior_ = VAR_config$TVPVAR$prior
    if (prior_=="MinnesotaPrior") {
      prior = MinnesotaPrior(gamma=VAR_config$TVPVAR$gamma, k=k, nlag=nlag)
    } else if (prior_=="UninformativePrior") {
      prior = UninformativePrior(k=k, nlag=nlag)
    } else if (prior_=="BayesPrior") {
	   print("#####到prior上面######")
      prior = BayesPrior(x=x, size=window.size, nlag=nlag)
	  print("#####到prior######")
    } else {
      stop("Error Prior choice")
    }
    configuration = list(l=c(VAR_config$TVPVAR$kappa1, VAR_config$TVPVAR$kappa2), nlag=nlag, prior=prior)
    var_model = TVPVAR
  } else if (model=="DCC-GARCH") {
    ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=TRUE),
                             variance.model=list(garchOrder=c(1,1), model="sGARCH"),
                             distribution.model="norm")
    mgarch.spec = rmgarch::dccspec(uspec=rugarch::multispec(replicate(k, ugarch.spec)),
                          dccOrder=c(1,1), distribution="mvnorm")
  } else {
    stop("Model does not exist")
  }

  message("Estimating model")
  if (model=="TVP-VAR") {
    fit = var_model(x, configuration=configuration)
    B_t = fit$B_t
    Q_t = fit$Q_t
	print("##############到fit#####################")
  } else if (model=="DCC-GARCH") {
    fit = rmgarch::dccfit(mgarch.spec, data=x)
    fevd = VFEVD(fit, nfore=nfore, standardize=DCC_config$standardize)$FEVD
    Q_t = fevd
  } else {
    Q_t = array(NA, c(k, k, t0))
    B_t = array(NA, c(k, k*nlag, t0))
    pb = progress_bar$new(total=t0)
    for (i in 1:t0) {
	
      fit = var_model(x[i:(i+window.size-1),], configuration=configuration)
	   
	  #print(fit$B)
	  
      B_t[,,i] = fit$B
      Q_t[,,i] = fit$Q
      pb$tick()
    }
  }
  DATE = as.character(zoo::index(x))
  date = DATE[(length(DATE)-dim(Q_t)[3]+1):length(DATE)]
  dimnames(Q_t)[[1]] = dimnames(Q_t)[[2]] = NAMES
  dimnames(Q_t)[[3]] = as.character(date)
  
  message("Computing connectedness measures")
  if (connectedness=="Time") {
    generalized = Connectedness_config$TimeConnectedness$generalized
    if (model=="DCC-GARCH") {
      dca = TimeConnectedness(FEVD=fevd, corrected=corrected)
      message("The DCC-GARCH connectedness approach is implemented according to:\n Gabauer, D. (2020). Volatility impulse response analysis for DCC-GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.")
    } else {
      dca = TimeConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore,
                              generalized=generalized,
                              corrected=corrected)
      if (model=="VAR" && !generalized) {
        message("The (orthogonalized) VAR connectedness approach is implemented according to:\n Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal, 119(534), 158-171.")
      } else if (model=="VAR" && generalized) {
        message("The (generalized) VAR connectedness approach is implemented according to:\n Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting, 28(1), 57-66.")
      } else if (model=="TVP-VAR") {
        message("The TVP-VAR connectedness approach is implemented according to:\n Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management, 13(4), 84.")
      } else if (model=="QVAR") {
        message("The QVAR connectedness approach is implemented according to:\n Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters, 204, 109891.")
      } else if (model=="LASSO" || model=="Ridge" || model=="Elastic") {
        message("The Elastic Net and its restricted models, namely, the LASSO and Ridge VAR connectedness approach are implemented according to:\n Gabauer, D., Gupta, R., Marfatia, H., & Miller, S. (2020). Estimating US Housing Price Network Connectedness: Evidence from Dynamic Elastic Net, Lasso, and Ridge Vector Autoregressive Models (No. 202065). University of Pretoria, Department of Economics.")
      }
    }
  } else if (connectedness=="Frequency") {
    dca = FrequencyConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore,
                                 partition=Connectedness_config$FrequencyConnectedness$partition,
                                 generalized=Connectedness_config$FrequencyConnectedness$generalized,
                                 scenario=Connectedness_config$FrequencyConnectedness$scenario, 
                                 corrected=corrected)
    if (model=="VAR") {
      message("The VAR frequency connectedness approach is implemented according to:\n Barunik, J., & Krehlik, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics, 16(2), 271-296.")
    } else if (model=="TVP-VAR") {
      message("The TVP-VAR frequency connectedness approach is implemented according to:\n Chatziantoniou, I., Gabauer, D., & Gupta, R. (2021). Integration and Risk Transmission in the Market for Crude Oil: A Time-Varying Parameter Frequency Connectedness Approach (No. 202147).")
    } else if (model=="QVAR") {
      message("The QVAR frequency connectedness approach is implemented according to:\n Chatziantoniou, I., Aikins Abakah, E. J., Gabauer, D., & Tiwari, A. K. (2021). Quantile time-frequency price connectedness between green bond, green equity, sustainable investments and clean energy markets: Implications for eco-friendly investors. Available at SSRN 3970746.")
    }
  } else if (connectedness=="Joint") {
    dca = JointConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore)
    if (model=="VAR") {
      message("The VAR joint connectedness approach is implemented according to:\n Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.")
    }
  } else if (connectedness=="Extended Joint") {
    dca = ExtendedJointConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore)
    if (model=="VAR" || model=="TVP-VAR") {
      message("The VAR and TVP-VAR extended joint connectedness approach is implemented according to:\n Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy, 73, 102219.")
    } else if (model=="QVAR") {
      message("The QVAR extended joint connectedness approach is implemented according to:\n Cunado, J, Chatziantoniou, I., Gabauer, D., Hardik, M., & de Garcia, F.P. (2022). Dynamic spillovers across precious metals and energy realized volatilities: Evidence from quantile extended joint connectedness measures.")
    }
  } else {
    stop("Connectedness approach does not exist")
  }
  dca
}







#' @title Time-varying parameter vector autoregression
#' @description Estimate TVP-VAR model
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param prior List of prior VAR coefficients and variance-covariance matrix
#' @param l forgetting factors (kappa1, kappa2)
#' @param configuration model configuration
#' @return Estimate TVP-VAR model
#' @examples
#' \donttest{
#' data(dy2012)
#' prior = BayesPrior(dy2012, nlag=1)
#' fit = TVPVAR(dy2012, configuration=list(nlag=1, prior=prior, l=c(0.99,0.99)))
#' }
#' @importFrom MASS ginv
#' @importFrom stats cov
#' @importFrom progress progress_bar
#' @references
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions. European Economic Review, 71, 101-116.\\
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management, 13(4), 84.
#' @author David Gabauer
#' @export
TVPVAR = function(x, configuration=list(l=c(0.99,0.99), nlag=1, prior=NULL)){
  l = as.numeric(configuration$l)
  nlag = configuration$nlag
  prior = configuration$prior
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  if (l[1] <= 0 || l[1] >= 1) {
    stop("kappa1 needs to be within 0 and 1")
  }
  if (l[2] <= 0 || l[2] >= 1) {
    stop("kappa2 needs to be within 0 and 1")
  }
  k = ncol(x)
  if (is.null(prior)) {
    prior = UninformativePrior(k, nlag)
  }
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  beta_0.mean = prior$bprior
  beta_0.var = prior$Vprior
  Q = prior$Q
  if (is.null(Q)) {
    Q = cov(x)
  }
  create_RHS_NI = function(templag, r, nlag, t){
    K = nlag*(r^2)
    x_t = matrix(0, (t-nlag)*r, K)
    for (i in 1:(t-nlag)){
      ztemp = NULL
      for (j in 1:nlag) {
        xtemp = templag[i,((j-1)*r+1):(j*r)]
        xtemp = t(kronecker(diag(r),xtemp))
        ztemp = cbind(ztemp, xtemp)
      }
      x_t[((i-1)*r+1):(i*r),] = ztemp
    }
    return=list(x_t=x_t, K=K)
  }
  x = scale(x,TRUE,FALSE)
  YX = cbind(x,x)
  r = p = n = ncol(x)
  m = nlag*(r^2)
  k = nlag*r
  t = nrow(x)
  q = n + p

  # Initialize matrices
  beta_0_prmean = beta_0.mean
  beta_0_prvar = beta_0.var

  beta_pred = matrix(0,m,t)
  beta_update = matrix(0,m,t)

  Rb_t = array(0,c(m,m,t))
  Sb_t = array(0,c(m,m,t))

  beta_t = array(0, c(k,k,t))
  Q_t = array(0, c(r,r,t), dimnames=list(NAMES, NAMES, as.character(zoo::index(x))))
  
  # Decay and forgetting factors
  l_2 = l[2]
  l_4 = l[1]

  # Define lags of the factors to be used in the state (VAR) equation
  yy = x[(nlag+1):t,]
  xx = embed(x, nlag+1)[,-c(1:r)]
  templag = embed(x, nlag+1)[,-c(1:r)]
  RHS1 = create_RHS_NI(templag,r,nlag,t)
  Flagtemp = RHS1$x_t
  Flag = rbind(matrix(0, k,m), Flagtemp)

  ###-----| 1. KALMAN FILTER
  pb = progress_bar$new(total=t)
  for (irep in 1:t) {
    #-----| Update the state covariances
    # 1. Get the variance of the factor
    # Update Q[t]
    if (irep==1) {
      Q_t[,,irep] = Q
    } else if (irep > 1) {
      if (irep <= (nlag+1)) {
        Gf_t = 0.1*(t(matrix(x[irep,],nrow=1))%*%(x[irep,]))
      } else {
        Gf_t = t(yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k])) %*% (yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k]))
      }
      Q_t[,,irep] = l_2*Q_t[,,(irep-1)] + (1-l_2)*Gf_t[1:r,1:r]
    }
    # -for beta
    if (irep <= (nlag+1)) {
      beta_pred[,irep] = beta_0_prmean
      beta_update[,irep] = beta_pred[,irep]
      Rb_t[,,irep] = beta_0_prvar
    } else if (irep > (nlag+1)) {
      beta_pred[,irep] = beta_update[,(irep-1)]
      Rb_t[,,irep] = (1/l_4)*Sb_t[,,(irep-1)]
    }

    # -for beta
    if (irep >= (nlag+1)) {
      # 2/ Update VAR coefficients conditional on Principal Componets estimates
      Rx = Rb_t[,,irep]%*%t(Flag[((irep-1)*r+1):(irep*r),])
      KV_b = Q_t[,,irep] + Flag[((irep-1)*r+1):(irep*r),]%*%Rx
      KG = Rx%*%ginv(KV_b)
      beta_update[,irep] = matrix(beta_pred[,irep], ncol=1) + (KG%*%(t(matrix(x[irep,], nrow=1))-Flag[((irep-1)*r+1):(irep*r),]%*%matrix(beta_pred[,irep], ncol=1)) )
      Sb_t[,,irep] = Rb_t[,,irep] - KG%*%(Flag[((irep-1)*r+1):(irep*r),]%*%Rb_t[,,irep])
    }

    # Assign coefficients
    bb = matrix(beta_update[,irep], ncol=1)
    splace = 0
    biga = matrix(0, r,r*nlag)
    for (ii in 1:nlag) {
      for (iii in 1:r) {
        biga[iii,((ii-1)*r+1):(ii*r)] = t(bb[(splace+1):((splace+r)),1])
        splace = splace + r
      }
    }
    B = rbind(biga, cbind(diag(r*(nlag-1)), matrix(0, nrow=r*(nlag-1), ncol=r)))

    if ((max(abs(eigen(B)$values))<=1)||(irep==1)){
      beta_t[,,irep] = B
    } else {
      beta_t[,,irep] = beta_t[,,(irep-1)]
      beta_update[,irep] = 0.99*beta_update[,(irep-1)]
    }
    pb$tick()
  }
  B_t = beta_t[1:ncol(Q_t),,]
  return = list(B_t=B_t, Q_t=Q_t)
}


#' @title Bayes Prior
#' @description Get Bayes prior
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param size Sample size used to calculate prior parameters
#' @return Get Bayes Prior
#' @examples
#' data(dy2012)
#' prior = BayesPrior(dy2012, nlag=1)
#' @references Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy. The Review of Economic Studies, 72(3), 821-852.
#' @author David Gabauer
#' @export
BayesPrior = function(x, size=NULL, nlag) {
  if (class(x)!="zoo") {
    stop("x needs to be a zoo matrix")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (is.null(size)) {
    size = nrow(x)
  }
  if ((size-nlag)<=0) {
    stop("size needs to be larger than nlag")
  }
  
  p = ncol(x)
  size = size - nlag
  yt = t(x[(nlag+1):(size+nlag),])
  m = p + nlag*(p^2)
  Zt = NULL
  for (i in (nlag+1):(size+nlag)) {
  
     print(paste0("第一个for","  ",i))
    ztemp = diag(p)
    for (j in 1:nlag) {
      xlag = x[(i-j),1:p]
      xtemp = matrix(0,p,p*p)
      for (jj in 1:p) {
        xtemp[jj,((jj-1)*p+1):(jj*p)] = xlag
      }
      ztemp = cbind(ztemp, xtemp)
    }
    Zt = rbind(Zt, ztemp)
  }
  vbar = matrix(0,m,m)
  xhy = matrix(0,m,1)
  for (i in 1:size) {
  print(paste0("第二个for","  ",i))
    zhat1 = Zt[((i-1)*p+1):(i*p),]
    vbar = vbar + t(zhat1)%*%zhat1
    xhy = xhy + t(zhat1)%*%as.matrix(yt[,i])
  }
  vbar = MASS::ginv(vbar)
  aols = vbar%*%xhy
  
  sse2 = matrix(0,p,p)
  for (i in 1:size) {
  print(paste0("第三个for","  ",i))
    zhat1 = Zt[((i-1)*p+1):(i*p),]
    sse2 = sse2 + (yt[,i] - zhat1%*%aols)%*%t(yt[,i] - zhat1%*%aols)
  }
  hbar = sse2/size
  return(list(bprior=aols[-c(1:p)],Vprior=vbar[-c(1:p),-c(1:p)],Q=hbar))
}



#' @title Aggregated Connectedness Measures
#' @description This function results in aggregated connectedness measures.
#' @param dca Dynamic connectedness object
#' @param groups List of at least two group vectors
#' @param start Start index
#' @param end End index
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Gabauer and Gupta (2018)
#' data("gg2018")
#' dca = ConnectednessApproach(gg2018, 
#'                             nlag=1, 
#'                             nfore=10, 
#'                             window.size=200,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.99, 
#'                             prior="BayesPrior")))
#' ac = AggregatedConnectedness(dca, groups=list("US"=c(1,2,3,4), "JP"=c(5,6,7,8)))
#' }
#' @references Chatziantoniou, I., Gabauer, D., & Stenfor, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.
#' @author David Gabauer
#' @export
AggregatedConnectedness = function(dca, groups, start=NULL, end=NULL) {
  corrected = dca$config$corrected
  message("Aggregated connectedness measures are introduced accoring to:\n Chatziantoniou, I., Gabauer, D., & Stenfor, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.")
  if (is.null(start)) {
    start = 1
  }
  if (is.null(end)) {
    end = dim(dca$CT)[3]
  }
  NAMES = dimnames(dca$NET)[[2]]
  k = length(NAMES)
  m = length(groups)
  CT = dca$CT
  t = dim(CT)[3]
  weights = NULL
  for (i in 1:m) {
    weights[i] = length(groups[i][[1]])
  }
  
  if (is.null(names(groups))) {
    NAMES_group = paste0("GROUP", 1:m)
  } else {
    NAMES_group = names(groups)
  }
  date = as.character(as.Date(dimnames(CT)[[3]]))
  
  if (length(groups) <= 1) {
    stop("groups need to consist of at least 2 vectors")
  }
  
  if (dca$config$approach == "Joint") {
    stop(paste("Aggregated connectedness measures are not implemented for", 
               dca$config$approach, "connectedness"))
  } else if (dca$config$approach == "Frequency") {
    mn = dim(CT)[4]
    TABLE = list()
    horizons = dimnames(CT)[[4]]
    TCI_ = TCI = array(NA, c(t,mn), dimnames=list(date,horizons))
    FROM = TO = NET = array(NA, c(t,m,mn), dimnames=list(date, NAMES_group,horizons))
    CT_ = NPDC = INFLUENCE = array(NA, c(m,m,t,mn), dimnames=list(NAMES_group, NAMES_group, date,horizons))
    for (jl in 1:mn) {
      for (il in 1:t) {
        ct0 = ct = CT[,,il,jl]
        for (i in 1:m) {
          for (j in 1:m) {
            if (i==j) {
              ct0[groups[i][[1]], groups[j][[1]]] = 0
            }
          }
        }
        ct1 = array(0, c(m, m), dimnames=list(NAMES_group, NAMES_group))
        for (i in 1:m) {
          for (j in 1:m) {
            ct1[i,j] = sum(ct0[groups[i][[1]], groups[j][[1]]])# / length(groups[j][[1]])
          }
        }
        for (i in 1:m) {
          ct1[i,i] = sum(ct[i,])-sum(ct1[i,])
        }
        CT_[,,il,jl] = ct1
        dca_ = ConnectednessTable(ct1)
        if (corrected) {
          TCI[il,jl] = dca_$cTCI
          TCI_[il,jl] = sum(dca_$TO * (k-weights)/(k-1))
        } else {
          TCI[il,jl] = dca_$TCI
          TCI_[il,jl] = sum(dca_$TO * (k-weights)/k)
        }
        TO[il,,jl] = dca_$TO
        FROM[il,,jl] = dca_$FROM
        NET[il,,jl] = dca_$NET
        NPDC[,,il,jl] = dca_$NPDC
      }
      TABLE[[jl]] = ConnectednessTable(CT_[,,,jl])$TABLE
    }
    names(TABLE) = horizons
  } else {
    TCI_ = TCI = array(NA, c(t,1), dimnames=list(date, "TCI"))
    FROM = TO = NET = array(NA, c(t,m), dimnames=list(date, NAMES_group))
    CT_ = PDC = INFLUENCE = array(NA, c(m,m,t), dimnames=list(NAMES_group, NAMES_group, date))
    for (il in 1:t) {
      ct0 = ct = CT[,,il]
      for (i in 1:m) {
        for (j in 1:m) {
          if (i==j) {
            ct0[groups[i][[1]], groups[j][[1]]] = 0
          }
        }
      }
      ct1 = array(0, c(m, m), dimnames=list(NAMES_group, NAMES_group))
      for (i in 1:m) {
        for (j in 1:m) {
          ct1[i,j] = sum(ct0[groups[i][[1]], groups[j][[1]]]) / length(groups[j][[1]])
        }
      }
      for (i in 1:m) {
        ct1[i,i] = 1-sum(ct1[i,])
      }
      CT_[,,il] = ct1
      dca_ = ConnectednessTable(ct1)
      if (corrected) {
        TCI[il, ] = dca_$cTCI
        TCI_[il,] = sum(dca_$TO * (k-weights)/(k-1))
      } else {
        TCI[il, ] = dca_$TCI
        TCI_[il,] = sum(dca_$TO * (k-weights)/k)
      }
      TO[il,] = dca_$TO
      FROM[il,] = dca_$FROM
      NET[il,] = dca_$NET
      NPDC[,,il] = dca_$NPDC
    }
    TABLE = ConnectednessTable(CT_)$TABLE
  }
  config = list(approach="Aggregated")
  return = list(TABLE=TABLE, TCI_ext=TCI_, TCI=TCI, 
                TO=TO, FROM=FROM, NPT=NULL, NET=NET, 
                NPDC=NPDC, INFLUENCE=NULL, PCI=NULL, config=config)
}



#' @title Bivariate DCC-GARCH
#' @description This function multiple Bivariate DCC-GARCH models that captures more accurately conditional covariances and correlations
#' @param x zoo dataset
#' @param spec A cGARCHspec A cGARCHspec object created by calling cgarchspec.
#' @param copula "mvnorm" or "mvt" (see, rmgarch package)
#' @param method "Kendall" or "ML" (see, rmgarch package)
#' @param transformation "parametric", "empirical" or "spd" (see, rmgarch package)
#' @param time.varying Boolean value to either choose DCC-GARCH or CCC-GARCH
#' @param asymmetric Whether to include an asymmetry term to the DCC model (thus estimating the aDCC).
#' @return Estimate Bivariate DCC-GARCH
#' @importFrom rmgarch cgarchspec
#' @importFrom rmgarch cgarchfit
#' @importFrom rmgarch rcor
#' @importFrom rmgarch rcov
#' @author David Gabauer
#' @export
BivariateDCCGARCH = function(x, spec, copula="mvt", method="Kendall", transformation="parametric", time.varying=TRUE, asymmetric=FALSE) {
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  t = nrow(x)
  k = ncol(x)
  NAMES = colnames(x)
  Z_t = NULL
  H_t = R_t = array(NA, c(k,k,t), dimnames=list(NAMES, NAMES, as.character(rownames(x))))
  for (i in 1:k) {
    for (j in 1:k) {
      if (i>j) {
        mgarch.spec = rmgarch::cgarchspec(uspec=multispec(c(spec[i], spec[j])), dccOrder=c(1,1), asymmetric=asymmetric,
                                          distribution.model=list(copula=copula, method=method, time.varying=time.varying, transformation=transformation))
        copula_fit = rmgarch::cgarchfit(mgarch.spec, data=x[,c(i,j)], solver=c("hybrid", "solnp"), fit.control=list(eval.se=FALSE))
        r = rcor(copula_fit)
        R_t[c(i,j),c(i,j),] = r
        h = rcov(copula_fit)
        H_t[c(i,j),c(i,j),] = h
      }
    }
    Z_t = cbind(Z_t, copula_fit@mfit$Z[,1])
  }
  return = list(H_t=H_t, R_t=R_t, Z_t=Z_t)
}



#' @title Kroner and Ng (1998) optimal bivariate portfolio weights
#' @description This function calculates the optimal portfolio weights according to Kroner and Ng (1998)
#' @param x zoo return matrix (in percentage)
#' @param H Residual variance-covariance, correlation or pairwise connectedness matrix
#' @param method Cumulative sum or cumulative product
#' @param long Allow only long portfolio position
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @return Get bivariate portfolio weights
#' @importFrom stats var.test
#' @importFrom onewaytests bf.test
#' @importFrom onewaytests homog.test
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' bpw = BivariatePortfolio(g2020, fit$Q, method="cumsum", statistics="Fisher")
#' bpw$TABLE
#' @references
#' Kroner, K. F., & Ng, V. K. (1998). Modeling asymmetric comovements of asset returns. The Review of Financial Studies, 11(4), 817-844.\\
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.\\
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
BivariatePortfolio = function(x, H, method=c("cumsum", "cumprod"), long=TRUE, statistics=c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), digit=2) {
  message("The optimal bivariate portfolios are computed according to:\n Kroner, K. F., & Ng, V. K. (1998). Modeling asymmetric comovements of asset returns. The Review of Financial Studies, 11(4), 817-844.")
  message("Hedging effectiveness is calculated according to:\n Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.")
  message("Statistics of the hedging effectiveness measure are implemented according to:\n Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.")
  
  method = match.arg(method)
  statistics = match.arg(statistics)
  x = x / 100
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  t = nrow(x)
  date = as.character(rownames(x))
  NAMES = colnames(x)

  summary = NULL
  portfolio_weights = array(0.5, c(k, k, t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:k) {
    for (j in 1:k) {
      pw = (H[j,j,]-H[i,j,]) / (H[i,i,]-2*H[i,j,]+H[j,j,])
      pw[which(is.na(pw))] = 0.5
      if (long) {
        pw = ifelse(pw>1,1,pw)
        pw = ifelse(pw<0,0,pw)
      }
      portfolio_weights[j,i,] = pw
      portfolio_weights[i,j,] = 1 - pw
      x_ = as.matrix(portfolio_weights[j,i,])
      summary_ = matrix(NA, nrow=ncol(x_), ncol=4)
      for (ij in 1:ncol(x_)){
        summary_[ij,] = matrix(c(mean(x_[,ij]), stats::sd(x_[,ij]), stats::quantile(x_[,ij],0.05), stats::quantile(x_[,ij],0.95)), nrow=1)
      }
      colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
      rownames(summary_) = paste0(NAMES[i], "/", NAMES[j])
      summary = rbind(summary, summary_)
    }
  }

  pvalue = HE = array(NA, c(k,k), dimnames=list(NAMES, NAMES))
  portfolio_return = cumulative_portfolio_return = cumulative_asset_return = array(NA,c(k,k,t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:k) {
    for (j in 1:k) {
      portfolio_return[j,i,] = portfolio_weights[j,i,]*x[,i] + (1-portfolio_weights[j,i,])*x[,j]
      HE[j,i] = 1 - var(portfolio_return[j,i,])/var(x[,i])
      df = rbind(data.frame(val=portfolio_return[j,i,], group="A"), data.frame(val=x[,i], group="B"))
      if (statistics=="Fisher") {
        pvalue[i,j] = stats::var.test(x=portfolio_return[j,i,],y=x[,i],ratio=1)$p.value
      } else if (statistics=="Bartlett") {
        pvalue[i,j] = onewaytests::homog.test(val~as.character(group), data=df, method="Bartlett", verbose=F)$p.value
      } else if (statistics=="Fligner-Killeen") {
        pvalue[i,j] = onewaytests::homog.test(val~as.character(group), data=df, method="Fligner", verbose=F)$p.value
      } else if (statistics=="Levene") {
        pvalue[i,j] = onewaytests::homog.test(val~as.character(group), data=df, method="Levene", verbose=F)$p.value
      } else if (statistics=="Brown-Forsythe") {
        pvalue[i,j] = onewaytests::bf.test(val~as.character(group), data=df, verbose=F)$p.value
      } else {
        stop("No valid hedging effectiveness statistics have been chosen.")
      }

      if (method=="cumsum") {
        cumulative_asset_return[j,i,] = cumsum(x[,i])
        cumulative_portfolio_return[j,i,] = cumsum(portfolio_return[j,i,])
      } else if (method=="cumprod") {
        cumulative_asset_return[j,i,] = cumprod(1+x[,i])
        cumulative_portfolio_return[j,i,] = cumprod(1+portfolio_return[j,i,])
      }
    }
  }
  TABLE = cbind(summary,c(HE),c(pvalue))
  TABLE = TABLE[-which(TABLE[,1]==0.5),]
  colnames(TABLE) = c("Mean","Std.Dev.","5%","95%","HE","p-value")

  return = list(TABLE=format(round(TABLE,digit),nsmall=digit), portfolio_weights=portfolio_weights, portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return)
}


#' @title ConditionalConnectedness
#' @description This function computes the conditional connectedness measures.
#' @param dca Dynamic connectedness object
#' @param group Group vector
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Chatzianzoniou, Gabauer and Stenfors (2022)
#' data(cgs2022)
#' dca = ConnectednessApproach(cgs2022, 
#'                             nlag=1, 
#'                             nfore=10, 
#'                             window.size=250,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.99, 
#'                             prior="BayesPrior")))
#' cc = ConditionalConnectedness(dca, group=c(1,4,7,10,13,16))
#' }
#' @references Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.
#' @author David Gabauer
#' @export
ConditionalConnectedness = function(dca, group=c(1,2,3), start=NULL, end=NULL) {
  corrected = dca$config$corrected
  message("Conditional connectedness measures are implemented according to:\n Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.")
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Conditional connectedness measures are not implemented for",dca$config$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = dca$CT[group,group,start:end,drop=FALSE]
    k = length(group)
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    t = length(date)
    
    TCI = array(NA, c(t,1), dimnames=list(date,"TCI"))
    NPT = TO = FROM = NET = array(NA, c(t,k), dimnames=list(date,NAMES))
    INFLUENCE = PCI = FEVD = NPDC = array(NA, c(k,k,t), dimnames=list(NAMES,NAMES,date))
    for (i in 1:t) {
      cc = ConnectednessTable(ct[,,i]/rowSums(ct[,,i]))
      FEVD[,,i] = cc$FEVD
      TO[i,] = cc$TO
      FROM[i,] = cc$FROM
      NET[i,] = cc$NET
      NPT[i,] = cc$NPT
      NPDC[,,i] = cc$NPDC
      PCI[,,i] = cc$PCI
      INFLUENCE[,,i] = cc$INFLUENCE
      if (corrected) {
        TCI[i,] = cc$cTCI
      } else {
        TCI[i,] = cc$TCI
      }
    }
    TABLE = ConnectednessTable(FEVD/100)$TABLE
    config = list(approach="Conditional")
    return = list(TABLE=TABLE, FEVD=FEVD, TCI=TCI, NET=NET, TO=TO, FROM=FROM, 
                  NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
  }
}

#' @title Partial Conditional Correlations
#' @description Compute partial conditional correlations
#' @param Q Variance-covariance matrix of dimension
#' @return Get partial conditional correlations
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' pcc = ConditionalCorrelation(fit$Q)
#' @author David Gabauer
#' @export
ConditionalCorrelation = function (Q) {
  if (length(dim(Q))<=1) {
    stop("Q needs to be at least a 2-dimensional matrix")
  }
  k = dim(Q)[1]
  NAMES = colnames(Q)
  if (length(dim(Q))==2) {
    Q = array(Q, c(k,k,1), dimnames=list(NAMES,NAMES))
  }
  R = Q
  for (i in 1:k) {
    for (j in 1:k) {
      R[i,j,] = Q[i,j,] / (sqrt(Q[i,i,])*sqrt(Q[j,j,]))
    }
  }
  return(R)
}

.onAttach <- 
  function(libname, pkgname) {
    packageStartupMessage("\nPlease cite as: \n")
    packageStartupMessage(" Gabauer, David (2022). ConnectednessApproach.")
    packageStartupMessage(" R package version 1.0.0. https://CRAN.R-project.org/package=ConnectednessApproach \n")
  }
  
  
#' @title Connectedness table
#' @description This function provides standard connectedness table.
#' @param FEVD Forecast error variance decomposition
#' @param digit Number of decimal places
#' @return Get connectedness table
#' @examples
#' \donttest{
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' fevd = FEVD(Phi=fit$B, Sigma=fit$Q, nfore=10, type="time", generalized=TRUE)$FEVD
#' dca = ConnectednessTable(fevd)
#' }
#' @references
#' Chatziantoniou, I., & Gabauer, D. (2021). EMU risk-synchronisation and financial fragility through the prism of dynamic connectedness. The Quarterly Review of Economics and Finance, 79, 1-14.\\
#' Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management, 60, 100680.
#' @export
ConnectednessTable = function(FEVD, digit=2) {
  if (length(dim(FEVD))<=1) {
    stop("FEVD needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(FEVD)
  k = dim(FEVD)[1]
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  CT = apply(FEVD,1:2,mean)*100 # spillover from others to one specific
  OWN = diag(diag(CT))
  TO = colSums(CT-OWN)
  FROM = rowSums(CT-OWN)
  NET = TO-FROM
  TCI = mean(TO) # Total connectedness index
  cTCI = TCI*k/(k-1)
  NPDC = CT-t(CT) # Net pairwise dynamic connectedness
  NPT = rowSums(NPDC<0)
  INFLUENCE = 100*abs(NPDC/t(t(CT)+CT))
  table = format(round(cbind(CT,FROM),digit),nsmall=digit)
  to = c(format(round(c(TO,sum(TO)),digit),nsmall=digit))
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "cTCI/TCI")
  tci = paste0(format(round(cTCI,digit),nsmall=digit),"/",format(round(TCI,digit),nsmall=digit))
  net = c(format(round(NET,digit),nsmall=digit))
  net = c(net, tci) 
  npt = c(format(round(NPT,digit),nsmall=digit), "")
  
  TABLE = rbind(table,to,inc,net,npt)
  colnames(TABLE) = c(NAMES,"FROM")
  rownames(TABLE) = c(NAMES,"TO","Inc.Own","NET","NPT")
  PCI = matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      PCI[i,j] = 200*(CT[i,j]+CT[j,i])/(CT[i,i]+CT[i,j]+CT[j,i]+CT[j,j])
    }
  }
  return = list(FEVD=CT, TCI=TCI, cTCI=cTCI, PCI=PCI,
                TO=TO, FROM=FROM, NET=NET, NPDC=NPDC, TABLE=TABLE,
                NPT=NPT, INFLUENCE=INFLUENCE)
}



#' @title DCC-GARCH selection specification
#' @description This function calculates the optimal DCC-GARCH specification
#' @param x zoo data matrix
#' @param distributions Vector of distributions
#' @param models Vector of GARCH models
#' @param ar AR(p)
#' @param ma MA(q)
#' @param prob The quantile (coverage) used for the VaR.
#' @param conf.level Confidence level of VaR test statistics
#' @param lag Lag length of weighted Portmanteau statistics
#' @return Get best DCC-GARCH
#' @importFrom stats bartlett.test coef fitted fligner.test integrate qnorm quantile residuals sd sigma var.test
#' @references
#' Ghalanos, A. (2014). rugarch: Univariate GARCH models, R package version 1.3-3.
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.
#' @author David Gabauer
#' @export
DCCGARCHselection = function(x, distributions=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH"), prob=0.05, conf.level=0.90, lag=20, ar=0, ma=0) {
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  for (i in 1:k) {
    print(colnames(x)[i])
    sel = GARCHselection(x=x[,i], distributions=distributions, models=models, 
                         prob=prob, conf.level=conf.level, lag=lag, ar=ar, ma=ma)
    if (i==1) {
      mspec = sel$best_ugarch
    } else {
      mspec = c(mspec, sel$best_ugarch)
    }
  }
  mspec
}


#' @title Elastic Net vector autoregression
#' @description Estimation of a VAR using equation-by-equation LASSO, Ridge or Elastic Net regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param nfolds N-fold cross validation
#' @param loss Loss function
#' @param alpha LASSO is alpha equal 1 and Ridge if alpha equal 0
#' @param delta_alpha Steps between 0 and 1. If alpha is NULL alpha is estimated based upon loss and nfolds
#' @param configuration Model configuration
#' @return Estimate VAR model
#' @examples
#' \donttest{
#' data(dy2012)
#' fit = ElasticNetVAR(dy2012, configuration=list(nlag=1, alpha=1, nfolds=10, loss="mae"))
#' }
#' @import glmnet
#' @importFrom stats predict
#' @references
#' Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor, J., & Tibshirani, R. J. (2012). Strong rules for discarding predictors in lasso‐type problems. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 74(2), 245-266.\\
#' Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1), 55-67.\\
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.\\
#' Demirer, M., Diebold, F. X., Liu, L., & Yilmaz, K. (2018). Estimating global bank network connectedness. Journal of Applied Econometrics, 33(1), 1-15.\\
#' Gabauer, D., Gupta, R., Marfatia, H., & Miller, S. M. (2020). Estimating US Housing Price Network Connectedness: Evidence from Dynamic Elastic Net, Lasso, and Ridge Vector Autoregressive Models. Lasso, and Ridge Vector Autoregressive Models (July 26, 2020).
#' @author David Gabauer
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet predict.glmnet
#' @export
ElasticNetVAR = function(x, configuration=list(nlag=1, nfolds=10, loss="mae", alpha=NULL, delta_alpha=0.1)) {
  nlag = configuration$nlag
  alpha = configuration$alpha
  nfolds = configuration$nfolds
  n_alpha = configuration$delta_alpha
  loss = configuration$loss
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (nfolds<=0) {
    stop("nfolds needs to be a positive integer")
  }
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  NAMES = colnames(x)
  k = ncol(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  if (is.null(alpha) ) {
    alpha = seq(0, 1, n_alpha)
  }
  alpha_ = Res = B = NULL
  for (i in 1:k) {

  
    MAE = NULL
    for (j in 1:length(alpha)) {
      z = embed(x, nlag+1)
      X = z[,-c(1:k)]
      y = z[,i]
      fit = glmnet::cv.glmnet(X, y, alpha=alpha[j], type.measure=loss, nfolds=nfolds)
      y_pred = predict(fit, s=fit$lambda.min, newx=X)
      MAE[j] = mean(abs(y - y_pred))
    }
    fit = glmnet::cv.glmnet(X, y, alpha=alpha[which(MAE==min(MAE))], type.measure=loss, nfolds=nfolds)
    y_pred = predict(fit, s=fit$lambda.min, newx=X)
    Res = cbind(Res, y - y_pred)
    b = predict(fit, type="coefficients", s=fit$lambda.min)[-1]
    B = rbind(B, b)
    alpha_[i] = alpha[which(MAE==min(MAE))]
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(as.character(zoo::index(x)),1)))
  results = list(B=B, Q=Q, alpha=alpha_)
}


#' @title Exclusive Connectedness Measures
#' @description This function results in exclusive connectedness measures
#' @param dca Dynamic connectedness object
#' @param group Vector of group indices
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Chatziantoniou, et al. (2022)
#' data("cegg2022")
#' dca = ConnectednessApproach(cegg2022,
#'                             nlag=1,
#'                             nfore=20,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             corrected=TRUE,
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, 
#'                             kappa2=0.99, prior="BayesPrior")))
#' exc = ExclusiveConnectedness(dca, group=c(1,2,3))
#' }
#' @references Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.
#' @author David Gabauer
#' @export
ExclusiveConnectedness = function(dca, group=c(1,2), start=NULL, end=NULL) {
  message("Partial connectedness measures are implemented according to:\n Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.")
  corrected = dca$config$corrected
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Partial connectedness measures are not implemented for",dca$config$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = dca$CT[,,start:end]
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    k = dim(ct)[1]
    t = dim(ct)[3]
    
    CT = ct
    for (i in group) {
      CT[,i,] = ct[,i,]*0
      CT[i,,] = ct[i,,]*0
    }
  
    TCI = array(NA, c(t,1), dimnames=list(as.character(date), "TCI"))
    NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
    NPDC = PCI = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
    for (i in 1:t) {
      dca_ = ConnectednessTable(CT[,,i])
      NPDC[,,i] = dca_$NPDC
      PCI[,,i] = dca_$PCI
      infl = dca_$INFLUENCE
      infl[which(is.nan(infl), arr.ind=TRUE)] = 0
      INFLUENCE[,,i] = infl
      TO[i,] = dca_$TO
      FROM[i,] = dca_$FROM
      NET[i,] = dca_$NET
      NPT[i,] = dca_$NPT
      if (corrected) {
        TCI[i,] = dca_$cTCI
      } else {
        TCI[i,] = dca_$TCI
      }
    }
    TABLE = ConnectednessTable(CT)$TABLE
    config = list(approach="Exclusive")
    return = list(TABLE=TABLE, TCI=TCI, NET=NET, TO=TO, FROM=FROM, NPT=NPT,
                  NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
  }
}


#' @title Balcilar et al. (2021) extended joint connectedness approach
#' @description This function provides extended joint connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Balcilar et al. (2021)
#' data("bgu2021")
#' prior = MinnesotaPrior(0.1, k=ncol(bgu2021), nlag=1)
#' fit = TVPVAR(bgu2021, configuration=list(l=c(0.99,0.99), nlag=1, prior=prior))
#' dca = ExtendedJointConnectedness(Phi=fit$B_t, Sigma=fit$Q_t, nfore=20)
#' dca$TABLE
#' }
#' @references
#' Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy, 73, 102219.
#' @author David Gabauer
#' @export
ExtendedJointConnectedness = function(Phi, Sigma, nfore) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(Sigma)
  if (length(dim(Phi))==2) {
    Phi = array(Phi, c(nrow(Phi),ncol(Phi),1))
  }
  if (length(dim(Sigma))==2) {
    Sigma = array(Sigma, c(nrow(Sigma),ncol(Sigma),1))
  }

  k = ncol(Sigma)
  t = dim(Sigma)[3]

  if (is.null(NAMES)) {
    NAMES = 1:k
  }

  date = dimnames(Sigma)[[3]]
  TCI = array(NA, c(t,1), dimnames=list(as.character(date), "TCI"))
  NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(as.character(date), NAMES))
  CT = PCI = NPDC = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, as.character(date)))

  pb = progress_bar$new(total=t)
  for (ij in 1:t) {
    # calculate the gFEVD
    gSOT = 100*FEVD(Phi[,,ij], Sigma[,,ij], nfore=nfore, type="time",
                    generalized=TRUE)$FEVD
    gSOI = mean(rowSums(gSOT-diag(diag(gSOT))))

    # calculate Xi (the forecast error covariance matrix)
    A = Wold(Phi[,,ij], nfore)  # the VMA coefficient matrices
    Xi_h = array(0,dim=c(k,k,nfore))
    for (h in 1:nfore) {
      Xi_h[,,h] = A[,,h]%*%Sigma[,,ij]%*%t(A[,,h]) # calculated Xi at each h
    }
    Xi = rowSums(Xi_h, dims=2) # sum them along THIRD dimension to form Xi  (note: because this is a row sum, dims=2, actually sums along the third dimension)
    I_K = diag(1,nrow=k,ncol=k)

    # Calculate the elimination matrices.
    M = array(0,dim=c(k,k-1,k))
    for (i in 1:k){
      M[,,i] = I_K[,-i] # calculate the elimination matrices
    }
    S_jnt_numerator_h = array(0,dim=c(k,nfore))
    for (i in 1:k) {
      for (h in 1:nfore){
        S_jnt_numerator_h[i,h] = I_K[i,]%*%A[,,h]%*%Sigma[,,ij]%*%M[,,i]%*%(ginv(t(M[,,i])%*%Sigma[,,ij]%*%M[,,i]))%*%t(M[,,i])%*%Sigma[,,ij]%*%t(A[,,h])%*%I_K[,i] #calculate the numerator of S_jnt at each h
      }
    }

    S_jnt_numerator = array(0,dim=c(k))
    for (i in 1:k) {
      S_jnt_numerator[i] = sum(S_jnt_numerator_h[i,]) # calculate the numerator of j_jnt  (sum over h)
    }

    S_jnt=array(0,dim=c(k))
    for (i in 1:k) {
      S_jnt[i] = (100*S_jnt_numerator[i])/Xi[i,i]
    }

    # calculate the joint spillover index (jSOI)
    gSOT_diag = gSOT
    diag(gSOT_diag) = 0
    jSOI = mean(S_jnt)
    lambda = S_jnt / apply(gSOT_diag, 1, sum)
    jSOT = gSOT
    colnames(jSOT)=rownames(jSOT)=NAMES
    for (i in 1:k) {
      jSOT[i,] = gSOT[i,]*lambda[i]
    }
    jSOT_diag = jSOT
    diag(jSOT_diag) = 0
    from_jnt = rowSums(jSOT_diag)
    to_jnt = colSums(jSOT_diag)
    jSOI = mean(to_jnt)
    diag(jSOT_diag) = 100 - from_jnt

    dca = ConnectednessTable(jSOT_diag/100)
    CT[,,ij] = dca$FEVD
    TO[ij,] = dca$TO
    FROM[ij,] = dca$FROM
    NET[ij,] = dca$NET
    NPDC[,,ij] = dca$NPDC
    TCI[ij,] = dca$TCI
    PCI[,,ij] = dca$PCI
    NPT[ij,] = dca$NPT
    INFLUENCE[,,ij] = dca$INFLUENCE
    pb$tick()
  }

  TABLE = ConnectednessTable(CT/100)$TABLE
  config = list(nfore=nfore, approach="Extended Joint", generalized=TRUE, corrected=FALSE)
  return = list(TABLE=TABLE, CT=CT/100, TCI=TCI, TO=TO, FROM=FROM,
                NET=NET, NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
}


#' @title External Connectedness Measures
#' @description This function provides external connectedness measures
#' @param dca Dynamic connectedness object
#' @param groups List of at least two group vectors
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("gg2018")
#' dca = ConnectednessApproach(gg2018, model="TVP-VAR",
#'                             connectedness="Time",
#'                             nlag=1, nfore=10, window.size=200,
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, 
#'                             kappa2=0.99, prior="BayesPrior")))
#' ext = ExternalConnectedness(dca, groups=list("US"=c(1,2,3,4), "JP"=c(5,6,7,8)))
#' }
#' @references Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.
#' @author David Gabauer
#' @export
ExternalConnectedness = function(dca, groups=list(c(1), c(2:ncol(dca$NET))), start=NULL, end=NULL) {
  message("The decomposed connectedness measures are implemented according to:\n Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.")
  corrected = dca$config$corrected
  if (length(groups)<=1) {
    stop("groups need to consist of at least 2 vectors")
  }
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Decomposed connectedness measures are not implemented for",dca$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = 100*dca$CT[,,start:end]
    NAMES = colnames(ct)
    k = dim(ct)[2]
    if (length(dim(ct))==2) {
      ct = array(ct, c(k,k,1),dimnames=list(NAMES,NAMES))
    }
    ct_inter = ct_wo = ct
    date = as.character(dimnames(ct)[[3]])
    t = dim(ct)[3]
    
    m = length(groups)
    NAMES_group = names(groups)
    if (is.null(NAMES_group)) {
      NAMES_group = paste0("GROUP", 1:m)
    }
    
    for (i in 1:m) {
      group_1 = groups[[i]]
      ct_wo[group_1,group_1,] = 0
    }
    
    TCI_wo = array(NA, c(t, 1), dimnames=list(date, c("TCI")))
    INFLUENCE_wo = PCI_wo = NPDC_wo = array(NA, c(k, k, t), dimnames=list(NAMES,NAMES,date))
    TO_wo = FROM_wo = NET_wo = NPT_wo = array(NA, c(t, k), dimnames=list(date, NAMES))
    for (i in 1:t) {
      dca_ = ConnectednessTable(ct_wo[,,i]/100)
      TO_wo[i,] = dca_$TO
      FROM_wo[i,] = dca_$FROM
      NET_wo[i,] = dca_$NET
      NPT_wo[i,] = dca_$NPT
      NPDC_wo[,,i] = dca_$NPDC
      PCI_wo[,,i] = dca_$PCI
      infl = dca_$INFLUENCE
      infl[which(is.nan(infl), arr.ind=TRUE)] = 0
      INFLUENCE_wo[,,i] = infl
      if (corrected) {
        TCI_wo[i,] = dca_$cTCI
      } else {
        TCI_wo[i,] = dca_$TCI
      }
    }
    if (corrected) {
      m_ = (k-1)
    } else {
      m_ = k
    }
    
    TCI_group = array(NA, c(t,m), dimnames=list(date, NAMES_group))
    for (i in 1:m) {
      group = groups[i][[1]]
      TCI_group[,i] = rowSums(TO_wo[,group,drop=FALSE])/m_
    }
    
    TABLE = ConnectednessTable(ct_wo/100)$TABLE
    config = list(approach="External")
    return = list(TABLE=TABLE, gTCI=TCI_group, TCI=TCI_wo, TO=TO_wo, FROM=FROM_wo, NPT=NPT_wo,
                  NET=NET_wo, NPDC=NPDC_wo, PCI=PCI_wo, INFLUENCE=INFLUENCE_wo, config=config)
  }
}


#' @title Forecast error variance decomposition
#' @description This function computes the orthogonalized/generalized forecast error variance decomposition
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param type Time or Frequency connectedness approach
#' @param range Partition range for frequency approach only.
#' @param generalized Generalized or orthogonalized FEVD
#' @param orth Orthogonalized shocks
#' @return Orthogonalized/generalized time/frequency forecast error variance decomposition
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' fevd = FEVD(Phi=fit$B, Sigma=fit$Q, nfore=10, type="time", generalized=TRUE)$FEVD
#' @references
#' Stiassny, A. (1996). A spectral decomposition for structural VAR models. Empirical Economics, 21(4), 535-555.\\
#' Koop, G., Pesaran, M. H., & Potter, S. M. (1996). Impulse response analysis in nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.\\
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. Economics Letters, 58(1), 17-29.
#' @importFrom stats fft
#' @export
FEVD = function (Phi, Sigma, nfore=100, type=c("time","frequency"), generalized=TRUE, orth=FALSE, range=NULL) {
  type = match.arg(type)
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  
  if (length(dim(Sigma))>2) {
    Sigma = Sigma[,,1]
  }
  irf_list = IRF(Phi=Phi, Sigma=Sigma, nfore=nfore, orth=orth)
  irf = irf_list$irf
  if (type=="time") {
    Phi = lapply(1:nfore, function(j) sapply(irf, function(i) i[j,]))
    # Phi is actually PSI here, which means the Moving Average coefficient matrics 
    denom = diag(Reduce("+", lapply(Phi, function(i) i %*% Sigma %*% t(i))))
    # denom: the total forecast error variance of each variable
    if (generalized) {
      enum = Reduce("+", lapply(Phi, function(i) (i %*% Sigma)^2))
      tab = sapply(1:nrow(enum), function(j) enum[j,]/(denom[j] * diag(Sigma)))
      FEVD = t(apply(tab, 2, function(i) i/sum(i)))
    } else {
      enum = Reduce("+", lapply(Phi, function(i) (chol(Sigma) %*% t(i))^2))
      FEVD = t(sapply(1:ncol(enum), function(i) enum[,i]/denom[i]))
    }
  } else if (type=="frequency") {
    fftir = lapply(irf, function(i) apply(i, 2, fft))
    fftir = lapply(1:(nfore+1), function(j) sapply(fftir, function(i) i[j,]))
    denom = diag(Re(Reduce("+", lapply(fftir, function(i) i %*% Sigma %*% t(Conj(i))/nfore)[range])))
    if (generalized) {
      enum = lapply(fftir, function(i) (abs(i %*% Sigma))^2/(nfore + 1))
      tab = lapply(enum, function(i) sapply(1:nrow(i), function(j) i[j, ]/(denom[j] * diag(Sigma))))
      tot = apply(Reduce("+", tab[range]), 2, sum)
      FEVD = lapply(tab, function(i) t(i)/tot)
    } else {
      enum = lapply(fftir, function(i) (abs(i %*% t(chol(Sigma))))^2/(nfore + 1))
      FEVD = lapply(enum, function(i) t(sapply(1:nrow(i), function(j) i[j,]/(denom[j]))))
    }
  }
  return = list(IRF=irf, FEVD=FEVD)
}

#' @title Baruník and Křehlík (2018) frequency connectedness approach
#' @description This function calculates the Baruník and Křehlík (2018) frequency connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param partition Frequency spectrum
#' @param generalized Orthorgonalized/generalized FEVD
#' @param scenario ABS or WTH
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param orth Orthorgonalized shocks
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("dy2012")
#' partition = c(pi+0.00001, pi/4, 0)
#' fit = VAR(dy2012, configuration=list(nlag=4))
#' dca = FrequencyConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=100, partition=partition)
#' }
#' @import frequencyConnectedness
#' @references
#' Baruník, J., & Křehlík, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics, 16(2), 271-296.
#' @author David Gabauer
#' @export
FrequencyConnectedness = function(Phi, Sigma, nfore=100, partition=c(pi,pi/2,0), generalized=TRUE, orth=FALSE, scenario="ABS", corrected=FALSE) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(Sigma)
  if (length(dim(Phi))==2) {
    Phi = array(Phi, c(nrow(Phi),ncol(Phi),1))
  }
  if (length(dim(Sigma))==2) {
    Sigma = array(Sigma, c(nrow(Sigma),ncol(Sigma),1))
  }

  k = dim(Sigma)[1]
  t = dim(Sigma)[3]
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  
  periods = round(pi/partition)
  period_names = NULL
  for (i in 1:(length(periods)-1)) {
    period_names = c(period_names, paste0(periods[i], "-", periods[i+1]))
  }
  period_names = c("Total",period_names)
  date = as.character(dimnames(Sigma)[[3]])
  interval = length(period_names)
  new_p = frequencyConnectedness::getPartition(partition, nfore)
  range = sort(unique(do.call(c, new_p)))
  
  TCI = array(0, c(t,interval), dimnames=list(date, period_names))
  PCI = INFLUENCE = CT = NPDC = array(0, c(k, k, t, interval), dimnames=list(NAMES, NAMES, date, period_names))
  NET = FROM = TO = array(0, c(t, k, interval), dimnames=list(date, NAMES, period_names))
  NPT = array(0, c(t, k, interval-1), dimnames=list(date, NAMES, period_names[-1]))
  PCI = INFLUENCE = array(0, c(k, k, t, interval-1), dimnames=list(NAMES, NAMES, date, period_names[-1]))
  pb = progress_bar$new(total=t)
  for (i in 1:t) {
    decomp = FEVD(Phi=Phi[,,i], Sigma=Sigma[,,i], nfore=nfore, generalized=generalized, type="frequency", range=range)$FEVD
    for (ij in 1:length(decomp)) {
      rownames(decomp[[ij]]) = colnames(decomp[[ij]]) = 1:ncol(Sigma)
    }
    tables = lapply(new_p, function(j) Reduce('+', decomp[j]))
    for (j in 2:(interval)) {
      if (scenario=="ABS") {
        dca = ConnectednessTable(tables[[j-1]])
        CT[,,i,j] = dca$FEVD
        TO[i,,j] = dca$TO
        FROM[i,,j] = dca$FROM
        NET[i,,j] = dca$NET
        NPDC[,,i,j] = dca$NPDC
        PCI[,,i,j-1] = dca$PCI
        INFLUENCE[,,i,j-1] = dca$INFLUENCE
        NPT[i,,j-1] = dca$NPT
        if (corrected) {
          TCI[i,j] = dca$cTCI
        } else {
          TCI[i,j] = dca$TCI
        }
      } else if (scenario=="WTH") {
        dca = ConnectednessTable(tables[[j-1]]/sum(sum(tables[[j-1]]))*k)
        CT[,,i,j] = dca$FEVD
        TO[i,,j] = dca$TO
        FROM[i,,j] = dca$FROM
        NET[i,,j] = dca$NET
        NPDC[,,i,j] = dca$NPDC
        PCI[,,i,j-1] = dca$PCI
        INFLUENCE[,,i,j-1] = dca$INFLUENCE
        NPT[i,,j-1] = dca$NPT
        if (corrected) {
          TCI[i,j] = dca$cTCI
        } else {
          TCI[i,j] = dca$TCI
        }
      }
    }
    pb$tick()
  }
  CT[,,,1] = apply(CT,1:3,sum)
  TCI[,1] = apply(TCI,1,sum)
  TO[,,1] = apply(TO,1:2,sum)
  FROM[,,1] = apply(FROM,1:2,sum)
  NET[,,1] = apply(NET,1:2,sum)
  NPDC[,,,1] = apply(NPDC,1:3,sum)

  TABLE = array(NA,c(k+4,k+1,interval), dimnames=list(c(NAMES, "TO", "Inc.Own", "Net", "NPDC"), c(NAMES, "FROM"), period_names))
  for (i in 1:interval) {
    TABLE[,,i] = ConnectednessTable(CT[,,,i]/100)$TABLE
  }
  config = list(partition=partition, nfore=nfore, generalized=generalized, orth=orth, scenario=scenario, corrected=corrected, approach="Frequency")
  return = list(TABLE=TABLE, CT=CT/100, TCI=TCI, TO=TO, FROM=FROM,
                NET=NET, NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
}


#' @title Univariate GARCH selection criterion
#' @description This function estimates and evaluates a combination of GARCH models with different distributions and suggests the best GARCH models among all alternatives given some test statistics
#' @param x zoo data matrix
#' @param distributions Vector of distributions
#' @param models Vector of GARCH models
#' @param ar AR(p)
#' @param ma MA(q)
#' @param prob The quantile (coverage) used for the VaR.
#' @param conf.level Confidence level of VaR test statistics
#' @param lag Lag length of weighted Portmanteau statistics
#' @return Get optimal univariate GARCH model specification
#' @importFrom stats bartlett.test coef fitted fligner.test integrate qnorm quantile residuals sd sigma var.test
#' @references
#' Ghalanos, A. (2014). rugarch: Univariate GARCH models, R package version 1.3-3.
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.
#' @author David Gabauer
#' @export
GARCHselection = function(x, distributions=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH"), prob=0.05, conf.level=0.90, lag=20, ar=0, ma=0) {
  message("A dynamic version of the optimal univariate GARCH selection procedure is implemented according to:\n Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.")
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  GARCH_IC = matrix(Inf, nrow=length(distributions), ncol=length(models))
  colnames(GARCH_IC) = models
  rownames(GARCH_IC) = distributions
  spec_list = list()
  table_list = list()
  for (i in 1:length(models)) {
    spec_list[[i]] = list()
    table_list[[i]] = list()
    for (j in 1:length(distributions)) {
      spec_list[[i]][[j]] = list()
    }
    names(spec_list[[i]]) = distributions
  }
  names(spec_list) = names(table_list) = models
  
  for (j in 1:length(models)) {
    message(paste0("-",models[j]))
    for (i in 1:length(distributions)) {
      message(paste0("--",distributions[i]))
      if (models[j] %in% c("AVGARCH","TGARCH","APARCH","NAGARCH","NGARCH","ALLGARCH")) {
        ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                          variance.model=list(model="fGARCH", submodel=models[j], garchOrder=c(1,1)),
                                          distribution.model=distributions[i])
      } else {
        ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                          variance.model=list(model=models[j], garchOrder=c(1,1)),
                                          distribution.model=distributions[i])
      }
      ugarch.fit = rugarch::ugarchfit(ugarch.spec, data=x, solver="hybrid", solver.list=list(outer.iter=10, inner.iter=1000, eval.se=FALSE, tol=1e-12))
      if (ugarch.fit@fit$convergence==0) {
        fit = GARCHtests(ugarch.fit, prob=prob, conf.level=conf.level, lag=lag)
        GARCH_IC[i,j] = fit$InformationCriterion
        spec_list[[models[j]]][[distributions[i]]] = ugarch.spec
        table_list[[j]][[distributions[i]]] = fit
      }
    }
  }
  GARCH_selection = which(GARCH_IC==min(GARCH_IC),arr.ind=TRUE)
  best_ugarch = spec_list[[GARCH_selection[2]]][[GARCH_selection[1]]]
  best_table = table_list[[GARCH_selection[2]]][[GARCH_selection[1]]]
  return = list(best_ugarch=best_ugarch, best_table=best_table, GARCH_IC=GARCH_IC, spec_list=spec_list, table_list=table_list)
}


#' @title Univariate GARCH test statistics
#' @description This function provides the results of multiple univariate GARCH test statistics
#' @param fit Fitted univariate GARCH
#' @param prob The quantile (coverage) used for the VaR.
#' @param conf.level Confidence level of VaR test statistics
#' @param lag Lag length of weighted Portmanteau statistics
#' @return Get best univariate GARCH
#' @references
#' Ghalanos, A. (2014). rugarch: Univariate GARCH models, R package version 1.3-3.
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.
#' @author David Gabauer
#' @export
#' @importFrom WeightedPortTest Weighted.Box.test
#' @importFrom rugarch qdist
#' @importFrom rugarch VaRDurTest
#' @importFrom stats qnorm
#' @importFrom utils tail
#' @importFrom xts as.xts
#' @importFrom zoo as.zoo
GARCHtests = function(fit, lag=20, prob=0.05, conf.level=0.90){
  distribution = fit@model$modeldesc$distribution
  model = fit@model$modeldesc$vmodel
  submodel = fit@model$modeldesc$vsubmodel
  ar = fit@model$modelinc['ar']
  ma = fit@model$modelinc['ma']
  if (is.null(submodel)) {
    ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                             variance.model=list(model=model, garchOrder=c(1,1)),
                             distribution.model=distribution)
  } else {
    ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                             variance.model=list(model=model, submodel=submodel, garchOrder=c(1,1)),
                             distribution.model=distribution)
  }
  x = xts::as.xts(zoo::as.zoo(fit@model$modeldata$data, order.by=fit@model$modeldata$index))
  t = length(x)

  VaR = rugarch::fitted(fit) + rugarch::sigma(fit)*rugarch::qdist(fit@model$modeldesc$distribution, p=prob, mu=fit@fit$matcoef[1,1], sigma=1, skew=ifelse(is.na(rugarch::coef(fit)["skew"]),0,rugarch::coef(fit)["skew"]), shape=rugarch::coef(fit)["shape"])
  var_test = rugarch::VaRTest(prob, actual=x, VaR=as.numeric(VaR), conf.level=conf.level)
  i = 99
  while (is.nan(var_test$uc.LRp)) {
    n = round(t*i/100)
    var_test = rugarch::VaRTest(prob, actual=utils::tail(x,n), VaR=utils::tail(as.numeric(VaR),n), conf.level=conf.level)
    i = i - 1
  }
  vardur_test = rugarch::VaRDurTest(prob, x, VaR,conf.level=conf.level)
  f = function(x) {
    rugarch::qdist(fit@model$modeldesc$distribution, p=x, mu=0,sigma=1,skew=ifelse(is.na(rugarch::coef(fit)["skew"]),0,rugarch::coef(fit)["skew"]), shape=rugarch::coef(fit)["shape"])
  }
  ES = rugarch::fitted(fit) + rugarch::sigma(fit)*stats::integrate(f,0,prob)$value/prob
  ES = rugarch::ESTest(prob, x, ES, VaR, boot=TRUE, n.boot=1000, conf.level=conf.level)
  sign.bias = rugarch::signbias(fit)[1,][1:2]
  warch = WeightedPortTest::Weighted.Box.test(rugarch::residuals(fit), type="Ljung-Box", lag=lag, sqrd.res=TRUE)

  statistics = c(sign.bias[[1]], warch$statistic, var_test$uc.LRstat, vardur_test$rLL, ES$actual.exceed/ES$expected.exceed)
  pvalues = c(sign.bias[2]$prob,warch$p.value, var_test$uc.LRp, round(ES$boot.p.value,3), vardur_test$LRp)
  if (length(which(is.nan(pvalues)))>0) {
    pvalues[which(is.nan(pvalues))] = 1.00
  }
  TABLE = rbind(statistics, pvalues)
  colnames(TABLE) = c("SignBias", paste0("WARCH(",lag,")"),"VaR","CVaR","VaR Dur.")

  qprob = stats::qnorm(1-prob)
  loss = sum(abs(fit@fit$robust.tval[-c(1:2)])<=qprob)
  if (is.na(fit@fit$robust.matcoef[1,2])==FALSE) {
    if ("skew" %in% rownames(fit@fit$robust.matcoef)) {
      upper = fit@fit$robust.matcoef["skew",1] + qprob*fit@fit$robust.matcoef["skew",2]
      lower = fit@fit$robust.matcoef["skew",1] - qprob*fit@fit$robust.matcoef["skew",2]
      if (upper>1 && lower<1) {
        loss = loss + 100
      }
    }
  }
  IC = -2*rugarch::likelihood(fit) + loss*log(t)
  IC = IC + sum(TABLE[2,]<0.10)*10^5
  IC = ifelse(is.na(IC), Inf, IC)
  
  return = list(InformationCriterion=IC, TABLE=TABLE)
}


#' @title Kroner and Sultan (1993) hedge ratios
#' @description This function calculates the hedge ratios of Kroner and Sultan (1993)
#' @param x zoo return matrix (in percentage)
#' @param H Residual variance-covariance, correlation or pairwise connectedness matrix
#' @param method Cumulative sum or cumulative product
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @return Get hedge ratios
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' hr = HedgeRatio(g2020, fit$Q)
#' hr$TABLE
#' @references
#' Kroner, K. F., & Sultan, J. (1993). Time-varying distributions and dynamic hedging with foreign currency futures. Journal of Financial and Quantitative Analysis, 28(4), 535-551.\\
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.\\
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
HedgeRatio = function(x, H, method=c("cumsum","cumprod"), statistics=c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), digit=2) {
  message("Hedge ratios are implemented according to:\n Kroner, K. F., & Sultan, J. (1993). Time-varying distributions and dynamic hedging with foreign currency futures. Journal of Financial and Quantitative Analysis, 28(4), 535-551.")
  message("Hedging effectiveness is calculated according to:\n Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.")
  message("Statistics of the hedging effectiveness measure are implemented according to:\n Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.")
  
  method = match.arg(method)
  statistics = match.arg(statistics)
  x = x / 100
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  t = nrow(x)
  date = as.character(rownames(x))
  NAMES = colnames(x)

  HR = array(NA, c(k, k, t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:k) {
    for (j in 1:k) {
      HR[i,j,] = H[i,j,] / H[j,j,]
    }
  }

  summary = NULL
  for (i in 1:k) {
    for (j in 1:k) {
      x_ = as.matrix(HR[i,j,])
      summary_ = matrix(NA, nrow=ncol(x_), ncol=4)
      for (ij in 1:ncol(x_)){
        summary_[ij,] = matrix(c(mean(x_[,ij]), stats::sd(x_[,ij]), stats::quantile(x_[,ij],0.05), stats::quantile(x_[,ij],0.95)), nrow=1)
      }
      colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
      rownames(summary_) = paste0(NAMES[i], "/", NAMES[j])
      summary = rbind(summary, summary_)
    }
  }
  colnames(summary) = c("Mean","Std.Dev.","5%","95%")

  HE = pvalue = array(NA,c(k,k), dimnames=list(NAMES,NAMES))
  portfolio_return = cumulative_portfolio_return = array(NA,c(k,k,t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:k) {
    for (j in 1:k) {
      portfolio_return[i,j,] = as.numeric(x[,i] - HR[i,j,]*x[,j])
      HE[j,i] = 1 - var(portfolio_return[i,j,])/var(x[,i])
      df = rbind(data.frame(val=portfolio_return[i,j,], group="A"), data.frame(val=x[,i], group="B"))
      if (statistics=="Fisher") {
        pvalue[i,] = var.test(x=portfolio_return[i,j,], y=x[,i],ratio=1)$p.value
      } else if (statistics=="Bartlett") {
        pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Bartlett", verbose=F)$p.value
      } else if (statistics=="Fligner-Killeen") {
        pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Fligner", verbose=F)$p.value
      } else if (statistics=="Levene") {
        pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Levene", verbose=F)$p.value
      } else if (statistics=="Brown-Forsythe") {
        pvalue[i,] = onewaytests::bf.test(val~as.character(group), data=df, verbose=F)$p.value
      } else {
        stop("No valid hedging effectiveness statistics have been chosen.")
      }

      if (method=="cumsum") {
        cumulative_portfolio_return[i,j,] = cumsum(portfolio_return[i,j,])
      } else if (method=="cumprod") {
        cumulative_portfolio_return[i,j,] = cumprod(1+portfolio_return[i,j,])-1
      }
    }
  }
  TABLE = cbind(summary,c(HE),c(pvalue))
  TABLE = TABLE[-which(TABLE[,1]==1),]
  colnames(TABLE)=c("Mean","Std.Dev.","5%","95%","HE","p-value")

  return = list(TABLE=format(round(TABLE,digit),nsmall=digit), hedge_ratio=HE, portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return)
}

#' @title Impulse response functions
#' @description This function calculates orthorgonalized/generalized impulse response functions of time or frequency domain.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual Variance-Covariance Matrix
#' @param nfore H-step ahead forecast horizon
#' @param orth Boolean
#' @return Orthorgonal/generalized time/frequency impulse response functions
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' irf = IRF(Phi=fit$B, Sigma=fit$Q, nfore=10, orth=TRUE)
#' @references
#' Stiassny, A. (1996). A spectral decomposition for structural VAR models. Empirical Economics, 21(4), 535-555.\\
#' Koop, G., Pesaran, M. H., & Potter, S. M. (1996). Impulse response analysis in nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.\\
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. Economics Letters, 58(1), 17-29.
#' @author David Gabauer
#' @export
IRF = function (Phi, Sigma, nfore=10, orth=TRUE) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  
  if (length(dim(Sigma))>2) {
    Sigma = Sigma[,,1]
  }
  p = 0
  k = 0
  if (length(Phi) > 0) { # Phi is a k-by-(lag*k) matrix
    k = dim(Phi)[1]
    k1 = dim(Phi)[2]
    p = floor(k1/k) # p = lag
  }
  if (is.null(Sigma)) {
    Sigma = diag(rep(1, k))
  }
  if (orth) {
    m1 = eigen(Sigma)
    v1 = sqrt(m1$values)
    vv = diag(v1)
    Pmtx = m1$vectors
    Sh = Pmtx %*% vv %*% t(Pmtx)
  }
  if (k < 1) {
    k = 1
  }
  PSI = diag(rep(1, k)) 
  # PSI is the combination of the N*N Moving Average coefficient matrices
  # where PSI_0 = I_N (Identity matrix), PSI_i = B_1 PSI_{i-1} + ... + B_p PSI_{i-p}
  # B is denoted as Phi here
  if (orth) {
    WGT = c(PSI %*% Sh)
  } else {
    WGT = c(PSI)
  }
  for (il in 1:nfore) {
    ilk = il * k
    tmp = matrix(0, k, k)
    if (p > 0) {
      iend = min(il, p)
      for (j in 1:iend) {
        jdx = (il - j)
        kdx = j * k
        tmp = tmp + Phi[, (kdx - k + 1):kdx] %*% PSI[,(jdx * k + 1):(jdx * k + k)]
      }
    }
    PSI = cbind(PSI, tmp)
    if (orth) {
      WGT = cbind(WGT, c(tmp %*% Sh))
    } else {
      WGT = cbind(WGT, c(tmp))
    }
  }
  wk1 = WGT
  for (i in 1:k^2) {
    wk1[i, ] = cumsum(WGT[i, ])
  }
  tdx = c(1:(nfore + 1)) - 1
  if (orth) {
    gmax = max(WGT)
    gmin = min(WGT)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
  } else {
    gmax = max(WGT)
    gmin = min(WGT)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
  }
  k = sqrt(nrow(WGT))
  irf = list()
  for (i in 1:k) {
    irf[[i]] = t(WGT[((i-1)*k+1):(i*k),])
  }
  return = list(irf=irf)
}


#' @title Inclusive Connectedness Measures
#' @description This function results in inclusive connectedness measures
#' @param dca Dynamic connectedness object
#' @param group Vector of group indices
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("cegg2022")
#' dca = ConnectednessApproach(cegg2022,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             nlag=1,
#'                             nfore=20,
#'                             corrected=TRUE,
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, 
#'                             kappa2=0.99, prior="BayesPrior")))
#' inc = InclusiveConnectedness(dca, group=c(1,2,3))
#' }
#' @references Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.
#' @author David Gabauer
#' @export
InclusiveConnectedness = function(dca, group=c(1,2), start=NULL, end=NULL) {
  message("The partial connectedness measures are implemented according to:\n Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.")
  corrected = dca$config$corrected
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Partial connectedness measures are not implemented for",dca$config$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = dca$CT[,,start:end]
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    k = dim(ct)[1]
    t = dim(ct)[3]
    
    CT = ct*0
    for (i in group) {
      CT[,i,] = ct[,i,]
      CT[i,,] = ct[i,,]
    }
  
    TCI = array(NA, c(t,1), dimnames=list(as.character(date), "TCI"))
    NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
    NPDC = PCI = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
    for (i in 1:t) {
      dca_ = ConnectednessTable(CT[,,i])
      NPDC[,,i] = dca_$NPDC
      PCI[,,i] = dca_$PCI
      infl = dca_$INFLUENCE
      infl[which(is.nan(infl), arr.ind=TRUE)] = 0
      INFLUENCE[,,i] = infl
      TO[i,] = dca_$TO
      FROM[i,] = dca_$FROM
      NET[i,] = dca_$NET
      NPT[i,] = dca_$NPT
      if (corrected) {
        TCI[i,] = dca_$cTCI
      } else {
        TCI[i,] = dca_$TCI
      }
    }
    TABLE = ConnectednessTable(CT)$TABLE
    config = list(approach="Inclusive")
    return = list(TABLE=TABLE, TCI=TCI, NET=NET, TO=TO, FROM=FROM, NPT=NPT,
                  NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
  }
}


#' @title Internal Connectedness Measures
#' @description This function provides internal connectedness measures
#' @param dca Dynamic connectedness object
#' @param groups List of at least two group vectors
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("gg2018")
#' dca = ConnectednessApproach(gg2018, 
#'                             nlag=1, 
#'                             nfore=10, 
#'                             window.size=200,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.99, 
#'                             prior="BayesPrior")))
#' int = InternalConnectedness(dca, groups=list("US"=c(1,2,3,4), "JP"=c(5,6,7,8)))
#' }
#' @references Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.
#' @author David Gabauer
#' @export
InternalConnectedness = function(dca, groups=list(c(1), c(2:ncol(dca$NET))), start=NULL, end=NULL) {
  corrected = dca$config$corrected
  message("The decomposed connectedness measures are implemented according to:\n Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.")
  if (length(groups)<=1) {
    stop("groups need to consist of at least 2 vectors")
  }
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Decomposed connectedness measures are not implemented for",dca$approach, "connectedness"))
  } else {
    if (dca$config$approach=="Extended Joint") {
      corrected = FALSE
    }
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = 100*dca$CT[,,start:end]
    NAMES = colnames(ct)
    k = dim(ct)[2]
    if (length(dim(ct))==2) {
      ct = array(ct, c(k,k,1),dimnames=list(NAMES,NAMES))
    }
    ct_wo = ct
    date = as.character(dimnames(ct)[[3]])
    t = dim(ct)[3]

    m = length(groups)
    NAMES_group = names(groups)
    if (is.null(NAMES_group)) {
      NAMES_group = paste0("GROUP", 1:m)
    }
  
    for (i in 1:m) {
      for (j in 1:m) {
        if (i>j) {
          group_1 = groups[[i]]
          group_2 = groups[[j]]
          ct_wo[group_1,group_2,] = 0
          ct_wo[group_2,group_1,] = 0
        }
      }
    }
  
    TCI_wo = array(NA, c(t, 1), dimnames=list(date, c("TCI")))
    PCI_wo = INFLUENCE_wo = NPDC_wo = array(NA, c(k, k, t), dimnames=list(NAMES,NAMES,date))
    TO_wo = FROM_wo = NET_wo = NPT_wo = array(NA, c(t, k), dimnames=list(date, NAMES))
    for (i in 1:t) {
      dca_ = ConnectednessTable(ct_wo[,,i]/100)
      TO_wo[i,] = dca_$TO
      FROM_wo[i,] = dca_$FROM
      NET_wo[i,] = dca_$NET
      NPT_wo[i,] = dca_$NPT
      NPDC_wo[,,i] = dca_$NPDC
      inf = dca_$INFLUENCE
      inf[which(is.nan(inf),arr.ind=TRUE)] = 0
      INFLUENCE_wo[,,i] = inf
      PCI_wo[,,i] = dca_$PCI
      if (corrected) {
        TCI_wo[i,] = dca_$cTCI
      } else {
        TCI_wo[i,] = dca_$TCI
      }
    }
    if (corrected) {
      m_ = (k-1)
    } else {
      m_ = k
    }
    TCI_group = array(NA, c(t,m), dimnames=list(date, NAMES_group))
    for (i in 1:m) {
      group = groups[i][[1]]
      TCI_group[,i] = rowSums(TO_wo[,group,drop=FALSE])/m_
    }
    TABLE = ConnectednessTable(ct_wo/100)$TABLE
    config = list(approach="Internal")
    return = list(TABLE=TABLE, gTCI=TCI_group, TCI=TCI_wo, TO=TO_wo, FROM=FROM_wo, 
                  NET=NET_wo, NPDC=NPDC_wo, PCI=PCI_wo, INFLUENCE=INFLUENCE_wo, NPT=NPT_wo, config=config)
  }
}



#' @title Lastrapes and Wiesen (2021) joint connectedness approach
#' @description This function calculates the Lastrapes and Wiesen (2021) joint connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data(lw2021)
#' fit = VAR(lw2021, configuration=list(nlag=2))
#' dca = JointConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=30)
#' dca$TABLE
#' }
#' @references
#' Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.
#' @author David Gabauer
#' @export
JointConnectedness = function(Phi, Sigma, nfore) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(Sigma)
  if (length(dim(Phi))==2) {
    Phi = array(Phi, c(nrow(Phi),ncol(Phi),1))
  }
  if (length(dim(Sigma))==2) {
    Sigma = array(Sigma, c(nrow(Sigma),ncol(Sigma),1))
  }

  k = ncol(Sigma)
  t = dim(Sigma)[3]

  if (is.null(NAMES)) {
    NAMES = 1:k
  }

  date = as.character(dimnames(Sigma)[[3]])
  TCI = array(NA, c(t,1), dimnames=list(date, "TCI"))
  NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
  CT = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
  CT_ = array(NA, c(k+3, k+1, t))
  pb = progress_bar$new(total=t)
  for (ij in 1:t) {
    # calculate the gFEVD
    gSOT = 100*FEVD(Phi[,,ij], Sigma[,,ij], nfore=nfore, type="time", generalized=TRUE)$FEVD
    gSOI = mean(rowSums(gSOT-diag(diag(gSOT))))

    # calculate Xi (the forecast error covariance matrix)
    A = Wold(Phi[,,ij], nfore)  # the VMA coefficient matrices
    Xi_h = array(0,dim=c(k,k,nfore))
    for (h in 1:nfore) {
      Xi_h[,,h] = A[,,h]%*%Sigma[,,ij]%*%t(A[,,h]) # calculated Xi at each h
    }
    Xi = rowSums(Xi_h, dims=2) # sum them along THIRD dimension to form Xi  (note: because this is a row sum, dims=2, actually sums along the third dimension)
    I_K = diag(1,nrow=k,ncol=k)

    #Calculate the elimination matrices.
    #Usually denoted as a KxK-1 matrix M_i, here it is an array where M[,,1]=M_1, and in general M[,,i]=M_i.
    M = array(0,dim=c(k,k-1,k))
    for (i in 1:k){
      M[,,i] = I_K[,-i]  #calculate the elimination matrices
    }

    #Calculate the joint total spillovers from all others to variable i (S_jnt)
    #calculate the numerator of S_jnt
    S_jnt_numerator_h = array(0,dim=c(k,nfore))
    for (i in 1:k){
      for (h in 1:nfore){
        S_jnt_numerator_h[i,h] = I_K[i,]%*%A[,,h]%*%Sigma[,,ij]%*%M[,,i]%*%(solve(t(M[,,i])%*%Sigma[,,ij]%*%M[,,i]))%*%t(M[,,i])%*%Sigma[,,ij]%*%t(A[,,h])%*%I_K[,i]     #calculate the numerator of S_jnt at each h
      }
    }
    S_jnt_from = array(0,dim=c(k))
    for (i in 1:k){
      S_jnt_from[i] = 100*sum(S_jnt_numerator_h[i,])/Xi[i,i]
    }

    jSOI = mean(S_jnt_from)
    lambda = gSOI/jSOI
    TCI[ij,] = jSOI
    FROM[ij,] = S_jnt_from
    TO[ij,] = colSums(gSOT-diag(diag(gSOT)))/lambda
    NET[ij,] = TO[ij,] - FROM[ij,]
    CT[,,ij] = gSOT
    CT_[,,ij] = rbind(rbind(rbind(cbind(gSOT, FROM[ij,]), c(TO[ij,], sum(TO[ij,]))), c(colSums(gSOT), 0)), c(NET[ij,], TCI[ij,]))
    pb$tick()
  }
  TABLE = apply(CT_,1:2,mean)
  rownames(TABLE) = c(NAMES, "TO", "Inc.Own", "NET")
  colnames(TABLE) = c(NAMES, "FROM")
  TABLE = format(round(TABLE, 2), nsmall=2)
  TABLE[k+2,k+1] = "TCI"
  
  config = list(nfore=nfore, approach="Joint", generalized=TRUE, corrected=FALSE)
  return = list(TABLE=TABLE, CT=CT/100, TCI=TCI, TO=TO, FROM=FROM, NET=NET,
                NPDC=NULL, NPT=NULL, PCI=NULL, INFLUENCE=NULL, config=config)
}


#' @title Minimum connectedness portfolio
#' @description This function calculates the minimum connectedness portfolio
#' @param x zoo return matrix (in percentage)
#' @param H Pairwise connectedness matrix or alternatively variance-covariance or correlation matrix
#' @param method Cumulative sum or cumulative product 
#' @param long Allow only long portfolio position
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @return Get portfolio weights
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' dca = TimeConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=10, generalized=TRUE)
#' mcp = MinimumConnectednessPortfolio(g2020, dca$PCI, statistics="Fisher")
#' mcp$TABLE
#' @references
#' Broadstock, D. C., Chatziantoniou, I., & Gabauer, D. (2020). Minimum Connectedness Portfolios and the Market for Green Bonds: Advocating Socially Responsible Investment (SRI) Activity. Available at SSRN 3793771.\\
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.\\
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
MinimumConnectednessPortfolio = function(x, H, method=c("cumsum","cumprod"), long=TRUE, statistics=c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), digit=2) {
  message("The minimum connectedness portfolio is implemented according to:\n Broadstock, D. C., Chatziantoniou, I., & Gabauer, D. (2020). Minimum Connectedness Portfolios and the Market for Green Bonds: Advocating Socially Responsible Investment (SRI) Activity. Available at SSRN 3793771.")
  message("Hedging effectiveness is calculated according to:\n Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.")
  message("Statistics of the hedging effectiveness measure are implemented according to:\n Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.")
  
  method = match.arg(method)
  statistics = match.arg(statistics)
  x = x / 100
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  t = nrow(x)
  date = as.character(rownames(x))
  NAMES = colnames(x)
  I = matrix(1, k, 1)
  if (dim(H)[[3]]==1) {
    H = array(H, c(k,k,t), dimnames=list(NAMES,NAMES,date))
  }

  portfolio_weights = array(NA,c(t,k), dimnames=list(date, NAMES))
  for (i in 1:t) {
    V_inv = MASS::ginv(H[,,i])
    pw = (V_inv%*%I) / c(t(I)%*%V_inv%*%I)
    if (long) {
      pw = ifelse(pw<0,0,pw)
      pw = ifelse(pw>1,1,pw)
      pw = pw/sum(pw)
    }
    portfolio_weights[i,] = pw
  }
  
  summary = NULL
  for (i in 1:k) {
    x_ = as.matrix(portfolio_weights[,i])
    summary_ = matrix(NA, nrow=ncol(x_), ncol=4)
    for (ij in 1:ncol(x_)){
      summary_[ij,] = matrix(c(mean(x_[,ij]), stats::sd(x_[,ij]), stats::quantile(x_[,ij],0.05), stats::quantile(x_[,ij],0.95)), nrow=1)
    }
    colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
    summary = rbind(summary, summary_)
  }
  rownames(summary) = NAMES
  
  portfolio_return = array(NA,c(t,1), dimnames=list(date))
  for (i in 1:t) {
    portfolio_return[i,] = sum(portfolio_weights[i,]*as.numeric(x[i,]))
  }
  
  if (method=="cumsum") {
    cumulative_portfolio_return = cumsum(portfolio_return)
  } else if (method=="cumprod") {
    cumulative_portfolio_return = cumprod(1+portfolio_return)-1
  }
  
  HE = pvalue = array(NA,c(k,1), dimnames=list(NAMES))
  for (i in 1:k) {
    HE[i,] = 1 - var(portfolio_return)/var(x[,i])
    df = rbind(data.frame(val=portfolio_return, group="A"), data.frame(val=x[,i], group="B"))
    if (statistics=="Fisher") {
      pvalue[i,] = var.test(x=portfolio_return, y=x[,i],ratio=1)$p.value
    } else if (statistics=="Bartlett") {
      pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Bartlett", verbose=F)$p.value
    } else if (statistics=="Fligner-Killeen") {
      pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Fligner", verbose=F)$p.value
    } else if (statistics=="Levene") {
      pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Levene", verbose=F)$p.value
    } else if (statistics=="Brown-Forsythe") {
      pvalue[i,] = onewaytests::bf.test(val~as.character(group), data=df, verbose=F)$p.value
    } else {
      stop("No valid hedging effectiveness statistics have been chosen.")
    }
  }
  
  TABLE = cbind(summary,HE,pvalue)
  colnames(TABLE)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
  
  return = list(TABLE=format(round(TABLE,digit),nsmall=digit), portfolio_weights=portfolio_weights, HE=HE, pvalue=pvalue, portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return)
}


#' @title Minnesota Prior
#' @description Get Minnesota Prior
#' @param gamma Diagonal value of variance-covariance matrix
#' @param k Number of series
#' @param nlag Lag length
#' @return Get Minnesota Prior
#' @examples
#' prior = MinnesotaPrior(0.1, k=4, nlag=1)
#' @references Koop, G., & Korobilis, D. (2010). Bayesian multivariate time series methods for empirical macroeconomics. Now Publishers Inc.
#' @author David Gabauer
#' @export
MinnesotaPrior = function(gamma=0.1, k, nlag) {
  if (k<=1) {
    stop("k represents the number of series and needs to be an integer larger than 1")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  m = nlag*(k^2)
  bprior = c(cbind(0*diag(k), matrix(0, ncol=(nlag-1)*k, nrow=k)))
  V_i = matrix(0, nrow=(m/k), ncol=k)
  for (i in 1:k) {
    for (j in 1:(m/k)) {
      V_i[j,i] = gamma/(ceiling(j/k)^2)
    }
  }
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  return = list(bprior=bprior, Vprior=Vprior)
}

#' @title Partial Contemporaneous Correlations
#'
#' @description Get partial contemporaneous correlations
#' @param Q variance-covariance matrix
#' @return Get partial contemporaneous correlations
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' pcc = PartialCorrelations(fit$Q)
#' @references Dahlhaus, R., & Eichler, M. (2003). Causality and graphical models in time series analysis. Oxford Statistical Science Series, 115-137.
#' @author David Gabauer
#' @export
PartialCorrelations = function (Q) {
  message("Partial correlations are computed according to:\n Dahlhaus, R., & Eichler, M. (2003). Causality and graphical models in time series analysis. Oxford Statistical Science Series, 115-137.")
  if (length(dim(Q))<=1) {
    stop("Q needs to be at least a 2-dimensional matrix")
  }
  k = dim(Q)[1]
  NAMES = colnames(Q)
  if (length(dim(Q))==2) {
    Q = array(Q, c(k,k,1), dimnames=list(NAMES,NAMES))
  }
  pcc = Q
  for (l in 1:dim(Q)[3]) {
    precision = MASS::ginv(Q[,,l])
    theta = diag(1/sqrt(diag(precision)))
    pcc[,,l] = -theta %*% precision %*% theta
  }
  return(pcc)
}


#' @title Quantile vector autoregression
#' @description Estimation of a QVAR using equation-by-equation quantile regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param tau quantile between 0 and 1
#' @param configuration model configuration
#' @return Estimate QVAR model
#' @examples
#' #data(dy2012)
#' #fit = QVAR(dy2012, configuration=list(nlag=1, tau=0.5))
#' @importFrom quantreg rq
#' @references
#' White, H., Kim, T. H., & Manganelli, S. (2015). VAR for VaR: Measuring tail dependence using multivariate regression quantiles. Journal of Econometrics, 187(1), 169-188.\\
#' Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters, 204, 109891.
#' @author David Gabauer
#' @export
QVAR = function(x, configuration=list(nlag=1, tau=0.5)) {
  tau = as.numeric(configuration$tau)
  nlag = as.numeric(configuration$nlag)
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (tau <= 0 || tau >= 1) {
    stop("tau needs to be within 0 and 1")
  }
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }

  Res = B = NULL
  for (i in 1:k) {
    z = embed(x, nlag+1)
    fit = rq(z[,i] ~ z[,-c(1:k)], tau=tau)
    B = rbind(B, fit$coefficients[-1])
    Res = cbind(Res, fit$residuals)
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(zoo::index(rownames(x)),1)))
  results = list(B=B, Q=Q)
}


#' @title Diebold and Yilmaz (2009, 2012) connectedness approach
#' @description This function allows to calculate the Diebold and Yilmaz (2009, 2012) connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param generalized Orthorgonalized/generalized FEVD
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param FEVD Alternatively, to provide Phi and Sigma it is also possible to use FEVD directly.
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of DY2012
#' data("dy2012")
#' fit = VAR(dy2012, configuration=list(nlag=4))
#' dca = TimeConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=10, generalized=TRUE)
#' dca$TABLE
#' }
#' @references
#' Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal, 119(534), 158-171.\\
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting, 28(1), 57-66.
#' @author David Gabauer
#' @export
TimeConnectedness = function(Phi=NULL, Sigma=NULL, nfore=10, generalized=TRUE, corrected=FALSE, FEVD=NULL) {
  if ((is.null(Phi) || is.null(Sigma)) && is.null(FEVD)) {
    stop("Either Sigma and Phi need to be given or FEVD")
  }
  if (is.null(FEVD)) {
    if (nfore<=0) {
      stop("nfore needs to be a positive integer")
    }
    if (length(dim(Sigma))<=1) {
      stop("Sigma needs to be at least a 2-dimensional matrix")
    }
    if (length(dim(Phi))<=1) {
      stop("Phi needs to be at least a 2-dimensional matrix")
    }
    NAMES = colnames(Sigma)
    k = ncol(Sigma)
    if (length(dim(Phi))==2) {
      Phi = array(Phi, c(nrow(Phi),ncol(Phi),1))
    }
    if (length(dim(Sigma))==2) {
      Sigma = array(Sigma, c(nrow(Sigma),ncol(Sigma),1))
    }
    t = dim(Sigma)[3]
  
    if (is.null(NAMES)) {
      NAMES = 1:k
    }
    date = as.character(as.Date(dimnames(Sigma)[[3]]))
    CT = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, as.character(date)))
    for (i in 1:t) {
      CT[,,i] = FEVD(Phi=Phi[,,i], Sigma=Sigma[,,i], nfore=nfore, generalized=generalized, type="time")$FEVD
    }
  } else {
    CT = FEVD
    t = dim(CT)[3]
    NAMES = dimnames(CT)[[1]]
    k = length(NAMES)
    date = as.character(as.Date(dimnames(CT)[[3]]))
  }
  TCI = array(NA, c(t,1), dimnames=list(date, "TCI"))
  NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
  PCI = NPDC = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
  pb = progress::progress_bar$new(total=t)
  for (i in 1:t) {
    ct = ConnectednessTable(CT[,,i])
    TO[i,] = ct$TO
    FROM[i,] = ct$FROM
    NET[i,] = ct$NET
    NPT[i,] = ct$NPT
    NPDC[,,i] = ct$NPDC
    INFLUENCE[,,i] = ct$INFLUENCE
    PCI[,,i] = ct$PCI
    if (corrected) {
      TCI[i,] = ct$cTCI
    } else {
      TCI[i,] = ct$TCI
    }
    pb$tick()
  }
  TABLE = ConnectednessTable(CT)$TABLE
  config = list(nfore=nfore, approach="Time", generalized=generalized, corrected=corrected)
  return = list(TABLE=TABLE, CT=CT, TCI=TCI, TO=TO, FROM=FROM,
                NET=NET, NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
}


#' @title Uninformative Prior
#' @description Get Uninformative Prior
#' @param k Number of series
#' @param nlag Lag length
#' @return Get Uninformative Prior
#' @examples
#' prior = UninformativePrior(k=4, nlag=1)
#' @references Koop, G., & Korobilis, D. (2010). Bayesian multivariate time series methods for empirical macroeconomics. Now Publishers Inc.
#' @author David Gabauer
#' @export
UninformativePrior = function(k, nlag) {
  if (k<=1) {
    stop("k represents the number of series and needs to be an integer larger than 1")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  m = nlag*(k^2)
  bprior = c(cbind(0*diag(k), matrix(0, ncol=(nlag-1)*k, nrow=k)))
  V_i = matrix(4, nrow=(m/k), ncol=k)
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  return = list(bprior=bprior, Vprior=Vprior)
}


#' @title Vector autoregression
#' @description Estimation of a VAR using equation-by-equation OLS regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param configuration model configuration
#' @return Estimate VAR model
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' @references Sims, C. A. (1980). Macroeconomics and reality. Econometrica, 1-48.
#' @author David Gabauer
#' @importFrom stats lm
#' @importFrom stats embed
#' @export
VAR = function(x, configuration=list(nlag=1)) {
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  nlag = as.numeric(configuration$nlag)
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  
  Res = B = se = NULL
  for (i in 1:k) {
 
    z = stats::embed(x, nlag+1)
    fit = summary(stats::lm(z[,i] ~ z[,-c(1:k)]))
	
    B = rbind(B, fit$coefficients[-1,1])
    se = rbind(se, fit$coefficients[-1,2])
    Res = cbind(Res, fit$residuals)
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(as.character(zoo::index(x)),1)))
  results = list(B=B, Q=Q, se=se)
}


#' @title Generalized volatility forecast error variance decomposition and volatility impulse response functions
#' @description This function provides the volatility impulse responses and the forecast error variance decomposition of DCC-GARCH models.
#' @param fit Fitted DCC-GARCH model
#' @param nfore H-step ahead forecast horizon
#' @param standardize Boolean value whether GIRF should be standardized
#' @return Get volatility impulse response functions and forecast error variance decomposition
#' @references
#' Engle, R. (2002). Dynamic conditional correlation: A simple class of multivariate generalized autoregressive conditional heteroskedasticity models. Journal of Business & Economic Statistics, 20(3), 339-350.\\
#' Gabauer, D. (2020). Volatility impulse response analysis for DCC‐GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.
#' @author David Gabauer
#' @importFrom rmgarch rcor
#' @importFrom rmgarch rcov
#' @export
VFEVD = function(fit, nfore=100, standardize=FALSE) {
  if (class(fit)!="DCCfit") {
    stop("fit needs to be of class DCCfit")
  }
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  NAMES = fit@model$modeldata$asset.names
  H = rmgarch::rcov(fit)
  R = rmgarch::rcor(fit)
  R.bar = apply(R,1:2,mean)
  Q.bar = fit@mfit$Qbar
  t = dim(H)[3]
  k = dim(H)[1]
  alpha = array(0,c(k,k,nfore))
  alpha[,,1] = diag(fit@mfit$matcoef[c(seq(3,(4*k),4)),1])
  # fit@mfit$matcoef contains the ("mu","omega","alpha1","beta1") of each
  # time series and the ("dcca1","dccb1") of the Joint
  beta = diag(fit@mfit$matcoef[c(seq(4,(4*k),4)),1])
  ALPHA = fit@mfit$matcoef[(4*k+1),1] #dcca1
  BETA = fit@mfit$matcoef[(4*k+2),1] #dccb1

  H.hat = array(0,c(k,k,nfore+1))
  VIRF = H.hat.shock = H.hat.no_shock = array(0,c(k,k,t,nfore+1))
  e = diag(k)
  for (i in 1:t) {
    H.hat[,,1] = H[,,i]
    Q.hat = H.hat
    Q.hat[,,1] = fit@mfit$Q[[i]]
    for (j in 1:nfore) {
      H.hat[,,j+1] = (alpha[,,j])%*%e^2 + beta%*%H.hat[,,j]
      D = diag(diag(H.hat[,,j+1])^0.5)
      u = D%*%e
      if (j==1) {
        Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + ALPHA*crossprod(u) + BETA*H.hat[,,1]
      } else {
        Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar + (ALPHA+BETA)*Q.hat[,,j]
      }
      R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
      H.hat[,,j+1] = D%*%R.hat%*%D
    }
    H.hat.shock[,,i,] = H.hat
  }
  if (standardize) {
    e = 0*diag(k)
    for (i in 1:t) {
      H.hat[,,1] = H[,,i]
      Q.hat = H.hat
      Q.hat[,,1] = fit@mfit$Q[[i]]
      for (j in 1:nfore) {
        H.hat[,,j+1] = beta%*%H.hat[,,j]
        D = diag(diag(H.hat[,,j+1])^0.5)
        if (j==1) {
          Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + BETA*H.hat[,,1]
        } else {
          Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar+(ALPHA+BETA)*Q.hat[,,j]
        }
        R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
        H.hat[,,j+1] = D%*%R.hat%*%D
      }
      H.hat.no_shock[,,i,] = H.hat
    }
  }
  for (i in 1:t) {
    VIRF[,,i,] = H.hat.shock[,,i,] - H.hat.no_shock[,,i,]
  }
  date = dimnames(H)[[3]]
  VFEVD = array(NA, c(k,k,t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:t) {
    num = apply(VIRF[,,i,]^2,1:2,sum)
    den = c(apply(num,1,sum))
    fevd = t(num)/den
    VFEVD[,,i] = (fevd/apply(fevd, 1, sum))
  }
  return = list(IRF=VIRF, FEVD=VFEVD)
}

#' @title Wold representation theorem
#' @description Transform VAR to VMA coefficients
#' @param x VAR coefficients
#' @param nfore H-step ahead forecast horizon
#' @return Get VMA coefficients
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' wold = Wold(fit$B, nfore=10)
#' @references Sims, C. A. (1980). Macroeconomics and reality. Econometrica, 1-48.
#' @author David Gabauer
#' @export
Wold = function (x, nfore=10) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  nstep = abs(as.integer(nfore))
  K = nrow(x)
  p = floor(ncol(x)/K)
  A = array(0, c(K,K,nstep))
  for (i in 1:p){
    A[,,i]=x[,((i-1)*K+1):(i*K)]
  }
  Phi = array(0, dim = c(K, K, nstep + 1))
  Phi[,,1]=diag(K)
  Phi[,,2]=Phi[,,1] %*% A[,,1]
  if (nstep > 1) {
    for (i in 3:(nstep + 1)) {
      tmp1 = Phi[,,1] %*% A[,,i-1]
      tmp2 = matrix(0, nrow=K, ncol=K)
      idx = (i-2):1
      for (j in 1:(i-2)) {
        tmp2 = tmp2 + Phi[,,j+1] %*% A[,,idx[j]]
      }
      Phi[,,i] = tmp1 + tmp2
    }
  }
  return(Phi)
}



