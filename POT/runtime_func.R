# ---- Runtime Functions ----

# Function to standardize one column
standardize_col <- function(x, name) {
  if (grepl("^USD/", name) || name == "DXY") {
    # Already in USD/XXX form (or special case DXY) -> return unchanged
    list(values = x, name = name)
  } else if (grepl("/USD$", name)) {
    # In form XXX/USD -> flip to USD/XXX
    parts <- unlist(strsplit(name, "/"))
    new_name <- paste0("USD/", parts[1])
    list(values = 1/x, name = new_name)
  } else {
    # Unexpected format: leave unchanged but warn
    warning("Unrecognized format: ", name)
    list(values = x, name = name)
  }
}

# descriptive statistics
data_outline <- function(x){
  # x is a zoo object or a numeric vector
  # JB test
  JB <- tseries::jarque.bera.test(x)
  if(JB$p.value<=0.01){JB_conclusion = "1%，non-normal"
  }else if(JB$p.value<=0.05){JB_conclusion = "5%，non-normal"
  }else if(JB$p.value<=0.1){JB_conclusion = "10%，non-normal"
  }else{JB_conclusion = "normal"}
  
  # Augmented Dickey-Fuller test
  adf <- tseries::adf.test(x)
  if(adf$p.value<=0.01){adf_conclusion = "1%，stationary"
  }else if(adf$p.value<=0.05){adf_conclusion = "5%，stationary"
  }else if(adf$p.value<=0.1){adf_conclusion = "10%，stationary"
  }else{adf_conclusion = "non-stationary"}
  
  # PP test（Phillips-Perron test）
  pp <- tseries::pp.test(x)
  if(pp$p.value<=0.01){pp_conclusion = "1%，stationary"
  }else if(pp$p.value<=0.05){pp_conclusion = "5%，stationary"
  }else if(pp$p.value<=0.1){pp_conclusion = "10%，stationary"
  }else{pp_conclusion = "non-stationary"}
  
  # ZA test（Zivot-Andrews test)
  za <- urca::ur.za(x)
  za_critical<-c(za@cval[1], za@cval[2], za@cval[3])
  if(za@teststat<=za@cval[1]){
    za_conclusion = "1%，stationary"
  }else if(za@teststat<=za@cval[2]){
    za_conclusion = "5%，stationary"
  }else if(za@teststat<=za@cval[3]){
    za_conclusion = "10%，stationary"
  }else{
    za_conclusion = "non-stationary"
  }
  
  return(list(mean = mean(x), median = median(x), min = min(x), max = max(x), 
              std_dev = sd(x), skew = skewness(x), kurt = kurtosis(x), 
              JB_stat = as.numeric(JB$statistic), JB_pvalue = JB$p.value, JB_conclusion = JB_conclusion,
              adf_stat = adf$statistic, adf_pvalue = adf$p.value, adf_conclusion = adf_conclusion,
              pp_stat = pp$statistic, pp_pvalue = pp$p.value, pp_conclusion = pp_conclusion,
              za_stat = za@teststat, za_critical = za_critical, za_conclusion = za_conclusion, ts = x))
}

return_plot <- function(ret, curr.coin, price_plot = FALSE){

  dt = data.table(date = as.Date(names(ret)), value = unname(ret))
  
  if (price_plot){
    title = paste("Spot",curr.coin)
  }else{
    title = paste("Log return of spot", curr.coin)}
  
  P = ggplot(dt, aes(x = date,y = value)) +
    geom_line(color="black",linewidth=0.5) + 
    scale_x_date(
      date_breaks = "2 years",          # one scale per year
      date_labels = "%Y",              # only display year
      expand = c(0, 0)                 # no expansion
    ) +
    labs(
      title = title,
      x = "",
      y = ""
    ) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),    # remove all grids
      panel.border = element_rect(color = "black"), # keep borders
      plot.title = element_text(hjust = 0.5)
      #axis.text.x = element_text(angle = 45, hjust = 1) # Year label tilted 45 degrees to prevent overlapping
    )
  return(P)
} # end of function

# optimization result
optim_result <- function(raw_obs,obs,num_of_para,init_lower_xi_sigma,init_upper_xi_sigma,
                         lowerbound_para,upperbound_para,burnin_est,ktype,Sigma_type,tau_global,
                         date,step,curr.coin,NE){
  # Stage 0: Initial estimate of static GPD model
  temp_result <- c(); reptimes <- 100
  for(rep in 1: reptimes){
    # A random start point
    #print(paste("static",rep))
    start <- init_lower_xi_sigma + (init_upper_xi_sigma-init_lower_xi_sigma)*runif(length(init_lower_xi_sigma))
    # MLE
    temp <- nlminb(start, initEstUniGPDMarg, upper = init_upper_xi_sigma, lower = init_lower_xi_sigma,
                   y = obs, tau = tau_global)
    temp_result <- rbind(temp_result, c(temp$par, temp$objective)) 
  }
  init_xi_sigma <- temp_result[which.min(temp_result[, (length(init_lower_xi_sigma)+1)]), c(1:length(init_lower_xi_sigma))]
  
  ########################
  # Stage 1: Estimate of Dynamic GPD-POT
  loc_seq <- seq(0, 1, step)
  reptimes <- length(loc_seq)
  
  print(paste("Optimizing",curr.coin))
  # Parallel execution with foreach
  temp_result <- foreach(rep = 1:reptimes, .combine = rbind,
                         .export = c("curr.coin", "lowerbound_para", "upperbound_para",
                                     "tau_global", "obs", "raw_obs", "init_xi_sigma",
                                     "burnin_est", "ktype", "Sigma_type", "dGPD",
                                     "loglikUnivariate_univPOT_garch", "recoverXiSigma_univPOT_garch")) %dopar% {
                                       
                                       print(paste(curr.coin, "dynamic est.", rep))
                                       loc <- runif(1, 0, 1)
                                       start <- loc * lowerbound_para + (1 - loc) * upperbound_para
                                       
                                       temp <- nlminb(start, loglikUnivariate_univPOT_garch,
                                                      upper = upperbound_para, lower = lowerbound_para,
                                                      tau = tau_global, y_obs = obs, raw_obs = raw_obs,
                                                      init = init_xi_sigma, burnin = burnin_est,
                                                      ktype = ktype, Sigma_type = Sigma_type)
                                       
                                       c(temp$par, temp$objective)
                                     }
  
  temp_result<-temp_result[order(temp_result[,num_of_para + 1]),] # sort the results, ascending
  temp_result <- temp_result[!apply(temp_result, 1, function(x) any(is.infinite(x))),] # romove rows with Inf
  
  for(m in 1:nrow(temp_result)){
    temp_min <- temp_result[m,]
    est_para <- list(c(tau_global, 1), temp_min[1:((num_of_para)/2)],
                     temp_min[((num_of_para)/2+1):(num_of_para)])
    neg_loglik <- temp_min[num_of_para+1]
    grad <- grad(loglikUnivariate_univPOT_garch, x = c(est_para[[2]], est_para[[3]]), method = 'Richardson', method.args=list(d=1e-4),
                 tau = tau_global, y_obs = obs, raw_obs = raw_obs,
                 init = init_xi_sigma, burnin = burnin_est, ktype = ktype, Sigma_type = Sigma_type)
    hess <- hessian(loglikUnivariate_univPOT_garch, x = c(est_para[[2]], est_para[[3]]), method = 'Richardson', method.args=list(d=1e-4),
                    tau = tau_global, y_obs = obs, raw_obs = raw_obs,
                    init = init_xi_sigma, burnin = burnin_est, ktype = ktype, Sigma_type = Sigma_type)
    stderr <- sqrt(diag(solve(hess))) # solve(hess) returns the inverse of hess 
    if (any(is.na(stderr))) {
      print(paste0("pass ",m))
      next
    }
    # calc significance，Wald test
    wald_stat = (c(est_para[[2]],est_para[[3]])/stderr[1:6])^2
    df <- 1  # degree of freedom
    # calc p value
    p_values<- 1 - pchisq(wald_stat, df)
    break
  }
  # get the filtered time series of xi and sigma
  xiSigma_timeseries_df = data.frame(date=date, recoverXiSigma_univPOT_garch(para = est_para, raw_obs = raw_obs, init = init_xi_sigma, Sigma_type = Sigma_type), exceedance = obs*ifelse(NE,-1,1))
  # xiSigma_timeseries_df has four columns, date, xi, sigma, and PE/NE
  
  # GARCH
  sGARCHspec<-ugarchspec(mean.model = list(armaOrder=c(0,0),include.mean=FALSE),
                         variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                         distribution.model = "norm")
  sGARCHfit <- ugarchfit(sGARCHspec,data=raw_obs)
  if(any(is.null(sGARCHfit@fit$sigma))){print("Null in garch model!")}
  
  # Standardization
  std_sigma_pot=(xiSigma_timeseries_df$sigma-mean(xiSigma_timeseries_df$sigma))/sd(xiSigma_timeseries_df$sigma)
  std_sigma_garch <- (sGARCHfit@fit$sigma-mean(sGARCHfit@fit$sigma))/sd(sGARCHfit@fit$sigma)
  xiSigma_timeseries_df$std_sigma_pot = std_sigma_pot # standardized sigma (scale para) filtered by the POT Zhao(2021) model
  xiSigma_timeseries_df$std_sigma_garch = std_sigma_garch # standardized sigma filtered by the garch model
  
  # correlation
  cor_pearson <- cor(std_sigma_pot, std_sigma_garch, method = "pearson")
  cor_diff_pearson <- cor(diff(std_sigma_pot), diff(std_sigma_garch), method = "pearson")
  cor_kendall <- cor(std_sigma_pot, std_sigma_garch, method = "kendall")
  cor_diff_kendall <- cor(diff(std_sigma_pot), diff(std_sigma_garch), method = "kendall")
  cor_spearman <- cor(std_sigma_pot, std_sigma_garch, method = "spearman")
  cor_diff_spearman <- cor(diff(std_sigma_pot), diff(std_sigma_garch), method = "spearman")
  
  corr <- c(cor_pearson,cor_diff_pearson,cor_kendall,cor_diff_kendall,cor_spearman,cor_diff_spearman)
  
  return(list(est_para=est_para,neg_loglik=neg_loglik,stderr=stderr,wald_stat=wald_stat,p_values=p_values,xiSigma_timeseries_df=xiSigma_timeseries_df,corr=corr))
  
} # end of function

compress_trans <- function() {
  # Divide values below 0 by 50 to compress the scale
  trans <- function(x) { ifelse(x < 0, x / 5, x) }
  # Inverse transform
  inv <- function(x) { ifelse(x < 0, x *5, x) }
  trans_new("compress", trans, inv, breaks = pretty_breaks())
}


get_crisis_date <- function(range,xiSigma_timeseries_df){
  # range: something like c("2005-01-01","2005-9-30")
  # the crisis date: the date with the lowest tail index value during a specific time interval
  # the plot period: (crisis date - 2 months, crisis date + 2 months)
  start_date <- range[1]
  end_date <- range[2]
  df_selected <- xiSigma_timeseries_df %>% filter(date >= start_date & date <= end_date)
  min_tail_index <- min(df_selected$xi) # the lower the tail index, the higher the tail risk
  crisis_date<-df_selected$date[which.min(df_selected$xi)]
  return(list(crisis_date=crisis_date,min_tail_index=min_tail_index))
}

get_shadow_area <- function(range_list,xiSigma_timeseries_df){
  shadow_width <- 15 # one-sided shadow_width in days
  #event_label_position_y = ifelse(NE,0.4,1.2)
  xmi <- c() # xmin: start date array (the crisis date - 2 months)
  xmx <- c() # xmax: end date array (the crisis date + 2 months)
  #lab <- c() # c("event1","event2",...)
  xx <- c() # the crisis date
  #yy <- c() # c(0.2,0.2,0.2,...)
  yyy<-c() # tail index values of the corresponding crisis dates
  ymin <- -Inf
  ymax <- Inf
  for (ii in 1:length(range_list)) {
    range = range_list[[ii]]
    xx_temp = as.Date(get_crisis_date(range,xiSigma_timeseries_df)$crisis_date)
    yyy_temp = get_crisis_date(range,xiSigma_timeseries_df)$min_tail_index
    xmi_temp = xx_temp - days(shadow_width)
    xmx_temp = xx_temp + days(shadow_width)
    xmi <- as.Date(c(xmi, xmi_temp))
    xmx <- as.Date(c(xmx, xmx_temp))
    xx <- as.Date(c(xx, xx_temp))
    #yy <- c(yy, event_label_position_y)
    yyy<-c(yyy,yyy_temp)
    #lab <- c(lab, paste("event", ii))
  }
  shadow_area <- data.frame(xx = xx, yyy = yyy, xmi = xmi, xmx = xmx, ymi = rep(ymin, length(xmi)), ymx = rep(ymax, length(xmi)))
  return(shadow_area)
}


POT_plot <- function(xiSigma_timeseries_df,range_list,curr.coin,NE){
  # range_list = NULL
  if (NE){
    title = paste("Tail Index and Negative Exceedances of",curr.coin)} else {
    title = paste("Tail Index and Positive Exceedances of",curr.coin)
  }
  scaling = 10 # exc = exc / scaling
  # shadow_area = get_shadow_area(range_list,xiSigma_timeseries_df)
  # xx=shadow_area$xx;yyy=shadow_area$yyy;xmi=shadow_area$xmi;ymi=shadow_area$ymi;ymx=shadow_area$ymx
  P_xi_exc <- ggplot(xiSigma_timeseries_df, aes(x = date)) +
    geom_line(aes(y = xi, color = "Tail index (top)"), size = 0.6, linetype = "solid", na.rm = TRUE) +
    geom_line(aes(y = ifelse(NE,-1,1)*exceedance/scaling, color = "Exceedances (bottom)"), size = 0.5, linetype = "solid", na.rm = TRUE) +
    scale_color_manual(values = c("Tail index (top)" = "red", "Exceedances (bottom)" = 'black')) +
    #geom_rect(data = shadow_area, aes(xmin = xmi, xmax = xmx, ymin = ymi, ymax = ymx), fill ="#333", alpha = 0.2, inherit.aes = FALSE)+
    #annotate("text", x = xx, y = yy, label = lab, angle = 90, vjust = 0.3, hjust = 0, size = 4)+
    theme(
      legend.position = "bottom",            
      legend.justification = "center",        
      legend.direction = "horizontal",        
      legend.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA, color = NA),
      legend.title=element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1) ,
      panel.spacing = unit(2, "cm"),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )+
    ggtitle(title)+
    scale_x_date(expand = c(0, 0),date_breaks = "2 year",date_labels = "%Y")
  
  # std_sigma_pot and std_sigma_garch
  P_std_sigma <- ggplot(xiSigma_timeseries_df, aes(x = date)) +
    geom_line(aes(y = std_sigma_pot, color = "POT")) +
    geom_line(aes(y = std_sigma_garch, color = "GARCH")) +
    scale_color_manual(values = c("POT" = "blue","GARCH"="red")) +
    theme(
      legend.position = "bottom",            
      legend.justification = "center",        
      legend.direction = "horizontal",   
      legend.background = element_blank(),
      legend.title=element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1) ,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )+
    scale_x_date(expand = c(0, 0),date_breaks = "2 year",date_labels = "%Y")+
    ggtitle(paste("Scale Parameters (Standardized) of", curr.coin))
  
  return(list(P_xi_exc=P_xi_exc,P_std_sigma=P_std_sigma))
} # end of function