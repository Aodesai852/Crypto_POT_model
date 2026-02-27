rm(list=ls())
library(writexl)
library(ggplot2)
library(GenSA)
library(foreach)
library(doParallel)
library(data.table) 
library(xts) 
library(numDeriv)
library(xts)
library(fBasics)
library(tseries)
library(TSA)
library(rugarch)
library(rcompanion)
library(readxl)
library(scales)
library(lubridate)
library(openxlsx)
library(moments)
library(svglite)
library(dplyr)
library(zoo)
library(urca)
library(stringr)

setwd("D:/MyFiles/EVT/Crypto_POT")
source("./code/POT/Util_UnivPoT.R")
source("./code/POT/runtime_func.R")
source("./code/POT/events.R")
source("./code/Params/est_param_down.R")
source("./code/Params/est_param_up.R")
load("./data/panel_dt.RData")

set.seed(123)

q = 0.95 # main:0.95; rob:0.90/0.85
step <- 0.001

# Setup parallel backend
cl <- makeCluster(detectCores()/2)
registerDoParallel(cl)

# hyper-parameters
num_of_para <- 6
scale_factor <- 100
Sigma_type <- "Square_Sigma"; # "Log_Sigma"
ktype <- 'fixed'              # c('fxed', 'nonfixed') 
burnin_est <- 0               # estimiation burnin
init_upper_xi_sigma <- c(10, 5) # (xi,sigma)
init_lower_xi_sigma <- c(1, 0.001)

# (phi_0 (constant term),phi_1(beta term),phi_2(alpha term),psi_0(constant term),psi_1(beta term),psi_2(alpha term))
# phi_i is parameter for xi, psi_i is parameter for sigma
upperbound_para_default <- c(0.3, 0.999,     0.3,     0.3, 0.999,     0.3) # general para
lowerbound_para_default <- c(0,     0.4, 0.00001, 0.00001,   0.4, 0.00001) # general para
#upperbound_para_default <- c(0.1,  0.999,  0.5, 0.5,   0.99, 0.01) # another general para
#lowerbound_para_default <- c(-0.2, 0.8,0.001, 0.001, 0.7, 0.001) # another general para
#upperbound_para_default <- c(0.3,  0.8,  0.5,   0.3,   0.9, 0.001) # another general para
#lowerbound_para_default <- c(0,    0.6,    0.001, 0.001, 0.7, 0.00001) # another general para

panel.dt = na.omit(panel.dt)

# choose which coins you want to estimate
coin.list.est = colnames(panel.dt)[-1]
# coin.list.est = c("AXS","BAT","BDX","BSV","CAKE","CFX","COMP","CRV","DASH","DCR",
#                   "DEXE","FET","GNO","HNT","INJ","IOTA","JST","LUNC","MANA","MX",
#                   "NEO","NEXO","SAND","STX","TEL","THETA","TRAC","TWT","XTZ")

curr.coin = "AAVE" # for test
for (curr.coin in coin.list.est) {

  # start with defaults
  upperbound_para <- upperbound_para_default
  lowerbound_para <- lowerbound_para_default
  
  # override if coin has selected params
  if (curr.coin %in% names(est.param.up)) {
    upperbound_para_up <- est.param.up[[curr.coin]] * 1.1
    lowerbound_para_up <- est.param.up[[curr.coin]] * 0.9
  }
  
  if (curr.coin %in% names(est.param.down)) {
    upperbound_para_down <- est.param.down[[curr.coin]] * 1.1
    lowerbound_para_down <- est.param.down[[curr.coin]] * 0.9
  }
  
  curr.ts = panel.dt[[curr.coin]] # current time series
  names(curr.ts) <- panel.dt[["Date"]]
  ret = diff(log(curr.ts))*scale_factor # log return * 100
  date = as.Date(names(ret))
  
  ret.desc.stats = data_outline(ret) # to be saved
  ret_plot_handle = return_plot(ret,curr.coin) # to be saved
  price_plot_handle = return_plot(curr.ts,curr.coin,price_plot = TRUE) # to be saved
  
  raw_obs_neg = -unname(ret)
  raw_obs_pos = unname(ret)
  
  # calc hurdle (tau)
  tau_global_neg = as.numeric(quantile(raw_obs_neg,q)) # downside risk hurdle
  tau_global_pos = as.numeric(quantile(raw_obs_pos,q)) # upside risk hurdle
  
  # downside and upside Y = max(R_t - tau,0)
  obs_neg = raw_obs_neg - tau_global_neg
  obs_neg[obs_neg < 0] = 0
  obs_pos = raw_obs_pos - tau_global_pos
  obs_pos[obs_pos < 0] = 0
  
  # save the descriptive results
  save(ret.desc.stats, file = paste0("./result/ret_desc_stats/",curr.coin,".RData"))
  ggsave(paste0("./result/return_plot/", curr.coin, ".pdf"),       
         plot = ret_plot_handle,                     
         width = 8, height = 4)        # inches
  
  ggsave(paste0("./result/price_plot/", curr.coin, ".pdf"),       
         plot = price_plot_handle,                     
         width = 8, height = 4)        # inches
  
  try.result = try({
    
  # estimating lower tail
  neg_optim_result <- optim_result(raw_obs_neg,obs_neg,num_of_para,init_lower_xi_sigma,init_upper_xi_sigma,
                           lowerbound_para_down,upperbound_para_down,burnin_est,ktype,Sigma_type,tau_global_neg,
                           date,step,curr.coin,NE=TRUE)
  
  neg_POT_plot <- POT_plot(xiSigma_timeseries_df = neg_optim_result$xiSigma_timeseries_df, range_list = events[[curr.coin]],
                           curr.coin = curr.coin, NE = TRUE)
  
  save(neg_optim_result, file = paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/optim_result/",curr.coin,".RData"))
  
  ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/P_xi_exc/pdf/",curr.coin,".pdf"),       
         plot = neg_POT_plot$P_xi_exc,                     
         width = 8, height = 4)        # inches
  
  ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/P_std_sigma/pdf/",curr.coin,".pdf"),       
         plot = neg_POT_plot$P_std_sigma,                     
         width = 8, height = 4)        # inches
  
  }, silent = TRUE) # end of try
  
  if (inherits(try.result, "try-error")) {
    message(paste("Error in lowertail", curr.coin, "— skipping to next"))
    next
  }
  
  try.result = try({
    
  # estimating upper tail
  pos_optim_result <- optim_result(raw_obs_pos,obs_pos,num_of_para,init_lower_xi_sigma,init_upper_xi_sigma,
                                   lowerbound_para_up,upperbound_para_up,burnin_est,ktype,Sigma_type,tau_global_pos,
                                   date,step,curr.coin,NE=FALSE)
  
  pos_POT_plot <- POT_plot(xiSigma_timeseries_df = pos_optim_result$xiSigma_timeseries_df, range_list = events[[curr.coin]],
                           curr.coin = curr.coin, NE=FALSE)
  
  save(pos_optim_result, file = paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/optim_result/",curr.coin,".RData"))
  ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/P_xi_exc/pdf/",curr.coin,".pdf"),       
         plot = pos_POT_plot$P_xi_exc,                     
         width = 8, height = 4)        # inches
  
  ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/P_std_sigma/pdf/",curr.coin,".pdf"),       
         plot = pos_POT_plot$P_std_sigma,                     
         width = 8, height = 4)        # inches
  
  }, silent = TRUE) # end of try
  
  if (inherits(try.result, "try-error")) {
    message(paste("Error in uppertail", curr.coin, "— skipping to next"))
    next
  }
  
}

# Stop cluster
if (!plot.mode){
  stopCluster(cl)}