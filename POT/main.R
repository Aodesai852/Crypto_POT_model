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

plot.mode = FALSE
rob.mode = FALSE
rob.q = 0.85

# Setup parallel backend
if (!plot.mode){
  cl <- makeCluster(detectCores()/2)
  registerDoParallel(cl)
  }

setwd("D:/MyFiles/EVT/Crypto_POT")
source("./code/POT/Util_UnivPoT.R")
source("./code/POT/runtime_func.R")
source("./code/POT/events.R")
load("./data/panel_dt.RData")
load("./data/panel_dt_Feb25added.RData")
set.seed(123)

# hyper-parameters
num_of_para <- 6
scale_factor <- 100
if (rob.mode){
  q <- rob.q
} else {
  q <- 0.95
}
Sigma_type <- "Square_Sigma"; # "Log_Sigma"
ktype <- 'fixed'              # c('fxed', 'nonfixed') 
burnin_est <- 0               # estimiation burnin
step <- 0.001
init_upper_xi_sigma <- c(10, 5) # (xi,sigma)
init_lower_xi_sigma <- c(1, 0.001)

# (phi_0,phi_1,phi_2,psi_0,psi_1,psi_2)
# phi_i is parameter for xi, psi_i is parameter for sigma
upperbound_para <- c(0.3,  0.999,  0.3, 0.3,   0.999, 0.3) # general para
lowerbound_para <- c(-0.9, 0.4,0.00001, 0.00001, 0.4, 0.00001) # general para
upperbound_para <- c(0.00135,0.996,0.015,0.3,0.924,0.0085)
lowerbound_para <- upperbound_para*0.95
#upperbound_para <- c(0.1,  0.999,  0.5, 0.5,   0.99, 0.01) # another general para
#lowerbound_para <- c(-0.2, 0.8,0.001, 0.001, 0.7, 0.001) # another general para
#upperbound_para <- c(0.5, 0.8, 0.05, 0.5,   0.8, 0.05)
#lowerbound_para <- c(0, 0.6,0.02, 0.001, 0.6, 0.02)
# upperbound_para <- c(0.3,  0.8,  0.5,   0.3,   0.9, 0.001) # another general para
# lowerbound_para <- c(0,    0.6,    0.001, 0.001, 0.7, 0.00001) # another general para


# process panel.dt
panel.dt = na.omit(panel.dt)
name.map = 
c(
  "aave-token"       = "AAVE",
  "algorand"         = "ALGO",
  "avalanche"        = "AVAX",
  "binance-coin"     = "BNB",
  "bitcoin-cash"     = "BCH",
  "bitcoin"          = "BTC",
  "bitget-token-new" = "BGB",
  "cardano"          = "ADA",
  "chainlink"        = "LINK",
  "cronos"           = "CRO",
  "dogecoin"         = "DOGE",
  "ethereum-classic" = "ETC",
  "ethereum"         = "ETH",
  "filecoin"         = "FIL",
  "gatechaintoken"   = "GT",
  "hedera-hashgraph" = "HBAR",
  "kucoin"           = "KCS",
  "litecoin"         = "LTC",
  "monero"           = "XMR",
  "near-protocol"    = "NEAR",
  "okb"              = "OKB",
  "polkadot"         = "DOT",
  "quant-network"    = "QNT",
  "render-token"     = "RNDR",
  "ripple"           = "XRP",
  "shiba-inu"        = "SHIB",
  "solana"           = "SOL",
  "stellar"          = "XLM",
  "tron"             = "TRX",
  "uniswap"          = "UNI",
  "unus-sed-leo"     = "LEO",
  "vechain-token"    = "VET",
  "xdc-network"      = "XDC",
  "zcash"            = "ZEC"
)
setnames(
  panel.dt,
  old = intersect(names(panel.dt), names(name.map)),
  new = name.map[intersect(names(panel.dt), names(name.map))]
)
#common_cols <- intersect(colnames(panel.dt), colnames(panel.dt.Feb25added))
panel.dt = panel.dt[panel.dt.Feb25added, on = "Date"]
panel.dt = na.omit(panel.dt)

# choose which coins you want to estimate
coin.list.est = c("AXS","BAT","BDX","BSV","CAKE","CFX","COMP","CRV","DASH","DCR",
                  "DEXE","FET","GNO","HNT","INJ","IOTA","JST","LUNC","MANA","MX",
                  "NEO","NEXO","SAND","STX","TEL","THETA","TRAC","TWT","XTZ")

curr.coin = "AAVE" # for test
for (curr.coin in coin.list.est) {
  try.result = try({
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
  
  if (!plot.mode){
    neg_optim_result <- optim_result(raw_obs_neg,obs_neg,num_of_para,init_lower_xi_sigma,init_upper_xi_sigma,
                             lowerbound_para,upperbound_para,burnin_est,ktype,Sigma_type,tau_global_neg,
                             date,step,curr.coin,NE=TRUE)
    pos_optim_result <- optim_result(raw_obs_pos,obs_pos,num_of_para,init_lower_xi_sigma,init_upper_xi_sigma,
                                     lowerbound_para,upperbound_para,burnin_est,ktype,Sigma_type,tau_global_pos,
                                     date,step,curr.coin,NE=FALSE)
  }else{
    if (rob.mode){
      load(paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/optim_result/",gsub("/", "_", curr.coin),".RData")) # neg_optim_result
      load(paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/optim_result/",gsub("/", "_", curr.coin),".RData")) # pos_optim_result
    } else {
      load(paste0("./result/down_tail/optim_result/",gsub("/", "_", curr.coin),".RData")) # neg_optim_result
      load(paste0("./result/up_tail/optim_result/",gsub("/", "_", curr.coin),".RData")) # pos_optim_result
    }
    }
  
  if (plot.mode){
    print(paste("Ploting",curr.coin))
  }
  
  # reciprocal tail index
  #reci_neg_tail_index = 1 / neg_optim_result$xiSigma_timeseries_df$xi
  #reci_pos_tail_index = 1 / pos_optim_result$xiSigma_timeseries_df$xi
  
  #up.tail.desc.stats = data_outline(reci_pos_tail_index)
  #down.tail.desc.stats = data_outline(reci_neg_tail_index)
  
  #desc.stats = list(return = ret.desc.stats, up_tail_index = up.tail.desc.stats, down_tail_index = down.tail.desc.stats)
  
  neg_POT_plot <- POT_plot(xiSigma_timeseries_df = neg_optim_result$xiSigma_timeseries_df, range_list = events[[curr.coin]],
                           curr.coin = curr.coin, NE = TRUE)
  pos_POT_plot <- POT_plot(xiSigma_timeseries_df = pos_optim_result$xiSigma_timeseries_df, range_list = events[[curr.coin]],
                           curr.coin = curr.coin, NE=FALSE)
  
  # save the results
  if (!rob.mode){
    save(ret.desc.stats, file = paste0("./result/ret_desc_stats/",curr.coin,".RData"))
  }
  if (!plot.mode){
    if (rob.mode){
      save(neg_optim_result, file = paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/optim_result/",curr.coin,".RData"))
      save(pos_optim_result, file = paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/optim_result/",curr.coin,".RData"))
    } else{
      save(neg_optim_result, file = paste0("./result/down_tail/optim_result/",curr.coin,".RData"))
      save(pos_optim_result, file = paste0("./result/up_tail/optim_result/",curr.coin,".RData"))
    }
  }
  
  if (!rob.mode){
  ggsave(paste0("./result/return_plot/", curr.coin, ".pdf"),       
         plot = ret_plot_handle,                     
         width = 8, height = 4)        # inches
  
  ggsave(paste0("./result/price_plot/", curr.coin, ".pdf"),       
         plot = price_plot_handle,                     
         width = 8, height = 4)        # inches
  #save(ret_plot_handle, file = paste0("./result/plot_handle/return/", curr.coin, ".RData"))
  #save(price_plot_handle, file = paste0("./result/plot_handle/price/", curr.coin, ".RData"))
  }
  
  if (rob.mode){
    ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/P_xi_exc/pdf/",curr.coin,".pdf"),       
           plot = pos_POT_plot$P_xi_exc,                     
           width = 8, height = 4)        # inches
    
    ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/up_tail/P_std_sigma/pdf/",curr.coin,".pdf"),       
           plot = pos_POT_plot$P_std_sigma,                     
           width = 8, height = 4)        # inches
    
    ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/P_xi_exc/pdf/",curr.coin,".pdf"),       
           plot = neg_POT_plot$P_xi_exc,                     
           width = 8, height = 4)        # inches
    
    ggsave(paste0("./result/rob_test/",sprintf("%.2f", q),"/down_tail/P_std_sigma/pdf/",curr.coin,".pdf"),       
           plot = neg_POT_plot$P_std_sigma,                     
           width = 8, height = 4)        # inches
  } else {
    ggsave(paste0("./result/up_tail/P_xi_exc/",curr.coin,".pdf"),       
           plot = pos_POT_plot$P_xi_exc,                     
           width = 8, height = 4)        # inches
    
    ggsave(paste0("./result/up_tail/P_std_sigma/",curr.coin,".pdf"),       
           plot = pos_POT_plot$P_std_sigma,                     
           width = 8, height = 4)        # inches
    
    ggsave(paste0("./result/down_tail/P_xi_exc/",curr.coin,".pdf"),       
           plot = neg_POT_plot$P_xi_exc,                     
           width = 8, height = 4)        # inches
    
    ggsave(paste0("./result/down_tail/P_std_sigma/",curr.coin,".pdf"),       
           plot = neg_POT_plot$P_std_sigma,                     
           width = 8, height = 4)        # inches
    
    #save(pos_POT_plot, file = paste0("./result/plot_handle/up_tail/",curr.coin,".RData"))
    #save(neg_POT_plot, file = paste0("./result/plot_handle/down_tail/",curr.coin,".RData"))
  }
  }, silent = TRUE)
  
  if (inherits(try.result, "try-error")) {
    message(paste("Error in", curr.coin, "— skipping to next"))
    next
  }
}

# Stop cluster
if (!plot.mode){
  stopCluster(cl)}