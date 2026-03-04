rm(list=ls())

library(data.table)
library(stats)
library(zoo)
library(rugarch)
library(rmgarch)
library(progress)

setwd("D:/MyFiles/EVT/Crypto_POT")
source("./code/DY/Risk_spillovers.R")

######## changable ########################
core.mode = FALSE
model = "VAR"
window.size = 200 # 150/200/250
addr_down = "./result/POT/percent95/down_tail/optim_result/"
addr_up = "./result/POT/percent95/up_tail/optim_result/"
target_folder = "D:/MyFiles/EVT/Crypto_POT/result/DY/percent95win200/"
coin.names = c("AAVE","ADA","ALGO","AVAX","AXS","BAT","BCH","BNB","BTC","CAKE","CRO","CRV","DASH",
               "DCR","DEXE","DOGE","DOT","ETC","ETH","FET","FIL","GNO","GT","HBAR","HNT","INJ","IOTA","JST",
               "KCS","LEO","LINK","LTC","MANA","NEAR","NEO","OKB","QNT","RNDR","SOL","STX","TEL","TRX",
               "TWT","UNI","VET","XLM","XMR","XRP","XTZ","ZEC")
core.coin.names = c("BCH","BNB","BTC","DOGE","ETC","ETH","SOL","TRX","XRP","ZEC")
###########################################

if (core.mode){
  names = core.coin.names
}else{
  names = coin.names
}

zoo.up.combined.list = list()
zoo.down.combined.list = list()

for (curr.coin in names){
  load(paste0(addr_down,curr.coin,".RData"))
  load(paste0(addr_up,curr.coin,".RData"))
  uptail_df = pos_optim_result$xiSigma_timeseries_df[,c("date","xi")] # upside tail index time series
  uptail_df$xi = 1 / uptail_df$xi
  uptail_zoo <- zoo(uptail_df$xi, order.by = uptail_df$date)
  downtail_df = neg_optim_result$xiSigma_timeseries_df[,c("date","xi")] # downside tail index time series
  downtail_df$xi = 1 / downtail_df$xi
  downtail_zoo <- zoo(downtail_df$xi, order.by = downtail_df$date)
  zoo.up.combined.list[[curr.coin]] = uptail_zoo
  zoo.down.combined.list[[curr.coin]] = downtail_zoo
}

zoo.up.combined <- na.omit(do.call(merge, c(zoo.up.combined.list, all = TRUE)))
zoo.down.combined <- na.omit(do.call(merge, c(zoo.down.combined.list, all = TRUE)))
# zoo.down.combined and zoo.up.combined are the reciprocals of the tail index of all currencies

dca_up = ConnectednessApproach(zoo.up.combined,
                            nlag=1,
                            nfore=20,
                            model=model,
                            corrected=FALSE,
                            window.size=window.size,
                            connectedness="Time",
                            VAR_config=list(
                              QVAR=list(tau=0.2),
                              ElasticNet=list(nfolds=4, alpha=0, loss="mae", delta_alpha=0.1),
                              TVPVAR=list(kappa1=0.99, kappa2=0.99, prior="BayesPrior", gamma=0.01)),
                            DCC_config=list(standardize=FALSE),
                            Connectedness_config = list(
                              TimeConnectedness=list(generalized=TRUE),
                              FrequencyConnectedness=list(partition=c(pi,pi/2,0), generalized=TRUE, scenario="ABS")
                            ))

dca_down = ConnectednessApproach(zoo.down.combined,
                               nlag=1,
                               nfore=20,
                               model=model,
                               corrected=FALSE,
                               window.size=window.size,
                               connectedness="Time",
                               VAR_config=list(
                                 QVAR=list(tau=0.2),
                                 ElasticNet=list(nfolds=4, alpha=0, loss="mae", delta_alpha=0.1),
                                 TVPVAR=list(kappa1=0.99, kappa2=0.99, prior="BayesPrior", gamma=0.01)),
                               DCC_config=list(standardize=FALSE),
                               Connectedness_config = list(
                                 TimeConnectedness=list(generalized=TRUE),
                                 FrequencyConnectedness=list(partition=c(pi,pi/2,0), generalized=TRUE, scenario="ABS")
                               ))

# save the results

if (core.mode){
  save(dca_up,dca_down, file = paste0(target_folder,"dca_core.RData"))
}else{
  save(dca_up,dca_down, file = paste0(target_folder,"dca.RData"))
}


