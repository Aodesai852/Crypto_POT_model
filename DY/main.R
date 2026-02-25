rm(list=ls())
library(data.table)
library(stats)
library(zoo)
library(rugarch)
library(rmgarch)
library(progress)
setwd("D:/MyFiles/EVT/BTC_POT/code")
source("./Risk_spillovers.R")
core.mode = FALSE # only estimate core currencies
rob.mode = FALSE # robustness test, if = FALSE, q = 0.95
rob.q = 0.85 # 0.90 or 0.85
win.size.rob = FALSE # alter the window size to see whether it's robust
# if win.size.rob = FALSE, window.size must be 200
window.size = 200 # 150, 200 or 250

if (rob.mode){
  addr_down = paste0("./result/rob_test/",sprintf("%.2f",rob.q),"/down_tail/optim_result/")
  addr_up = paste0("./result/rob_test/",sprintf("%.2f",rob.q),"/up_tail/optim_result/")
} else {
  addr_down = "./result/down_tail/optim_result/"
  addr_up = "./result/up_tail/optim_result/"
}

coin.names = c(
"AAVE","ALGO","AVAX","BNB","BCH","BTC","ADA","LINK","CRO","DOGE","ETC","ETH","FIL","HBAR",
"KCS","XMR","NEAR","OKB","DOT","QNT","RNDR","SOL","XLM","TRX","UNI","LEO","VET","ZEC"
)

core.coin.names = c(
  "BTC","ETH"
)

zoo.up.combined.list = list()
zoo.down.combined.list = list()

if (core.mode){
  names = core.coin.names
}else{
  names = coin.names
}

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

dca_up = ConnectednessApproach(zoo.up.combined,
                            nlag=1,
                            nfore=20,
                            model="VAR",
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
                               model="VAR",
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
if (rob.mode){
  if (core.mode){
    save(zoo.up.combined,zoo.down.combined,dca_up,dca_down, file = paste0("./result/DY/rob_test/",sprintf("%.2f", rob.q),"/dca_core.RData"))
  }else{
    save(zoo.up.combined,zoo.down.combined,dca_up,dca_down, file = paste0("./result/DY/rob_test/",sprintf("%.2f", rob.q),"/dca.RData"))
  }
} else {
  if (win.size.rob){
    if (core.mode){
      save(zoo.up.combined,zoo.down.combined,dca_up,dca_down, file = paste0("./result/DY/win_size_rob/",sprintf("%d",window.size),"/dca_core.RData"))
    }else{
      save(zoo.up.combined,zoo.down.combined,dca_up,dca_down, file = paste0("./result/DY/win_size_rob/",sprintf("%d",window.size),"/dca.RData"))
    }
  } else {
    if (core.mode){
      save(zoo.up.combined,zoo.down.combined,dca_up,dca_down, file = "./result/DY/dca_core.RData")
    }else{
      save(zoo.up.combined,zoo.down.combined,dca_up,dca_down, file = "./result/DY/dca.RData")
    }
  }
}

# zoo.down.combined and zoo.up.combined are the reciprocals of the tail index of all currencies
# dca_up and dca_down are the static and dynamic DY results 

