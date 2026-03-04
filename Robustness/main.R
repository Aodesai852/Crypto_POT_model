rm(list=ls())

setwd("D:/MyFiles/EVT/Crypto_POT")
library("igraph")
library("data.table")
library("ggplot2")
source("./code/Robustness/func_robustness_plot.R")

dca <- new.env()
load("./result/DY/percent95win200/dca.RData", env = dca)
dca$net_down = dca$dca_down$TABLE["NET",]
dca$net_down = dca$net_down[-length(dca$net_down)]
dca$net_up = dca$dca_up$TABLE["NET",]
dca$net_up = dca$net_up[-length(dca$net_up)]

dca_rob90 <- new.env()
load("./result/DY/percent90win200/dca.RData", env = dca_rob90)
dca_rob90$net_down = dca_rob90$dca_down$TABLE["NET",]
dca_rob90$net_down = dca_rob90$net_down[-length(dca_rob90$net_down)]
dca_rob90$net_up = dca_rob90$dca_up$TABLE["NET",]
dca_rob90$net_up = dca_rob90$net_up[-length(dca_rob90$net_up)]

dca_rob85 <- new.env()
load("./result/DY/percent85win200/dca.RData", env = dca_rob85)
dca_rob85$net_down = dca_rob85$dca_down$TABLE["NET",]
dca_rob85$net_down = dca_rob85$net_down[-length(dca_rob85$net_down)]
dca_rob85$net_up = dca_rob85$dca_up$TABLE["NET",]
dca_rob85$net_up = dca_rob85$net_up[-length(dca_rob85$net_up)]

dca_win150 <- new.env()
load("./result/DY/percent95win150/dca.RData", env = dca_win150)
dca_win150$net_down = dca_win150$dca_down$TABLE["NET",]
dca_win150$net_down = dca_win150$net_down[-length(dca_win150$net_down)]
dca_win150$net_up = dca_win150$dca_up$TABLE["NET",]
dca_win150$net_up = dca_win150$net_up[-length(dca_win150$net_up)]

dca_win250 <- new.env()
load("./result/DY/percent95win250/dca.RData", env = dca_win250)
dca_win250$net_down = dca_win250$dca_down$TABLE["NET",]
dca_win250$net_down = dca_win250$net_down[-length(dca_win250$net_down)]
dca_win250$net_up = dca_win250$dca_up$TABLE["NET",]
dca_win250$net_up = dca_win250$net_up[-length(dca_win250$net_up)]


lst_up <- list(bm=dca$net_up,rob90=dca_rob90$net_up,rob85=dca_rob85$net_up,
               win150=dca_win150$net_up,win250=dca_win250$net_up)
lst_down <- list(bm=dca$net_down,rob90=dca_rob90$net_down,rob85=dca_rob85$net_down,
               win150=dca_win150$net_down,win250=dca_win250$net_down)

DT_up <- vecs_to_dt(lst_up)
DT_down <- vecs_to_dt(lst_down)

DT_up[, (names(DT_up)[-1]) := lapply(.SD, as.numeric), .SDcols = -1]
DT_down[, (names(DT_down)[-1]) := lapply(.SD, as.numeric), .SDcols = -1]

DT_up[, (names(DT_up)[-1]) := lapply(.SD, frank, ties.method = "average"),
   .SDcols = -1]
DT_down[, (names(DT_down)[-1]) := lapply(.SD, frank, ties.method = "average"),
      .SDcols = -1]



lst_TCI_up <- list(dca = dca$dca_up$TCI, dca_rob90 = dca_rob90$dca_up$TCI, dca_rob85 = dca_rob85$dca_up$TCI,
                   dca_win150 = dca_win150$dca_up$TCI, dca_win250 = dca_win250$dca_up$TCI)
lst_TCI_down <- list(dca = dca$dca_down$TCI, dca_rob90 = dca_rob90$dca_down$TCI, dca_rob85 = dca_rob85$dca_down$TCI,
                   dca_win150 = dca_win150$dca_down$TCI, dca_win250 = dca_win250$dca_down$TCI)

DT_TCI_up <- na.omit(matrices_to_dt_by_time(lst_TCI_up))
DT_TCI_down <- na.omit(matrices_to_dt_by_time(lst_TCI_down))


plotup_rob90 <- rank45_plot(DT_up[,c(1,2,3)])
plotup_rob85 <- rank45_plot(DT_up[,c(1,2,4)])
plotup_win150 <- rank45_plot(DT_up[,c(1,2,5)])
plotup_win250 <- rank45_plot(DT_up[,c(1,2,6)])

plotdown_rob90 <- rank45_plot(DT_down[,c(1,2,3)])
plotdown_rob85 <- rank45_plot(DT_down[,c(1,2,4)])
plotdown_win150 <- rank45_plot(DT_down[,c(1,2,5)])
plotdown_win250 <- rank45_plot(DT_down[,c(1,2,6)])

TCI_plot_up_rob = plot_multi_series(DT_TCI_up[,c("time","dca_TCI","dca_rob90_TCI","dca_rob85_TCI")])
TCI_plot_up_win = plot_multi_series(DT_TCI_up[,c("time","dca_TCI","dca_win150_TCI","dca_win250_TCI")])
TCI_plot_down_rob = plot_multi_series(DT_TCI_down[,c("time","dca_TCI","dca_rob90_TCI","dca_rob85_TCI")])
TCI_plot_down_win = plot_multi_series(DT_TCI_down[,c("time","dca_TCI","dca_win150_TCI","dca_win250_TCI")])


# change labels

TCI_plot_up_win <- TCI_plot_up_win + labs(
  x = NULL,
  y = "Total risk spillover index"
) + scale_color_discrete(
  labels = c(
    "window size = 200 days",
    "window size = 150 days",
    "window size = 250 days"
  )
) 

TCI_plot_down_win <- TCI_plot_down_win + labs(
  x = NULL,
  y = "Total risk spillover index"
) + scale_color_discrete(
  labels = c(
    "window size = 200 days",
    "window size = 150 days",
    "window size = 250 days"
  )
) 

TCI_plot_up_rob <- TCI_plot_up_rob + labs(
  x = NULL,
  y = "Total risk spillover index"
) + scale_color_discrete(
  labels = c(
    "threshold = 0.95",
    "threshold = 0.90",
    "threshold = 0.85"
  )
) 

TCI_plot_down_rob <- TCI_plot_down_rob + labs(
  x = NULL,
  y = "Total risk spillover index"
) + scale_color_discrete(
  labels = c(
    "threshold = 0.95",
    "threshold = 0.90",
    "threshold = 0.85"
  )
) 

plotup_rob90 <- plotup_rob90 + labs(
  x = "threshold = 0.90",
  y = "threshold = 0.95" # benchmark
)

plotup_rob85 <- plotup_rob85 + labs(
  x = "threshold = 0.85",
  y = "threshold = 0.95" # benchmark
)

plotup_win150 <- plotup_win150 + labs(
  x = "window size = 150",
  y = "window size = 200" # benchmark
)

plotup_win250 <- plotup_win250 + labs(
  x = "window size = 250",
  y = "window size = 200" # benchmark
)

plotdown_rob90 <- plotdown_rob90 + labs(
  x = "threshold = 0.90",
  y = "threshold = 0.95" # benchmark
)

plotdown_rob85 <- plotdown_rob85 + labs(
  x = "threshold = 0.85",
  y = "threshold = 0.95" # benchmark
)

plotdown_win150 <- plotdown_win150 + labs(
  x = "window size = 150",
  y = "window size = 200" # benchmark
)

plotdown_win250 <- plotdown_win250 + labs(
  x = "window size = 250",
  y = "window size = 200" # benchmark
)

# save the plots
ggsave("./result/Robustness/TCI_plot/rob_down.pdf",       
       plot = TCI_plot_down_rob,                     
       width = 10, height = 3)        # inches
ggsave("./result/Robustness/TCI_plot/rob_up.pdf",       
       plot = TCI_plot_up_rob,                     
       width = 10, height = 3)        # inches
ggsave("./result/Robustness/TCI_plot/win_down.pdf",       
       plot = TCI_plot_down_win,                     
       width = 10, height = 3)        # inches
ggsave("./result/Robustness/TCI_plot/win_up.pdf",       
       plot = TCI_plot_up_win,                     
       width = 10, height = 3)        # inches

ggsave("./result/Robustness/rank45plot/rob90down.pdf",
       plot = plotdown_rob90,
       width = 8, height = 4)
ggsave("./result/Robustness/rank45plot/rob85down.pdf",
       plot = plotdown_rob85,
       width = 8, height = 4)
ggsave("./result/Robustness/rank45plot/win150down.pdf",
       plot = plotdown_win150,
       width = 8, height = 4)
ggsave("./result/Robustness/rank45plot/win250down.pdf",
       plot = plotdown_win250,
       width = 8, height = 4)

ggsave("./result/Robustness/rank45plot/rob90up.pdf",
       plot = plotup_rob90,
       width = 8, height = 4)
ggsave("./result/Robustness/rank45plot/rob85up.pdf",
       plot = plotup_rob85,
       width = 8, height = 4)
ggsave("./result/Robustness/rank45plot/win150up.pdf",
       plot = plotup_win150,
       width = 8, height = 4)
ggsave("./result/Robustness/rank45plot/win250up.pdf",
       plot = plotup_win250,
       width = 8, height = 4)






