rm(list=ls())
setwd("D:/MyFiles/EVT/BTC_POT/code")
library("igraph")
library("zoo")
library("data.table")
source("./Risk_spillovers_plot.R")
threshold = 0.20
threshold_color = 0.50
rob.mode = FALSE

if (rob.mode) {
dca_rob85 <- new.env()
load("./result/DY/rob_test/0.85/dca.RData", envir = dca_rob85)
dca_rob85$reci.downtail.index <- data.table(
  time = index(dca_rob85$zoo.down.combined),
  as.data.frame(coredata(dca_rob85$zoo.down.combined))
)
dca_rob85$reci.uptail.index <- data.table(
  time = index(dca_rob85$zoo.up.combined),
  as.data.frame(coredata(dca_rob85$zoo.up.combined))
)


dca_rob90 <- new.env()
load("./result/DY/rob_test/0.90/dca.RData", envir = dca_rob90)
dca_rob90$reci.downtail.index <- data.table(
  time = index(dca_rob90$zoo.down.combined),
  as.data.frame(coredata(dca_rob90$zoo.down.combined))
)
dca_rob90$reci.uptail.index <- data.table(
  time = index(dca_rob90$zoo.up.combined),
  as.data.frame(coredata(dca_rob90$zoo.up.combined))
)

dca_win250 <- new.env()
load("./result/DY/win_size_rob/250/dca.RData", envir = dca_win250)

dca_win150 <- new.env()
load("./result/DY/win_size_rob/150/dca.RData", envir = dca_win150)
}

dca <- new.env()
load("./result/DY/dca.RData", envir = dca)
dca$reci.downtail.index <- data.table(
  time = index(dca$zoo.down.combined),
  as.data.frame(coredata(dca$zoo.down.combined))
)
dca$reci.uptail.index <- data.table(
  time = index(dca$zoo.up.combined),
  as.data.frame(coredata(dca$zoo.up.combined))
)

suppressWarnings({
PlotNetwork(dca$dca_down, path = "./result/down_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNetwork(dca$dca_up, path = "./result/up_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNET(dca$dca_down, path = "./result/down_tail/NET_plots", separate = TRUE)
PlotNET(dca$dca_up, path = "./result/up_tail/NET_plots", separate = TRUE)

if (rob.mode) {
PlotNetwork(dca_rob85$dca_down, path = "./result/rob_plot/0.85/down_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNetwork(dca_rob85$dca_up, path = "./result/rob_plot/0.85/up_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNET(dca_rob85$dca_down, path = "./result/rob_plot/0.85/down_tail/NET_plots", separate = TRUE)
PlotNET(dca_rob85$dca_up, path = "./result/rob_plot/0.85/up_tail/NET_plots", separate = TRUE)

PlotNetwork(dca_rob90$dca_down, path = "./result/rob_plot/0.90/down_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNetwork(dca_rob90$dca_up, path = "./result/rob_plot/0.90/up_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNET(dca_rob90$dca_down, path = "./result/rob_plot/0.90/down_tail/NET_plots", separate = TRUE)
PlotNET(dca_rob90$dca_up, path = "./result/rob_plot/0.90/up_tail/NET_plots", separate = TRUE)

PlotNetwork(dca_win150$dca_down, path = "./result/win_plot/150/down_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNetwork(dca_win150$dca_up, path = "./result/win_plot/150/up_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNET(dca_win150$dca_down, path = "./result/win_plot/150/down_tail/NET_plots", separate = TRUE)
PlotNET(dca_win150$dca_up, path = "./result/win_plot/150/up_tail/NET_plots", separate = TRUE)

PlotNetwork(dca_win250$dca_down, path = "./result/win_plot/250/down_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNetwork(dca_win250$dca_up, path = "./result/win_plot/250/up_tail/", threshold=threshold, threshold_color = threshold_color)
PlotNET(dca_win250$dca_down, path = "./result/win_plot/250/down_tail/NET_plots", separate = TRUE)
PlotNET(dca_win250$dca_up, path = "./result/win_plot/250/up_tail/NET_plots", separate = TRUE)
}
})