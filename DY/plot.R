rm(list=ls())

setwd("D:/MyFiles/EVT/Crypto_POT")
library("igraph")
library("zoo")
library("data.table")
source("./code/DY/Risk_spillovers_plot.R")

############## changable #################
threshold = 0.50
threshold_color = 0.75
load_addr = "./result/DY/percent95win250/dca.RData"
target_folder = "./result/DY/percent95win250/"
##########################################

load(load_addr)

suppressWarnings({
PlotNetwork(dca_down, path = paste0(target_folder,"network_plot_down.pdf"), threshold=threshold, threshold_color = threshold_color)
PlotNetwork(dca_up, path = paste0(target_folder,"network_plot_up.pdf"), threshold=threshold, threshold_color = threshold_color)
PlotNET(dca_down, path = paste0(target_folder,"NET_plot_down"), separate = TRUE)
PlotNET(dca_up, path = paste0(target_folder,"NET_plot_up"), separate = TRUE)
PlotTCI(dca_down, path = paste0(target_folder,"TCI_plot_down.pdf"), ylim=c(60,100))
PlotTCI(dca_up, path = paste0(target_folder,"TCI_plot_up.pdf"), ylim=c(60,100))
})