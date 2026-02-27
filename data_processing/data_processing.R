rm(list = ls())
setwd("D:/MyFiles/EVT/Crypto_POT/data")
library("data.table")

files <- list.files("./ori_data", pattern = "\\.csv$", full.names = TRUE)

coins <- sub(".csv", "", basename(files))

dt_list <- lapply(files, function(f) {
  coin <- sub(".csv", "", basename(f))  # coin name from filename
  
  dt <- fread(f, select = c("End", "Close"))
  
  # End -> Date
  dt[, End := as.Date(End, format = "%d/%m/%Y")]
  
  # Close -> numeric
  # (fread usually reads numeric already; this enforces it)
  dt[, Close := as.numeric(Close)]
  
  # rename Close column to coin name, keep End as key
  setnames(dt, "Close", coin)
  dt
})

panel.dt <- Reduce(function(x, y) merge(x, y, by = "End", all = TRUE), dt_list)

setorder(panel.dt, End)
setnames(panel.dt, "End", "Date")

save(panel.dt, file = "panel_dt.RData")