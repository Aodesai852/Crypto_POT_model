rm(list = ls())
setwd("D:/MyFiles/EVT/BTC_POT/code/data")
library("data.table")

stable.coins <- c(
  "tether",
  "usd-coin",
  "dai",
  "paypal-usd",
  "tether-gold",
  "pax-gold",
  "usds",
  "global-dollar-usdg",
  "ethena-usde",
  "onus-sed-leo",
  "usd1-wlfi",
  "decentralized-usd-stablecoin"
)

files <- list.files("./ori_data", pattern = "\\.csv$", full.names = TRUE)

coins <- sub("_.*", "", basename(files))

start.dates <- as.Date(
  sub(".*_(\\d{4}-\\d{2}-\\d{2})_.*", "\\1", basename(files))
)

keep <- start.dates < as.Date("2021-01-01") &
  !coins %in% stable.coins

files.filtered <- files[keep]

dt_list <- lapply(files.filtered, function(f) {
  coin <- sub("_.*", "", basename(f))  # coin name from filename
  
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