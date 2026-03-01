# check if two folders have the same files
rm(list = ls())
setwd("D:/MyFiles/EVT/Crypto_POT")

folderA <- "./result/POT/percent95/down_tail/optim_result"
folderB <- "./result/POT/percent95/up_tail/optim_result"

filesA <- list.files(folderA)
filesB <- list.files(folderB)

onlyA <- setdiff(filesA, filesB)

onlyB <- setdiff(filesB, filesA)

cat("A only:\n")
print(onlyA)

cat("\nB only:\n")
print(onlyB)