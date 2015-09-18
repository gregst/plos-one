########################################################################################
## Counts how many times specific features occured in the results from stabilityExp.R ##
########################################################################################

fname = "exp3\\stabilityExp3"
stabfile = stabfile <- paste0("c:\\temp\\resultsTemp\\", fname)

ntf = 5
i<-1
resSE <- readRDS(paste0(stabfile, ".NTF.", ntf, ".i.",i,".rds"))[["se"]]
resMIN <- readRDS(paste0(stabfile, ".NTF.", ntf, ".i.",i,".rds"))[["min"]]
for(i in 2:1000) {
  resTemp = NULL
  try(resTemp <- readRDS(paste0(stabfile, ".NTF.", ntf, ".i.",i,".rds")))
  resSE <- rBind(resSE, resTemp[["se"]])
  resMIN <- rBind(resMIN, resTemp[["min"]])
}

resSE <- as.matrix(resSE)
resMIN <- as.matrix(resMIN)

rankPosSE <- table(names(resSE[resSE>0,]))
max <- rankPosSE[order(rankPosSE, decreasing = T)][1]
ranked <- as.data.frame(rankPosSE[order(rankPosSE, decreasing = T)]/max)
colnames(ranked) <- "Final set presence (%)"
write.csv(ranked, file = "RankingPos1SE.csv")

rankNegSE <- table(names(resSE[resSE<0,]))
max <- rankNegSE[order(rankNegSE, decreasing = T)][1]
ranked <- as.data.frame(rankNegSE[order(rankNegSE, decreasing = T)]/max)
colnames(ranked) <- "Final set presence (%)"
write.csv(ranked, file = "RankingNeg1SE.csv")

rankPosMIN <- table(names(resMIN[resMIN>0,]))
max <- rankPosMIN[order(rankPosMIN, decreasing = T)][1]
ranked <- as.data.frame(rankPosMIN[order(rankPosMIN, decreasing = T)]/max)
colnames(ranked) <- "Final set presence (%)"
write.csv(ranked, file = "RankingPosMIN.csv")

rankNegMIN <- table(names(resMIN[resMIN<0,]))
max <- rankNegMIN[order(rankNegMIN, decreasing = T)][1]
ranked <- as.data.frame(rankNegMIN[order(rankNegMIN, decreasing = T)]/max)
colnames(ranked) <- "Final set presence (%)"
write.csv(ranked, file = "RankingNegMIN.csv")
