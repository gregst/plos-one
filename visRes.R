library(ggplot2)
library(plyr)

#############################################################
## Visualization of results stored in stabilityExp_all.log ##
#############################################################

data <- read.csv("c:/temp/resultsTemp/exp/stabilityExp_all.log", stringsAsFactors=FALSE)

# Merge C5.0 into single record (sum number of rules)
tempData <- data[data$Model=="C5.0",]
prevIter <- 1
prevNtf <- 5
sumRules <- 0
newData <- rbind(newData, tempData[1,])
for(i in 1:nrow(tempData)) {
  iter <- tempData[i,"Iteration"]
  ntf <- tempData[i,"NumToFind"]
  if(iter==prevIter & ntf==prevNtf) {
    sumRules <- sumRules + tempData$MComplexity[i]
  } else {
    prevIter <- iter
    prevNtf <- ntf 
    newData$MComplexity[nrow(newData)] <- sumRules
    if(newData$MComplexity[nrow(newData)] == 0) {
      newData[nrow(newData),] <- tempData[i,]
    } else {
      newData <- rbind(newData, tempData[i,])
    }
    sumRules <- 0
  }
}
newData$MComplexity[nrow(newData)] <- sumRules

newData <- rbind(newData, data[data$Model=="GLINTERNET" | data$Model=="iLASSO (1SE)" | data$Model=="iLASSO (MIN)",])

data <- newData

#data2 <- data[!is.na(data$AUC),c("Model", "NumToFind", "AUC")]
data2 <- data[!is.na(data$AUC),]
data2 <- data2[data2$Model %in% c("GLINTERNET", "iLASSO (MIN)", "iLASSO (1SE)", "C5.0"),]
data2$Model[data2$Model=="iLASSO (MIN)"] <- "OPT"
data2$Model[data2$Model=="iLASSO (1SE)"] <- "1SE"
data2$Model[data2$Model=="GLINTERNET"] <- "GLI"

colnames(data2)[14] <- "Complexity"

cdf3 <- ddply(data2, c("Model", "NumToFind"), summarise, auc.mean=mean(AUC), auc.sd=sd(AUC), auc.lq=quantile(AUC, 0.025), auc.hq=quantile(AUC, 0.975))
cdf3Size <- ddply(data2, c("Model", "NumToFind"), summarise, size.mean=mean(Complexity), size.lq=quantile(MComplexity, 0.025), size.hq=quantile(MComplexity, 0.975))


data2 <- data2[data2$NumToFind%%5==0,]
# AUC
ggplot(data2, aes(x=Model, y=AUC, fill=Model)) + geom_boxplot() +  facet_grid(. ~ NumToFind) + theme_bw() + ylim(0.72,0.79)
png("InteractionsAUC.png", width = 14, height = 7, units = 'in', res = 300)
ggplot(data2, aes(x=Model, y=AUC, fill=Model)) + geom_boxplot() +  facet_grid(. ~ NumToFind) + theme_bw() + ylim(0.72,0.79)
dev.off()


# Brier
ggplot(data2, aes(x=Model, y=Brier, fill=Model)) + geom_boxplot() +  facet_grid(. ~ NumToFind) + theme_bw()

# Size
ggplot(data2, aes(x=Model, y=Complexity, fill=Model)) + geom_boxplot() +  facet_grid(. ~ NumToFind) + theme_bw()
png("InteractionsComplexity.png", width = 10, height = 5, units = 'in', res = 300)
ggplot(data2, aes(x=Model, y=Complexity, fill=Model)) + geom_boxplot() +  facet_grid(. ~ NumToFind) + theme_bw()
dev.off()

