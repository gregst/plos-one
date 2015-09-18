library(glinternet)
library(glmnet)
library(pROC)
library(parallel)
library(doParallel)
library(C50)

#################################################################
## Classification performance test for GLI, OPT, 1SE and C5.0. ##
## Supports multi-threading for lower computational times.     ##
#################################################################

data <- readRDS("pediatricSID_CA_readmitYEAR.rds")

yR30 <- data[["y"]]$readmitBIN
data <- as.matrix(data[["x"]])

hids <- as.numeric(names(table(data[,"DSHOSPID"]))[table(data[,"DSHOSPID"])>149])
keep <- data[,"DSHOSPID"] %in% hids
data <- data[keep,]
yR30 <- yR30[keep]

# Keep only ICDs, AGEMONTH, FEMALE, LOS, NCHRONIC, NPR and TOTCHG
keep <- colnames(data) %in% c("AGEMONTH", "FEMALE", "LOS", "NCHRONIC", "NPR", "TOTCHG", "TOTCHG.NA", "TOTCHG_LOG", "LOS_LOG")
# Use only 213 (5%) of active ICD codes
keep[1:213] <- TRUE 
data <- data[,keep]

fname <- "exp\\stabilityExp"

runExp <- function(i) {
  logfile <- paste0("c:\\temp\\resultsTemp\\", fname, "_run", i, ".log")
  stabfile <- paste0("c:\\temp\\resultsTemp\\", fname)
  if (file.exists(logfile)) file.remove(logfile)
  cat(paste0("Timestamp, Model, Iteration, NumToFind, AUC, Threshold, ACC, Sensitivity, Specificity, PPV, NPV, Time, MComplexity, NumInteractions, Brier", "\n"), 
      file=logfile, append=T)
  
  set.seed(i)
  trainInt <- sample(1:nrow(data), as.integer(nrow(data)/3*2))
  X <- data[trainInt,]
  Y <- as.integer(yR30[trainInt])
  
  Xt <- as.matrix(data[-trainInt,])
  Yt <- as.integer(yR30[-trainInt])
  
  numLevels <- rep(1, ncol(X))
  for(j in 1:ncol(X)) {
    if(length(unique(X[,j])) == 2) {
      numLevels[j] <- 2
    }
  }
  
  ## Evaluate LASSO (via GLMNET)
  evalGLMNET <- function(X, Y, Xt, Yt, alpha, i, ntf, iLASSO) {
    ptm <- proc.time()
    # build model
    fit <- cv.glmnet(X, Y, family="binomial", type="auc", nfolds=5, alpha=alpha)    

    time <- proc.time()-ptm
    
    numSelected1SE <- fit$nzero[which(fit$lambda == fit$lambda.1se)]
    numSelectedMin <- fit$nzero[which(fit$lambda == fit$lambda.min)]
    #predict
    pred1SE <- predict(fit, type="response", Xt, s="lambda.1se")
    predMin <- predict(fit, type="response", Xt, s="lambda.min")
    
    # Count Int coefficients (for iLASSO)
    names <- names(coef(fit, s = "lambda.min")[as.vector(coef(fit, s = "lambda.min") != 0),])
    numIntMIN <- sum(substr(names, 1, 3)=="Int")
    names <- names(coef(fit, s = "lambda.1se")[as.vector(coef(fit, s = "lambda.1se") != 0),])
    numInt1SE <- sum(substr(names, 1, 3)=="Int")
    
    # Calculate Metrics 1SE
    brier <- mean((pred1SE-Yt)^2)
    mcom <- numSelected1SE
    result.roc <- roc(Yt, pred1SE)
    
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", 
                            ret=c("threshold", "accuracy", "sensitivity", "specificity", "ppv", "npv"))
    name <- "LASSO (1SE)"
    if(iLASSO) name <- "iLASSO (1SE)"
    str <- paste(Sys.time(), name, i, ntf, result.roc$auc, paste0(result.coords, collapse=","), time[1], mcom, numInt1SE, brier, sep=",")
    cat(paste0(str,"\n"), file=logfile, append=T)
    print(str)
    
    # Calculate Metrics MIN
    brier <- mean((predMin-Yt)^2)
    mcom <- numSelectedMin
    result.roc <- roc(Yt, predMin)
    
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", 
                            ret=c("threshold", "accuracy", "sensitivity", "specificity", "ppv", "npv"))
    name <- "LASSO (MIN)"
    if(iLASSO) name <- "iLASSO (MIN)"
    str <- paste(Sys.time(), name, i, ntf, result.roc$auc, paste0(result.coords, collapse=","), time[1], mcom, numIntMIN, brier, sep=",")
    cat(paste0(str,"\n"), file=logfile, append=T)
    print(str)
    
    # Dump coefficients for later feature selection stability calculations
    stab <- list(min=coef(fit, s = "lambda.min"), se=coef(fit, s = "lambda.1se"))
    saveRDS(stab, paste0(stabfile, ".NTF.", ntf, ".i.",i,".rds"))
    
    # Build decision tree
    ptm <- proc.time()
    treeModel <- C5.0(X, as.factor(Y), trials = 10, rules = T)
    
    time <- proc.time()-ptm
      
    numSelected <- treeModel$size
    pred <- predict(treeModel, type="prob", Xt)[,2]
      
    # Calculate Metrics for C5.0 tree
    brier <- mean((pred-Yt)^2)
    mcom <- numSelected
    numInt <- 0
    result.roc <- roc(Yt, pred)
      
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", 
                            ret=c("threshold", "accuracy", "sensitivity", "specificity", "ppv", "npv"))
      
    str <- paste(Sys.time(), "C5.0", i, ntf, result.roc$auc, paste0(result.coords, collapse=","), time[1], mcom, numInt, brier, sep=",")
    cat(paste0(str,"\n"), file=logfile, append=T)
    print(str)
      
    return(fit)
  }
  
  
  #  fitGLMNET <- NULL
  ## Evaluate GLINTERNET and InterLASSO
  ## for numToFind = [5,10,15,20]
  for(ntf in c(5,10,15,20)) {
    
    # Run GLINTERNET
    ptm <- proc.time()
    
    fit<-glinternet(X, Y, numLevels, family="binomial", lambdaMinRatio=0.01, verbose=T, numToFind=ntf, screenLimit=0.5)

    time <- proc.time()-ptm
    numL <- length(fit[["lambda"]])
    coefs = coef(fit)[[numL]]
    
    pred <- predict(fit, Xt, type="response", lambda = fit[["lambda"]][numL])
    
    # Calculate Metrics
    brier <- mean((pred-Yt)^2)
    mcom <- length(fit$betahat[[numL]])-1
    numInt <- sum(nrow(coefs$interactions$catcat)*4, 
                  nrow(coefs$interactions$contcont), 
                  nrow(coefs$interactions$catcont), na.rm=T)
    
    result.roc <- roc(Yt, pred[,ncol(pred)])
    
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", 
                            ret=c("threshold", "accuracy", "sensitivity", "specificity", "ppv", "npv"))
    str <- paste(Sys.time(), "GLINTERNET", i, ntf, result.roc$auc, paste0(result.coords, collapse=","), time[1], mcom, numInt, brier, sep=",")
    cat(paste0(str,"\n"), file=logfile, append=T)
    print(str)
    
    # Prepare data for InterLASSO
    catX <- catXt <- contX <- contXt <- NULL
    
    # Prepare main features from GLINTERNET
    if(!is.null(coefs$mainEffects$cat)) {
      catX <- X[,numLevels==2][,coefs$mainEffects$cat]  
      catXt <- Xt[,numLevels==2][,coefs$mainEffects$cat]      
    }
    if(!is.null(coefs$mainEffects$cont)) {
      contX <- X[,numLevels==1][,coefs$mainEffects$cont]  
      contXt <- Xt[,numLevels==1][,coefs$mainEffects$cont]  
    }
    
    # Prepare interactions from GLINTERNET
    if(!is.null(coefs$interactions$catcat)) {
      for(k in 1:nrow(coefs$interactions$catcat)) {
        x1 <- X[,numLevels==2][,coefs$interactions$catcat[k,1]]
        x2 <- X[,numLevels==2][,coefs$interactions$catcat[k,2]]
        xt1 <- Xt[,numLevels==2][,coefs$interactions$catcat[k,1]]
        xt2 <- Xt[,numLevels==2][,coefs$interactions$catcat[k,2]]
        
        name1 <- colnames(X[,numLevels==2])[coefs$interactions$catcat[k,1]]
        name2 <- colnames(X[,numLevels==2])[coefs$interactions$catcat[k,2]]
        
        name <- paste0("Int.", name1, "=1.", name2, "=0")
        catX <- cbind(catX, (x1==1) * (x2==0)) 
        catXt <- cbind(catXt, (xt1==1) * (xt2==0)) 
        colnames(catX)[ncol(catX)] <- colnames(catXt)[ncol(catXt)] <- name
        
        name <- paste0("Int.", name1, "=0.", name2, "=1")
        catX <- cbind(catX, (x1==0) * (x2==1)) 
        catXt <- cbind(catXt, (xt1==0) * (xt2==1)) 
        colnames(catX)[ncol(catX)] <- colnames(catXt)[ncol(catXt)] <- name      
        
        name <- paste0("Int.", name1, "=1.", name2, "=1")
        catX <- cbind(catX, (x1==1) * (x2==1)) 
        catXt <- cbind(catXt, (xt1==1) * (xt2==1)) 
        colnames(catX)[ncol(catX)] <- colnames(catXt)[ncol(catXt)] <- name      
      }
    }
    
    if(!is.null(coefs$interactions$contcont)) {
      for(k in 1:nrow(coefs$interactions$contcont)) {
        x1 <- X[,numLevels==1][,coefs$interactions$contcont[k,1]]
        x2 <- X[,numLevels==1][,coefs$interactions$contcont[k,2]]
        xt1 <- Xt[,numLevels==1][,coefs$interactions$contcont[k,1]]
        xt2 <- Xt[,numLevels==1][,coefs$interactions$contcont[k,2]]
        
        name1 <- colnames(X[,numLevels==1])[coefs$interactions$contcont[k,1]]
        name2 <- colnames(X[,numLevels==1])[coefs$interactions$contcont[k,2]]
        
        name <- paste0("Int.", name1, ".", name2)
        contX <- cbind(contX, x1 * x2) 
        contXt <- cbind(contXt, xt1 * xt2) 
        colnames(contX)[ncol(contX)] <- colnames(contXt)[ncol(contXt)] <- name
      }
    }
    
    if(!is.null(coefs$interactions$catcont)) {
      for(k in 1:nrow(coefs$interactions$catcont)) {
        x1 <- X[,numLevels==2][,coefs$interactions$catcont[k,1]]
        x2 <- X[,numLevels==1][,coefs$interactions$catcont[k,2]]
        xt1 <- Xt[,numLevels==2][,coefs$interactions$catcont[k,1]]
        xt2 <- Xt[,numLevels==1][,coefs$interactions$catcont[k,2]]
        
        name1 <- colnames(X[,numLevels==2])[coefs$interactions$catcont[k,1]]
        name2 <- colnames(X[,numLevels==1])[coefs$interactions$catcont[k,2]]
        
        name <- paste0("Int.", name1, ".", name2)
        catX <- cbind(catX, x1 * x2) 
        catXt <- cbind(catXt, xt1 * xt2) 
        colnames(catX)[ncol(catX)] <- colnames(catXt)[ncol(catXt)] <- name
      }
    }
    
    # Run InterLASSO
    evalGLMNET(cbind(catX, contX), Y, cbind(catXt, contXt), Yt, 1, i, ntf, TRUE)
  }
  
}

numWorkers <- 12
cl <- makeCluster(numWorkers)
registerDoParallel(cl)
system.time(foreach(i=1:1000, .packages=c('glmnet', 'pROC', 'glinternet', 'C50')) %dopar% runExp(i))  
stopCluster(cl)

allData <- data.frame()
for(i in 1:1000) {
  logfile <- paste0("c:\\temp\\resultsTemp\\", fname, "_run", i, ".log")
  allData <- rbind(allData, read.csv(logfile))  
}
write.csv(allData, paste0("c:\\temp\\resultsTemp\\", fname, "_all.log"))
