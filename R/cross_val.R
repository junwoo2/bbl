#' Cross-validation run
#' @export
crossval <- function(object, lambda=0.1, nfold=5, method='pseudo', 
                     eps=0.9, L=NULL, auc=TRUE, naive=FALSE, verbose=1, 
                     computeZ=FALSE, mf=TRUE, useC=TRUE, prior.count=TRUE){
  
  groups <- object@groups
  Ly <- length(groups)
  nsample <- NROW(object@data)
  y <- object@data[,which(object@y==colnames(object@data))]
  
  if(Ly!=2) auc <- FALSE
  
  if(method=='pseudo') reglist <- lambda
  else if(method=='mf') reglist <- eps
  else if(method=='nb') reglist <- eps <- 0
  else stop('Unknown method')
  
  res <- NULL
  for(reg in reglist){
    if(verbose>0){ 
      if(method=='pseudo') cat('Cross validation under lambda = ',sep='')
      else cat('Cross validation under epsilon = ',sep='')
        cat(reg,'\n',sep='')
    }
    pred <- NULL
    for(k in seq_len(nfold)){
      if(verbose>0) cat(' Fold no. ',k,'...\n',sep='')
      itrain <- ival <- NULL
      for(iy in seq_len(Ly)){
        idy <- which(y==groups[iy])
        ns <- length(idy)
        nval <- max(1,floor(ns/nfold))
        imax <- k*nval
        if(imax > ns) imax <- ns
        if(imax < 1) imax <- 1
        iyval <- idy[seq((k-1)*nval+1, imax)]
        iytrain <- idy[!idy %in% iyval]
        ival <- c(ival, iyval)
        itrain <- c(itrain, iytrain)
      }
      if(sum(is.na(ival)>0) | sum(is.na(itrain)>0)) stop('error in crossval')
      obval <- object[ival,]
      obtrain <- object[itrain,]
      if(method=='pseudo')
        obtrain <- train(object=obtrain, L=L, lambda=reg, method=method, 
                         prior.count=prior.count, verbose=verbose-1)
      else{
        obtrain <- train(object=obtrain, L=L, eps=reg, method=method, 
                         prior.count=prior.count, verbose=verbose-1)
        computeZ <- TRUE
      } 
      pr <- predict(object=obtrain, newdata=obval@data, logit=TRUE,
                    L=L, computeZ=computeZ, mf=mf, useC=useC)
      pred <- rbind(pred, cbind(data.frame(y=y[ival], pr)))
    }
    if(auc){ auc <- pROC::roc(response=pred$y, levels=groups, 
                               predictor=pred[,3], direction='<')$auc
      if(verbose>0) cat(' AUC = ',auc,'\n',sep='')
    }
    else{
      yhat <- groups[apply(pred[,-1],1,which.max)]
      score <- mean(pred$y==yhat)
      if(verbose>0) cat(' Prediction score = ',score,'\n',sep='')
    }
    if(method=='pseudo'){
      if(auc) rx <- data.frame(lambda=reg, auc=auc)
      else rx <- data.frame(lambda=reg, score=score)
    } else{
      if(auc) rx <- data.frame(epsilon=reg, auc=auc)
      else rx <- data.frame(epsilon=reg, score=score)
    }
    res <- rbind(res, rx)
  }
  
  return(res)
}