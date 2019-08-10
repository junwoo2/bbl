#' Cross-validation run
#' @export
crossval <- function(object, lambda=0, nfold=5, L=NULL, auc=TRUE, naive=FALSE,
                     verbose=1){
  
  groups <- object@groups
  Ly <- length(groups)
  nsample <- NROW(object@data)
  
  if(Ly!=2) auc <- FALSE
  
  res <- NULL
  for(xlambda in lambda){
    if(verbose>0){ 
      if(!naive)
        cat('Cross validation under lambda=',xlambda,' ...\n',sep='')
      else
        cat('Cross validation with naive Bayes ...\n',sep='')
    }
    pred <- NULL
    for(k in seq_len(nfold)){
      if(verbose>0) cat(' Fold no. ',k,'...\n',sep='')
      itrain <- ival <- NULL
      for(iy in seq_len(Ly)){
        idy <- which(object@data$y==groups[iy])
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
      obtrain <- train(object=obtrain, lambda=xlambda, L=L, naive=naive,
                       verbose=verbose-1)
      pr <- predict(object=obtrain, newdata=obval@data, logit=TRUE)
      pred <- rbind(pred, cbind(data.frame(y=obval@data$y), pr))
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
    if(auc) rx <- data.frame(lambda=xlambda, auc=auc)
    else rx <- data.frame(lambda=xlambda, score=score)
    res <- rbind(res, rx)
  }
  
  return(res)
}