#' Cross-Validation of BB Learning
#' 
#' Run multiple fittings of \code{bbl} model with training/validation
#' division of data
#' 
#' The \code{data} slot of \code{object} is split into training and validation 
#' subsets of (\code{nfold}-1):1 ratio. The model is trained with the
#' former and validated on the latter. Individual division/fold results are 
#' combined into validation result for all instances in the data set and
#' prediction score is evaluated using the known response group
#' identity.
#' 
#' @param object Object of class \code{bbl} containing data.
#' @param lambda Vector of L2 penalizer values for \code{method = 'pseudo'}. Inferences
#'        will be repeated for each value. Restricited to non-negative values.
#' @param lambdah L2 penalizer in \code{method = 'pseudo'} applied to
#'        parameter \code{h}. In contrast to \code{lambda}, 
#'        only a single value is allowed.
#' @param eps Vector of regularization parameters, \eqn{\epsilon\in[0,1]}, 
#'        for \code{method = 'mf'}. Inference will be repeeated
#'        for each value.
#' @param nfold Number of folds for training/validation split.
#' @param method \code{c('pseudo','mf')} for pseudo-likelihood maximization or
#'        mean field.
#' @param naive Naive Bayes (no interactions). Equivalent to \code{method = 'mf'}
#'        together with \code{eps = 0}.
#' @param use.auc Use AUC as the measure of prediction accuracy. Only works
#'        if response groups are binary. If \code{FALSE}, mean prediction group
#'        accuracy will be used as score.
#' @param verbose Verbosity level. Downgraded when relayed into \code{\link{train}}.
#' @param useC Use \code{C++} version in \code{\link{predict}} method of \code{bbl}.
#' @param prior.count Use prior count in \code{method = 'mf'}.
#' @param progress.bar Display progress bar in \code{\link{predict}}.
#' @param fixL Do not alter the levels of predictors in training step.
#' @param ... Other parameters to \code{\link{mlestimate}}.
#' @return Data frame of regularization parameter values and validation scores.
#' @examples
#' set.seed(513)
#' m <- 5
#' n <- 100
#' predictors <- list()
#' for(i in 1:m) predictors[[i]] <- c('a','c','g','t')
#' 
#' par0 <- randompar(predictors)
#' xi0 <- sample_xi(nsample=n, predictors=predictors, h=par0$h, J=par0$J)
#' par1 <- randompar(predictors, h0=0.1, J0=0.1)
#' xi1 <- sample_xi(nsample=n, predictors=predictors, h=par1$h, J=par1$J)
#' xi <- rbind(xi0, xi1)
#' dat <- cbind(xi, data.frame(y=c(rep('control',n),rep('case',n))))
#' model <- bbl(data=dat, groups=c('control','case'))
#' 
#' cv <- crossval(object=model, method='mf', eps=seq(0.1,0.9,0.1))
#' plot(cv, type='b')
#' @export
crossval <- function(object, lambda=0.1, lambdah=0,
                     eps=0.9, nfold=5, method='pseudo', 
                     naive=FALSE, use.auc=TRUE, verbose=1, useC=TRUE, 
                     prior.count=TRUE, progress.bar=FALSE, 
                     fixL=FALSE, ...){
  
  groups <- object@groups
  Ly <- length(groups)
  nsample <- NROW(object@data)
  y <- object@data[,which(object@y==colnames(object@data))]
  
  if(Ly!=2) use.auc <- FALSE
  
  if(naive){
    method <- 'mf'
    eps <- 0
  }
  
  if(method=='pseudo'){ 
    reglist <- lambda
    if(length(lambdah)>1) 
      stop('Only a single value of lambdah allowed')
  }
  else if(method=='mf'){ 
    if(!naive) reglist <- eps
    else reglist <- 0
  }
  else stop('Unknown method')
  
  res <- NULL
  for(reg in reglist){
    if(is.null(lambdah)) regh <- reg
    else regh <- lambdah
    if(verbose>0){ 
      if(method=='pseudo'){
        if(regh>0)
          cat('Cross validation under lambda = (',reg,
            ',',regh,')\n',sep='')
        else
          cat('Cross validation under lambda = ',reg,'\n',sep='')
      }
      else cat('Cross validation under epsilon = ',reg,'\n',sep='')
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
      obval <- '['(x=object, i=ival, remove.const=FALSE)
      obtrain <- '['(x=object, i=itrain, remove.const=FALSE)
      if(method=='pseudo')
        obtrain <- train(object=obtrain, method=method, lambda=reg,  
                         naive=naive, verbose=verbose-1, lambdah=regh,
                         fixL=fixL, ...)
      else
        obtrain <- train(object=obtrain, method=method, eps=reg,
                         naive=naive, verbose=verbose-1, lambdah=lambdah,
                         fixL=fixL, ...)
      pr <- predict(object=obtrain, newdata=obval@data, logit=!use.auc,
                    useC=useC, progress.bar=progress.bar, 
                    verbose=verbose-1, naive=naive)
      pred <- rbind(pred, cbind(data.frame(y=y[ival], pr)))
    }
    if(use.auc){ auc <- pROC::roc(response=pred$y, levels=groups, 
                               predictor=pred[,3], direction='<')$auc
      if(verbose>0) cat(' AUC = ',auc,'\n',sep='')
    }
    else{
      yhat <- groups[apply(pred[,-1],1,which.max)]
      score <- mean(pred$y==yhat)
      if(verbose>0) cat(' Prediction score = ',score,'\n',sep='')
    }
    if(method=='pseudo'){
      if(use.auc) rx <- data.frame(lambda=reg, auc=auc)
      else rx <- data.frame(lambda=reg, score=score)
    } else{
      if(use.auc) rx <- data.frame(epsilon=reg, auc=auc)
      else rx <- data.frame(epsilon=reg, score=score)
    }
    res <- rbind(res, rx)
  }
  
  return(res)
}