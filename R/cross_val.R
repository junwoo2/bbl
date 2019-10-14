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
#'        for \code{method = 'mf'}. Inference will be repeated
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
crossVal <- function(formula, data, freq=NULL,
                     lambda=1e-5, lambdah=0, eps=0.9, nfold=5, method='pseudo', 
                     naive=FALSE, use.auc=TRUE, verbose=1, useC=TRUE, 
                     prior.count=TRUE, progress.bar=FALSE, 
                     fixL=FALSE, ...){
  
  
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  term <- terms(formula, data=data)
  xlevels <- .getXlevels(term, m=data)
  m <- length(xlevels)
  idy <- attributes(term)$response
  resp <- all.vars(cl)[idy]
  y <- data[,resp]
  yx <- data[,c(resp,names(xlevels))]
  groups <- levels(factor(y))
  Ly <- length(groups)
  
  if(Ly==1) warning('Only one response group in data')
  if(Ly!=2) use.auc <- FALSE
  
  if(!is.null(freq)){
    if(!all.equal(freq, as.integer(freq))) stop('Non-integer freq')
    if(length(freq)!=NROW(data))
      stop('Length of freq does not match data')
    yx <- freq2raw(fdata=yx, freq=freq)
    y <- yx[,resp]
  }
  
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
  
  res <- NULL    # data frame of cv result
  bbopt <- NULL  # optimally trained bbl object
  maxscore <- -Inf
  regstar <- 0
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
        nval <- max(1,ceiling(ns/nfold))
        imax <- k*nval
        if(imax > ns) imax <- ns
        if(imax < 1) imax <- 1
        iyval <- idy[seq((k-1)*nval+1, imax)]
        iytrain <- idy[!idy %in% iyval]
        ival <- c(ival, iyval)
        itrain <- c(itrain, iytrain)
      }
      if(sum(is.na(ival)>0) | sum(is.na(itrain)>0)) stop('error in crossval')
      dval <- yx[ival, ,drop=FALSE]
      dtrain <- yx[itrain, ,drop=FALSE]
      if(method=='pseudo')
        obtrain <- bbl(formula, data=dtrain, method=method, lambda=reg, 
                       verbose=verbose-1, lambdah=regh, fixL=fixL, ...)
      else
        obtrain <- bbl(formula, data=dtrain, method=method, eps=reg, 
                       verbose=verbose-1, lambdah=lambdah, fixL=fixL, ...)
      pr <- predict(object=obtrain, newdata=dval, logit=!use.auc, useC=useC, 
                    progress.bar=progress.bar, verbose=verbose-1)
      pred <- rbind(pred, cbind(yx[ival,1,drop=FALSE], pr))
    }
    if(use.auc){
      score <- pROC::roc(response=pred[,1], levels=groups, 
                               predictor=pred[,3], direction='<')$auc
      if(verbose>0) cat(' AUC = ',score,'\n',sep='')
    }
    else{
      score <- mean(pred[,1]==pred$yhat)
      if(verbose>0) cat(' Prediction score = ',score,'\n',sep='')
    }
    if(method=='pseudo'){
      if(use.auc) rx <- data.frame(lambda=reg, auc=score)
      else rx <- data.frame(lambda=reg, score=score)
    } else{
      if(use.auc) rx <- data.frame(epsilon=reg, auc=score)
      else rx <- data.frame(epsilon=reg, score=score)
    }
    if(score > maxscore){
      maxscore <- score
      regstar <- reg
      bbopt <- obtrain
    }
    res <- rbind(res, rx)
  }
  
  cv <- c(bbopt, list(regstar=regstar, maxscore=maxscore, cvframe=res))
  class(cv) <- 'cv.bbl'
  
  return(cv)
}