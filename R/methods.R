#' @export
print.bbl <- function(x, showcoeff=TRUE, maxcoeff=3L, ...){
  
  cat('\nCall:\n', paste0(deparse(x$call), sep='',collapse='\n'),'\n')
  
  term <- x$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  nvar <- length(x$xlevels)
  predictors <- x$xlevels
  cat(' ',nvar, ' predictor states:\n',sep='')
  for(i in seq_len(min(3,length(predictors))))
    cat('  ',names(x$xlevels)[i],'=',predictors[[i]],'\n',sep=' ')
  if(nvar > 3) cat('   ...\n',sep='')
  cat(' Responses:\n  ', var[1],'=',x$groups, '\n',sep=' ')
  
  if(showcoeff){
    cat('\nCoefficients:\n')
    Ly <- length(x$groups)
    hj <- x$coefficients
    maxi <- min(maxcoeff, nvar)
    for(i in seq_len(maxi)){
      for(iy in seq_len(Ly)){
        cat('h^(',x$groups[iy],')_[',names(predictors)[i],']: \n',sep='')
        print(hj$h[[iy]][[i]])
        cat('\n')
      }
    }
    if(maxi < nvar) print('...') else cat('\n')
    for(i in seq_len(maxi-1)) for(j in seq(from=i+1,to=maxi)){
      if(!(x$qJ)[i,j]) next()
      for(iy in seq_len(Ly)){
        cat('J^(',x$groups[iy],')_[',names(predictors)[i],',',
            names(predictors)[j],']: \n',sep='')
        print(hj$J[[iy]][[i]][[j]])
        cat('\n')
      }
    }
    if(maxi < nvar) print('...') else cat('\n')
  }
}

#' Naive Bayes Summary
#' 
#' Estimate significant of predictor-group association using naive Bayes model
#' 
#' @param object Object of class \code{bbl}.
#' @return Object of class \code{summary.bbl}. 
#' @export
summary.bbl <- function(object, ...){
  
  naive <- sum(object$qJ)==0  # no interaction
  xlevels <- object$xlevels
  data <- object$model
  y <- data[,object$groupname]
  x <- data[,names(xlevels)]
  
  dev <- rep(0,length(xlevels))
  names(dev) <- names(xlevels)
  h0 <- vector('list',length(object$groups))
  names(h0) <- object$groups
  for(iy in object$groups){
    h0[[iy]] <- vector('list',length(xlevels))
    names(h0[[iy]]) <- names(xlevels)
  }
    
  for(i in seq_along(xlevels)){
    L <- 0
    for(iy in object$groups){
      ny <- sum(y==iy)
      f <- table(data[y==iy,names(xlevels)[i]])
      f <- f/sum(f)
      h0[[iy]][[i]] <- log(f[-1]/f[1])
      names(h0[[iy]][[i]]) <- xlevels[[i]][-1]
      L <- L+ ny*sum(f*log(f))
    }
    f <- table(data[,names(xlevels)[i]])
    f <- f/sum(f)
    dev[i] <- L - NROW(data)*sum(f*log(f))
  }
  df <- lengths(xlevels)-1
  pv <- pchisq(dev, df=df, lower.tail=F)
  
  ans <- c(object, list(h0=h0, chisq=dev, df=df, pvalue=pv))
  class(ans) <- 'summary.bbl'
  return(ans)
}

#' @export
print.summary.bbl <- function(x, ...){
  
  cat('\nCall:\n', paste0(deparse(x$call), sep='',collapse='\n'),'\n')
  
  term <- x$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  nvar <- length(x$xlevels)
  predictors <- x$xlevels
  cat(' ',nvar, ' predictor states:\n',sep='')
  for(i in seq_len(min(3,length(predictors))))
    cat('  ',names(x$xlevels)[i],'=',predictors[[i]],'\n',sep=' ')
  if(nvar > 3) cat('   ...\n',sep='')
  cat(' Responses:\n  ', var[1],'=',x$groups, '\n',sep=' ')
  
  cat('\nNaive Bayes coefficients:\n')
  Ly <- length(x$groups)
  
  for(i in seq_len(nvar)){
    cat(names(predictors)[i],': \n',sep='')
    H <- x$h0[[1]][[i]]
    for(iy in seq(2,Ly))
      H <- rbind(H, x$h0[[iy]][[i]])
    rownames(H) <- x$groups
    print(H)
    cat('chisq = ',x$chisq[i], ', df = ',x$df[i], ', Pr(>chisq) = ',
        x$pvalue[i],'\n\n',sep='')
  }
}

#' Predict Response Group Using \code{bbl} Model
#' 
#' Make prediction of response group identity based on trained model
#' 
#' Will use new data set for predictors and trained \code{bbl} model
#' parameters and compute posterior probabilities of response group 
#' identity.
#' 
#' @param object Object of class \code{bbl} containing trained model
#' @param newdata Data frame of new data for which prediction is to
#'        be made. Columns must contain all of those in \code{model@data}.
#'        If column names are present, the columns will be matched 
#'        based on them. Extra columns will be ignored. If column names
#'        are not provided, the columns should exactly match 
#'        \code{model@data} predictor parts. If \code{NULL}, replaced
#'        by \code{model@data} (self-prediction).
#' @param logit Return predictors whose logistic function gives probability;
#'              otherwise return probability itself.
#' @param useC Use \code{C++} version for posterior probability computation.
#' @param verbose Verbosity level
#' @param naive Naive Bayes. Skip all interaction terms.
#' @param progress.bar Display progress of response group probability. Useful
#'        for large samples.
#' @return Matrix of predictors/posterior proabilities with samples in rows
#'         and response groups in columns.
#' @examples
#' set.seed(154)
#' 
#' m <- 5
#' L <- 3
#' n <- 1000
#' 
#' predictors <- list()
#' for(i in 1:m) predictors[[i]] <- seq(0,L-1)
#' par0 <- randompar(predictors=predictors, h0=0, J0=0, dJ=0.5)
#' xi0 <- sample_xi(nsample=n, predictors=predictors, h=par0$h, J=par0$J) 
#'
#' par1 <- randompar(predictors=predictors, h0=0.1, J0=0.1, dJ=0.5)
#' xi1 <- sample_xi(nsample=n, predictors=predictors, h=par1$h, J=par1$J) 
#'
#' xi <- rbind(xi0,xi1)
#' y <- c(rep(0,n),rep(1,n))
#' dat <- cbind(data.frame(y=y),xi)
#' dat <- dat[sample(2*n),]
#' dtrain <- dat[seq(n),]
#' dtest <- dat[seq(n+1,2*n),]
#' ytest <- dtest[,'y']
#' 
#' model <- bbl(data=dtrain)
#' model <- train(model)
#' 
#' pred <- predict(object=model, newdata=dtest)
#' yhat <- apply(pred,1,which.max)-1
#' score <- mean(ytest==yhat)
#' score
#' 
#' auc <- pROC::roc(response=ytest, predictor=pred[,2], direction='<')$auc
#' auc
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
predict.bbl <- function(object, newdata, logit=TRUE, useC=TRUE, 
                        verbose=1, naive=FALSE, progress.bar=FALSE){
  
  if(verbose<=0) progress.bar <- FALSE
  
  if(missing(newdata)) data <- object$model
  else data <- newdata
  
  term <- object$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  y <- data[,var[idy]]
  x <- data[,var[-idy]]
  
  Ly <- length(object$groups)
  
  h <- object$coefficients$h
  J <- object$coefficients$J
  nsample <- NROW(data)
  m <- length(object$xlevels)
  xid <- matrix(0, nrow=nsample, ncol=m)
  
  for(i in seq_len(m)){
    if(sum(!levels(factor(x[,i])) %in% object$xlevels[[i]])>0)
      stop('Levels in test data not in trained model')
    xid[,i] <- match(x[,i],object$xlevels[[i]]) - 1
  }
  lz <- py <- rep(0, Ly)
  for(iy in seq_len(Ly)){
    if(length(h[[iy]])!=NCOL(xid) | length(J[[iy]])!=NCOL(xid))
      stop('Parameters and data sizes do not match')
    lz[iy] <- object$lz[iy]
    py[iy] <- sum(y==object$groups[iy])        # marginal distribution P(y)
  }
  py <- py/nsample
  
  ay <- matrix(0, nrow=nsample, ncol=Ly)
  if(verbose>1) cat(' Predicting group probabilities...\n')
  if(progress.bar) pb <- txtProgressBar(style=3)
  for(k in seq_len(nsample)){
    xk <- xid[k,]
    if(!useC){
      E <- rep(0, Ly)
      for(iy in seq_len(Ly))
        E[iy] <- ham(xk, h[[iy]], J[[iy]], naive=naive) - 
          lz[iy] + log(py[iy])
    }else
      E <- predict_class(xk, c(Ly), h, J, lz, py, c(naive))
    for(iy in seq_len(Ly))
      ay[k,iy] <- -log(sum(exp(E[-iy]-E[iy])))
    if(progress.bar) setTxtProgressBar(pb, k/nsample)
  }
  if(progress.bar) close(pb)
  
  if(!logit) ay <- 1/(1+exp(-ay))  # posterior probability
  rownames(ay) <- seq_len(nsample)
  colnames(ay) <- object$groups
  yhat <- object$groups[apply(ay,1,which.max)]
  
  return(list(prob=ay, yhat=yhat))    
}