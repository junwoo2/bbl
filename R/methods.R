#' Print Boltzmann Bayes Learning Fits
#' 
#' This method displays model structure and first elements of coefficients
#' 
#' Displays the call to \code{\link{bbl}}, response variable and its levels,
#' predictors and their levels, and the first few fit coefficients.
#' 
#' @param x An object of class \code{bbl}, usually dervied from a call to
#'          \code{\link{bbl}}.
#' @param showcoeff Display first few fit coefficients
#' @param maxcoeff Maximum number of coefficients to display
#' @param ... Further arguments passed to or from other methods
#' @import methods
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
        cat('dh_[',names(predictors)[i],']^(',x$groups[iy],'): \n',sep='')
        print(hj$h[[iy]][[i]]-hj$h0[[i]])
        cat('\n')
      }
    }
    if(maxi < nvar) cat('...\n') else cat('\n')
    for(i in seq_len(maxi-1)) for(j in seq(from=i+1,to=maxi)){
      if(!(x$qJ)[i,j]) next()
      for(iy in seq_len(Ly)){
        cat('dJ_[',names(predictors)[i],',',
            names(predictors)[j],']^(',x$groups[iy],'): \n',sep='')
        print(hj$J[[iy]][[i]][[j]]-hj$J0[[i]][[j]])
        cat('\n')
      }
    }
    if(maxi < nvar) cat('...\n') else cat('\n')
  }
}

#' Naive Bayes Summary
#' 
#' Estimate significant of predictor-group association using naive Bayes model
#' 
#' This \code{summary.bbl} method gives a rough overview of associations
#' within a \code{bbl} fit object via naive Bayes coefficients and test 
#' p-values. Note that naive Bayes results displayed ignore interactions 
#' even when interactions are present in the model being displayed. This
#' feature is because simple analytic results exist for naive Bayes 
#' coefficients and test p-values. The likelihood ratio test is with respect
#' to the null hypothesis that coefficients are identical for all response
#' groups.   
#' 
#' @param object Object of class \code{bbl}
#' @param prior.count Prior count to be used for computing coefficients
#'        and test results. If \code{0}, will produce \code{NA}s for 
#'        factor levels without data points.
#' @param ... Other arguments to methods.
#' @return Object of class \code{summary.bbl} extending \code{bbl} class; 
#'      a list with extra components
#'    \item{h}{List of bias coefficients of response groups under naive 
#'    Bayes approximation}
#'    \item{h0}{Bias coefficients of pooled group under naive Bayes}
#'    \item{chisqNaive}{Vector of chi-square statistics for likelihood ratio test
#'    for each predictor}
#'    \item{dfNaive}{Vector of degrees of freedom for likelihood ratio test for
#'    each predictor}
#'    \item{pvNaive}{Vector p-values for each predictor}
#' @export
summary.bbl <- function(object, prior.count=1, ...){
  
  naive <- sum(object$qJ)==0  # no interaction
  xlevels <- object$xlevels
  data <- object$model
  y <- data[,object$groupname]
  x <- data[,names(xlevels)]
  if(is.null(object$freq))
    freq <- rep(1L, length(y))
  else
    freq <- object$freq
  
  dev <- rep(0,length(xlevels))
  names(dev) <- names(xlevels)
  h <- vector('list',length(object$groups))
  names(h) <- object$groups
  h0 <- vector('list', length(xlevels))
  names(h0) <- names(xlevels)
  for(iy in object$groups){
    h[[iy]] <- vector('list',length(xlevels))
    names(h[[iy]]) <- names(xlevels)
  }
    
  ntot <- sum(freq)
  for(i in seq_along(xlevels)){
    Li <- length(xlevels[[i]])
    f0 <- rep(0, Li)  # pooled inference
    names(f0) <- xlevels[[i]]
    for(w in xlevels[[i]])
      f0[as.character(w)] <- sum((data[,names(xlevels)[i]]==w)*freq)
    f0 <- f0 + prior.count/Li
    f0 <- f0 /(ntot + prior.count)
    h0[[i]] <- log(f0[-1]/f0[1])
    
    L <- 0
    for(iy in object$groups){
      ny <- sum(freq[y==iy])
      fv <- rep(0,Li)
      names(fv) <- xlevels[[i]]
      for(w in xlevels[[i]])
        fv[as.character(w)] <- 
        sum((data[y==iy,names(xlevels)[i]]==w)*freq[y==iy])
      fv <- fv + prior.count/Li
      fv <- fv /(ny+prior.count)
      h[[iy]][[i]] <- log(fv[-1]/fv[1])
      names(h[[iy]][[i]]) <- xlevels[[i]][-1]
      L <- L+ ny*sum(fv*log(fv))
    }
    dev[i] <- 2*(L - ntot*sum(f0*log(f0)))
  }

  df <- (length(object$groups)-1)*(lengths(xlevels)-1)
  pv <- pchisq(dev, df=df, lower.tail=F)
  
  ans <- c(object, list(h=h, h0=h0, chisqNaive=dev, dfNaive=df, pvNaive=pv))
  class(ans) <- 'summary.bbl'
  return(ans)
}

#' Print Summary of Boltzmann Bayes Learning
#' 
#' This method prints the summary of \code{bbl} object
#' 
#' The naive Bayes summary of \code{summary.bbl} object is displayed.
#' 
#' @param x Object of class \code{summary.bbl}
#' @param ... Other arguments to methods
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
  cat('Fit method: ',x$method,'\n',sep='')
  
  
  cat('\nnaive Bayes coefficients:\n')
  Ly <- length(x$groups)
  
  for(i in seq_len(nvar)){
    cat('h_',names(predictors)[i],': \n',sep='')
    H <- t(x$h[[1]][[i]])
    if(Ly>1){
      for(iy in seq(2,Ly))
        H <- rbind(H, x$h[[iy]][[i]])
      H <- rbind(H, x$h0[[i]])
    }
    rownames(H) <- c(x$groups,'pooled')
    print(H)
    cat('chisq = ',x$chisqNaive[i], ', df = ',x$dfNaive[i], ', Pr(>chisq) = ',
        x$pvNaive[i],'\n\n',sep='')
  }
}

#' Log likelihood for bbl object
#' 
#' Compute log likelihood from a fitted \code{bbl} object
#' 
#' This method will use inferred parameters from calls to \code{bbl} and data
#' compute the log likelihood.
#' 
#' @param object Object of class \code{bbl}
#' @param ... Other arguments to methods
#' @return Log likelihood value
#' @export
logLik.bbl <- function(object, ...){

  if(object$method!='mf') stop('Object was not trained with mf method')
  term <- object$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  predictors <- object$xlevels
  nvar <- length(predictors)
  xvar <- names(predictors)
  
  h <- coef(object)$h
  J <- coef(object)$J
  y <- object$model[,idy]
  xdat <- object$model[,names(predictors)]

  E <- 0.0
  for(yi in object$groups){
    dat <- xdat[y==yi,]
    ny <- NROW(dat)
    for(k in seq_len(ny)){
      for(i in seq_len(nvar)){
        xi <- xvar[[i]]
        zi <- as.character(dat[k,xi])
        if(!(zi %in% predictors[[xi]][-1])) next()
        E <- E + h[[yi]][[xi]][zi]
        if(i==nvar) next()
        for(j in seq(i+1,nvar)){
          xj <- xvar[[j]]
          zj <- as.character(dat[k,xj])
          if(zj %in% predictors[[xj]][-1])
            E <- E + J[[yi]][[xi]][[xj]][zi,zj]
        }
      }
    }
    E <- E - ny*object$lz[yi]
  }
  return(as.numeric(E))
}

#' Plot bbl object
#' 
#' Visualize bias and interaction parameters
#' 
#' This method displays a barplot of bias parameters and heatmaps
#' (one per response group) of interaction parameters. All parameters are
#' offset by the pooled values (single group inference) unless missing. 
#' @param x Object of class \code{bbl}
#' @param layout Matrix of layouts for arrangment of linear and interaction 
#'        parameters. If \code{NULL}, the top half will be used for linear 
#'        parameter barplot and bottom half will be divided into interaction 
#'        heatmaps for each response group.
#' @param hcol Color for linear barplots. Grayscale if \code{NULL}.
#' @param Jcol Color for interaction heatmaps. Default (\code{NULL}) is 
#'        \code{RdBu} from \code{RColorBrewer}.
#' @param npal Number of color scales.
#' @param ... Other graphical parameters for \code{\link{plot}}.
#'        
#' @export
plot.bbl <- function(x, layout=NULL, hcol=NULL, Jcol=NULL, npal=100, ...){
  
  oldpar <- par(no.readonly=TRUE)
  predictors <- x$xlevels
  naive <- sum(x$qJ)==0
  nvar <- length(predictors)
  Ly <- length(x$groups)
  h <- coef(x)$h
  h0 <- coef(x)$h0
  if(is.null(h0)) warning('Pooled parameters missing; use bbl(..., testNull=TRUE)')
  
  name <- NULL
  for(i in seq_len(nvar))
    name <- c(name, paste0(names(predictors)[i],':',predictors[[i]][-1]))
  
  if(!naive){
    if(is.null(layout))
      layout(matrix(c(rep(1,Ly), seq(2,Ly+1)), nrow=Ly, ncol=2, byrow=TRUE))
    else
      layout(layout)
  }
  
  H <- NULL
  for(i in seq_len(nvar)){
    if(!is.null(h0))
      hi <- t(h[[1]][[i]]-h0[[i]])
    else
      hi <- t(h[[1]][[i]])
    if(Ly>1){
      for(iy in seq(2,Ly)){
        if(!is.null(h0))
          hi <- rbind(hi, h[[iy]][[i]]-h0[[i]])
        else
          hi <- rbind(hi, h[[iy]][[i]])
      }
    }
    H <- cbind(H, hi)
  }
  rownames(H) <- x$groups
  bp <- barplot(H, beside=TRUE, col=hcol, names.arg=rep('',length(H)),main='', 
                ylab=expression(Delta * italic(h)), las=1, cex.axis=0.9, ...)

  axis(side=1, at=colMeans(bp), labels=name, las=2, cex.axis=0.9)
  if(is.null(hcol)) col <- gray.colors(Ly)
  legend(x='topright', fill=col, legend=x$groups, xpd=NA, 
         title=x$groupname, cex=0.5)
  
  if(naive){
    par(oldpar)
    return(invisible(x))
  }
  
  if(is.null(Jcol)) 
    Jcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,'RdBu'))(npal))
  ndim <- sqrt(length(unlist(coef(x)$J[[1]])))
  J0 <- coef(x)$J0
  IJ <- list()
  z0 <- NULL
  for(iy in seq_len(Ly)){
    J <- coef(x)$J[[iy]]
    I <- matrix(0, nrow=ndim, ncol=ndim)
    idx <- 1
    for(i in seq_len(nvar)) for(li in predictors[[i]][-1]){
      jdx <- 1
      for(j in seq_len(nvar)) for(lj in predictors[[j]][-1]){
        tmp <- J[[i]][[j]][li,lj]
        if(!is.null(J0))
          tmp <- tmp - J0[[i]][[j]][li,lj]
        I[idx, jdx] <- tmp
        jdx <- jdx + 1
      }
      idx <- idx + 1
    }
    z0 <- max(z0, max(abs(I)))
    IJ[[iy]] <- I
  }
  for(iy in seq_len(Ly)){
    image(IJ[[iy]], col=Jcol, zlim=c(-z0,z0), axes=FALSE)
    axis(side=1, at=seq(0,1,length.out=ndim),labels=name,las=2, lwd=0, cex.axis=0.8)
    axis(side=2, at=seq(0,1,length.out=ndim),labels=name,las=2, lwd=0, cex.axis=0.8)
    title(adj=0.5, main=bquote(
      Delta*italic(J)*'('*.(x$groupname)*'='*.(x$groups[iy])*')'),cex.main=1.0)
  }
  
  x0 <- 1.15
  y0 <- x0
  ny <- npal/10
  dy <- y0/40
  y <- seq(from=y0, to=y0-(ny-1)*dy, by=-dy)
  rect(xleft=x0,xright=x0*1.1, ytop=y, ybottom=y-dy, 
       col=rev(Jcol[seq(1,npal,length.out=ny)]),xpd=NA,lwd=0)
  zt <- signif(z0,digits=2)
  text(x=x0*1.3,y=max(y),adj=1, label=zt, cex=0.6, xpd=NA)
  text(x=x0*1.3,y=min(y),adj=1, label=-zt, cex=0.6, xpd=NA)
  
  par(oldpar)
  return(invisible(x))
}

#' Predict Response Group Using \code{bbl} Model
#' 
#' Make prediction of response group identity based on trained model
#' 
#' This method uses a new data set for predictors and trained \code{bbl} model
#' parameters to compute posterior probabilities of response group 
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
#' @param verbose Verbosity level
#' @param naive Nnaive Bayes. Skip all interaction terms.
#' @param progress.bar Display progress of response group probability. Useful
#'        for large samples.
#' @param ... Other arguments to methods
#' @return Data frame of predicted posterior probabilities with samples in rows
#'         and response groups in columns. The last column is the predicted
#'         response group with maximum probability.
#' @examples
#' set.seed(154)
#' 
#' m <- 5
#' L <- 3
#' n <- 1000
#' 
#' predictors <- list()
#' for(i in 1:m) predictors[[i]] <- seq(0,L-1)
#' names(predictors) <- paste0('v',1:m)
#' par <- list(randompar(predictors=predictors, dJ=0.5),
#'             randompar(predictors=predictors, h0=0.1, J0=0.1, dJ=0.5))
#' dat <- randomsamp(predictors=predictors, response=c('ctrl','case'), par=par, 
#'                  nsample=n)
#' dat <- dat[sample(n),]
#' dtrain <- dat[seq(n/2),]
#' dtest <- dat[seq(n/2+1,n),]
#' 
#' model <- bbl(y ~ .^2, data=dtrain)
#' pred <- predict(model, newdata=dtest)
#' score <- mean(dtest$y==pred$yhat)
#' score
#' 
#' auc <- pROC::roc(response=dtest$y, predictor=pred$case, direction='<')$auc
#' auc
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
predict.bbl <- function(object, newdata, logit=TRUE, verbose=1, naive=FALSE, 
                        progress.bar=FALSE, ...){
  
  term <- object$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  y <- object$model[,var[idy]]  # y is from training data
  predictors <- object$xlevels
  nvar <- length(predictors)
  if(verbose<=0) progress.bar <- FALSE
  
  if(missing(newdata)) data <- object$model
  else data <- newdata
  if(sum(!colnames(data)%in% var)!=0) 
    stop('Variable names in data not in object')
  x <- data[,var[-idy]]
  
  Ly <- length(object$groups)
  
  h <- coef(object)$h
  J <- coef(object)$J
  nsample <- NROW(data)
  xid <- matrix(0, nrow=nsample, ncol=nvar)
  
  for(i in seq_len(nvar)){
    if(sum(!levels(factor(x[,i])) %in% predictors[[i]])>0)
      stop('Levels in test data not in trained model')
    xid[,i] <- match(x[,i],predictors[[i]]) - 1
  }
  lz <- py <- rep(0, Ly)
  for(iy in seq_len(Ly)){
    if(length(h[[iy]])!=NCOL(xid) | length(J[[iy]])!=NCOL(xid))
      stop('Parameters and data sizes do not match')
    lz[iy] <- object$lz[iy]
    if(!is.null(object$freq))
      py[iy] <- sum((y==object$groups[iy])*object$freq)
    else
      py[iy] <- sum(y==object$groups[iy])
                  # marginal distribution P(y)
  }
  if(!is.null(object$freq)) py <- py/sum(object$freq)
  else py <- py/length(y)

  ay <- matrix(0, nrow=nsample, ncol=Ly)
  if(verbose>1) cat(' Predicting group probabilities...\n')
  if(progress.bar) pb <- txtProgressBar(style=3)
  for(k in seq_len(nsample)){
    xk <- xid[k,]
    E <- predict_class(xk, c(Ly), h, J, lz, py, c(naive))
    for(iy in seq_len(Ly)){
#     ay[k,iy] <- -log(sum(exp(E[-iy]-E[iy])))
      e <- E[-iy] - max(E[-iy])
      f <- log(sum(exp(e))) + max(E[-iy]) 
      ay[k,iy] <- -f + E[iy]
    }
    if(progress.bar) setTxtProgressBar(pb, k/nsample)
  }
  if(progress.bar) close(pb)
  
  if(!logit) ay <- 1/(1+exp(-ay))  # posterior probability
  rownames(ay) <- seq_len(nsample)
  colnames(ay) <- object$groups
  yhat <- factor(object$groups[apply(ay,1,which.max)],levels=object$groups)
  prob <- data.frame(as.data.frame(ay), yhat=yhat)
  
  return(prob)
}

#' Display Cross-validation Result
#' 
#' Print cross-validation optimal result and data frame
#' 
#' This method prints \code{\link{crossVal}} object with the optimal
#' regularization condition and maximum accuracy score on top and
#' the entire score profile as a data frame below.
#' @param x Object of class \code{cv.bbl}
#' @param ... Other arguments to methods
#' @export
print.cv.bbl <- function(x, ...){
  
  if(x$method=='mf') cat('Optimal epsilon = ',x$regstar,'\n',sep='')
  else cat('Optimal lambda = ',x$regstar,'\n',sep='')
  cat('Max. score: ',x$maxscore,'\n\n',sep='')
  print(x$cvframe)
  
}

#' Predict using Cross-validation Object
#' 
#' Use the optimal fitted model from cross-validation run to make prediction
#' 
#' This method will use the fitted model with maximum accuracy score returned
#' by a call to \code{\link{crossVal}} to make prediction on new data
#' 
#' @param object Object of class \code{cv.bbl}.
#' @param ... Other parameters to \code{\link{predict.bbl}}.
#' @return Data frame of prediction; see \code{\link{predict.bbl}}.
#' @export
predict.cv.bbl <- function(object, ...){
  
  class(object) <- 'bbl'
  predict(object, ...)
  
}

#' Plot Cross-validation Outcome
#' 
#' Plot cross-validation score as a function of regularization parameter
#' 
#' This function will plot accuracy score as a function of regularization parameter
#' from a call to \code{\link{crossVal}}.
#' 
#' @param x Object of class \code{cv.bbl} from a call to 
#' \code{\link{crossVal}}
#' @param type Symbol type in \code{\link{plot}}, present here to set default.
#' @param log Log scale argument to \code{\link{plot}}.
#' @param ... Other arguments to \code{\link{plot}}.
#' @export
plot.cv.bbl <- function(x, type='b', log='x', ...){
  
  plot(x$cvframe, type=type, log=log, ...)
  segments(x0=x$regstar,x1=x$regstar,y0=min(x$cvframe[,2]),y1=x$maxscore,
           lty=2, col='red')
  
}
