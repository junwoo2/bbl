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
        cat('dh^(',x$groups[iy],')_[',names(predictors)[i],']: \n',sep='')
        print(hj$h[[iy]][[i]]-hj$h0[[i]])
        cat('\n')
      }
    }
    if(maxi < nvar) print('...') else cat('\n')
    for(i in seq_len(maxi-1)) for(j in seq(from=i+1,to=maxi)){
      if(!(x$qJ)[i,j]) next()
      for(iy in seq_len(Ly)){
        cat('dJ^(',x$groups[iy],')_[',names(predictors)[i],',',
            names(predictors)[j],']: \n',sep='')
        print(hj$J[[iy]][[i]][[j]]-hj$J0[[i]][[j]])
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
  for(iy in object$groups){
    h[[iy]] <- vector('list',length(xlevels))
    names(h[[iy]]) <- names(xlevels)
  }
    
  ntot <- sum(freq)
  for(i in seq_along(xlevels)){
    
    f0 <- rep(0, length(xlevels[[i]]))  # pooled inference
    names(f0) <- xlevels[[i]]
    for(w in xlevels[[i]])
      f0[w] <- sum((data[,names(xlevels)[i]]==w)*freq)
    f0 <- (f0 + prior.count)/(ntot + prior.count)
    h0 <- log(f0[-1]/f0[1])
    
    L <- 0
    for(iy in object$groups){
      ny <- sum(freq[y==iy])
      fv <- rep(0,length(xlevels[[i]]))
      names(fv) <- xlevels[[i]]
      for(w in xlevels[[i]])
        fv[w] <- sum((data[y==iy,names(xlevels)[i]]==w)*freq[y==iy])
#     f <- table(data[y==iy,names(xlevels)[i]])
      fv <- (fv + prior.count)/(ny+prior.count)
#     fv <- fv/sum(fv)
      h[[iy]][[i]] <- log(fv[-1]/fv[1]) - h0
      names(h[[iy]][[i]]) <- xlevels[[i]][-1]
      L <- L+ ny*sum(fv*log(fv))
    }
    dev[i] <- L - ntot*sum(f0*log(f0))
  }
  df <- lengths(xlevels)-1
  pv <- pchisq(dev, df=df, lower.tail=F)
  
  ans <- c(object, list(dhNaive=h, chisqNaive=dev, dfNaive=df, pvNaive=pv))
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
  cat('Fit method: ',x$method,'\n',sep='')
  
  
  cat('\nNaive Bayes coefficients:\n')
  Ly <- length(x$groups)
  
  for(i in seq_len(nvar)){
    cat('dH_',names(predictors)[i],': \n',sep='')
    H <- t(x$dhNaive[[1]][[i]])
    if(Ly>1){
      for(iy in seq(2,Ly))
        H <- rbind(H, x$dhNaive[[iy]][[i]])
    }
    rownames(H) <- x$groups
    print(H)
    cat('chisq = ',x$chisqNaive[i], ', df = ',x$dfNaive[i], ', Pr(>chisq) = ',
        x$pvNaive[i],'\n\n',sep='')
  }
}

#' Log likelihood for bbl object
#' @export
logLik.bbl <- function(x, ...){

  if(x$method!='mf') stop('Object was not trained with mf method')
  term <- x$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  predictors <- x$xlevels
  nvar <- length(predictors)
  xvar <- names(predictors)
  
  h <- coef(x)$h
  J <- coef(x)$J
  y <- x$model[,idy]
  xdat <- x$model[,names(predictors)]

  E <- 0.0
  for(yi in x$groups){
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
    E <- E - ny*x$lz[yi]
  }
  return(as.numeric(E))
}

#' Plot bbl object
#' @param x Object of class \code{bbl}
#' @param layout Matrix of layouts for arrangment of linear and interaction 
#'        parameters. If \code{NULL}, the top half will be used for linear parameter
#'        barplot and bottom half will be divided into interaction heatmaps
#'        for each response group.
#' @param hcol Color for linear barplots. Grayscale if \code{NULL}.
#' @param Jcol Color for interaction heatmaps. Default (\code{NULL}) is 
#'        \code{RdBu} from \code{RColorBrewer}.
#' @param npal Number of color scales.
#' @param mar Plot margins.
#' @param ... Other graphical parameters for \link{\code{plot}}.
#'        
#' @export
plot.bbl <- function(x, layout=NULL, hcol=NULL, Jcol=NULL, npal=100, 
                     mar=c(3,3,3,3), ...){
  
  oldpar <- par(no.readonly=TRUE)
  par(mar=mar)
  predictors <- x$xlevels
  nvar <- length(predictors)
  Ly <- length(x$groups)
  h <- coef(x)$h
  h0 <- coef(x)$h0
  
  name <- NULL
  for(i in seq_len(nvar))
    name <- c(name, paste0(names(predictors)[i],':',predictors[[i]][-1]))
  
  if(is.null(layout))
    layout(matrix(c(rep(1,Ly), seq(2,Ly+1)), nrow=Ly, ncol=2, byrow=TRUE))
  else
    layout(layout)
  
  H <- NULL
  for(i in seq_len(nvar)){
    hi <- t(h[[1]][[i]]-h0[[i]])
    if(Ly>1){
      for(iy in seq(2,Ly))
        hi <- rbind(hi, h[[iy]][[i]]-h0[[i]])
    }
    H <- cbind(H, hi)
  }
  rownames(H) <- x$groups
  bp <- barplot(H, beside=TRUE, col=hcol, names.arg=rep('',length(H)),main='', 
                ylab=expression(Delta * italic(h)), ...)

  axis(side=1, at=colMeans(bp), label=name, las=2)
  if(is.null(hcol)) col <- gray.colors(Ly)
  legend(x='topright', fill=col, legend=x$groups, xpd=NA, title=x$groupname,
         cex=0.7)
  
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
        I[idx, jdx] <- J[[i]][[j]][li,lj] - J0[[i]][[j]][li,lj]
        jdx <- jdx + 1
      }
      idx <- idx + 1
    }
    z0 <- max(z0, max(abs(I)))
    IJ[[iy]] <- I
  }
  for(iy in seq_len(Ly)){
    image(IJ[[iy]], col=Jcol, zlim=c(-z0,z0), axes=FALSE)
    axis(side=1, at=seq(0,1,length.out=ndim),label=name,las=2, lwd=0)
    axis(side=2, at=seq(0,1,length.out=ndim),label=name,las=2, lwd=0)
    title(adj=0.5, main=bquote(
      Delta*italic(J)*'('*.(x$groupname)*'='*.(x$groups[iy])*')'))
  }
  
  x0 <- 1.15
  y0 <- x0
  ny <- npal/10
  dy <- y0/40
  y <- seq(from=y0, to=y0-(ny-1)*dy, by=-dy)
  rect(xleft=x0,xright=x0*1.1, ytop=y, ybottom=y-dy, 
       col=rev(Jcol[seq(1,npal,length.out=ny)]),xpd=NA,lwd=0)
  zt <- signif(z0,digits=2)
  text(x=x0*1.2,y=max(y),adj=1, label=zt, cex=0.6, xpd=NA)
  text(x=x0*1.2,y=min(y),adj=1, label=-zt, cex=0.6, xpd=NA)
  
  par(oldpar)
  return(invisible(x))
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
  
  term <- object$terms
  idy <- attr(term,'response')
  var <- as.character(attr(term, 'variables')[-1])
  y <- object$model[,var[idy]]  # y is from training data
  predictors <- object$xlevels
  nvar <- length(predictors)
  if(verbose<=0) progress.bar <- FALSE
  
  if(missing(newdata)) data <- object$model
  else data <- newdata
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
  prob <- data.frame(as.data.frame(ay), yhat=yhat)
  
  return(prob)    
}

#' @export
print.cv.bbl <- function(x, ...){
  
  if(x$method=='mf') cat('Optimal epsilon = ',cv$regstar,'\n',sep='')
  else cat('Optimal lambda = ',cv$regstar,'\n',sep='')
  cat('Max. score: ',x$maxscore,'\n\n',sep='')
  print(x$cvframe)
  
}

#' @export
plot.cv.bbl <- function(x, type='b', log='x', ...){
  
  plot(x$cvframe, type=type, log=log, ...)
  segments(x0=x$regstar,x1=x$regstar,y0=min(x$cvframe[,2]),y1=x$maxscore,
           lty=2, col='red')
  
}