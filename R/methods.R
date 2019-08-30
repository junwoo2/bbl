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
setMethod('predict', 'bbl', 
          function(object, newdata=NULL, logit=TRUE, useC=TRUE, verbose=1, 
                   naive=FALSE, progress.bar=FALSE){

  if(verbose<=0) progress.bar <- FALSE
  iy <- which(object@y==colnames(object@data))
  y <- object@data[,iy]   # y is from training data
  xi0 <- object@data[,-iy]
  
  if(is.null(newdata)) 
    xi <- xi0           # self-prediction
  else{
    if(!is(newdata, 'data.frame')) stop('newdata must be a data frame')
    if(object@y %in% colnames(newdata))
      newdata <- newdata[,colnames(newdata)!=object@y]
    if(NCOL(newdata)!= NCOL(xi0)) 
      stop('Test data predictors do not match model')
    if(!is.null(colnames(newdata)))
      xi <- newdata[,match(colnames(newdata),colnames(xi0))]
  }
  
  Ly <- length(object@groups)
  m <- NCOL(object@data) - 1
  h <- object@h
  J <- object@J
  nsample <- NROW(xi)
  lz <- py <- rep(0, Ly)
  
  xid <- matrix(0, nrow=nsample, ncol=m)
  for(i in seq_len(m)){
    if(sum(!levels(factor(xi[,i])) %in% object@predictors[[i]])>0)
      stop('Levels in test data not in trained model')
    xid[,i] <- match(xi[,i],object@predictors[[i]]) - 1
  }
  
  for(iy in seq_len(Ly)){
    if(length(h[[iy]])!=NCOL(xi) | length(J[[iy]])!=NCOL(xi))
      stop('Parameters and data sizes do not match')
    lz[iy] <- object@lz[iy]
    py[iy] <- sum(y==object@groups[iy])        # marginal distribution P(y)
  }
  py <- py/nsample
  
  ay <- matrix(0, nrow=nsample, ncol=Ly)
  if(verbose>0) cat(' predicting group probabilities...\n')
  if(progress.bar) pb <- txtProgressBar(style=3)
    for(k in seq_len(nsample)){
    x <- xid[k,]
    if(!useC){
      E <- rep(0, Ly)
      for(iy in seq_len(Ly))
        E[iy] <- ham(x, h[[iy]], J[[iy]], naive=naive) - lz[iy] + log(py[iy])
    }else
      E <- predict_class(x, c(Ly), h, J, lz, py, c(naive))
    for(iy in seq_len(Ly))
      ay[k,iy] <- -log(sum(exp(E[-iy]-E[iy])))
    if(progress.bar) setTxtProgressBar(pb, k/nsample)
  }
  if(progress.bar) close(pb)

  if(!logit) ay <- 1/(1+exp(-ay))  # posterior probability
  rownames(ay) <- seq_len(nsample)
  colnames(ay) <- object@groups
  return(ay)
})

ham <- function(x, h, J, naive){

  m <- length(h)
  
  e <- 0
  for(i in seq_len(m)){
    if(x[i]==0) next()
    if(length(h[[i]])<x[i]) next()
    e <- e + h[[i]][x[i]]
    if(naive) next()
    for(j in seq_len(m)){
      if(j==i | x[j]==0) next()
      if(NROW(J[[i]][[j]])<x[i] | NCOL(J[[i]][[j]])<x[j]) next()
      e <- e + J[[i]][[j]][x[i],x[j]]/2
    }
  }
  return(e)  
}

#' @describeIn bbl
#' Subsetting of \code{bbl} object along rows (sample index)
#' @param x Object of class \code{bbl} to be subsetted
#' @param i Row index to keep
#' @param j Not used.
#' @param remove.const Remove predictor levels not found in data.
#' @export
setMethod('[', 'bbl', function(x,i,j, remove.const=TRUE){
  
  data <- x@data
  iy <- which(colnames(x@data)==x@y)
  if(!missing(i)){
    i <- as.vector(i)
    data <- data[i,]
  }
  if(!missing(j))
    stop('column subsetting not implemented')
  if(remove.const)
    x <- bbl(data=data, groups=x@groups, y=x@y)
  else
    x <- bbl(data=data, groups=x@groups, y=x@y, predictors=x@predictors)
  return(x)
})