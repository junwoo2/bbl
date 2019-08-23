#' Predict class from new data using bbl model
#' 
#' @param newdata List of names \code{xi} and \code{y}; new data for 
#'                which prediction is made
#' @param L Forces uniform predictor levels
#' @param logit Return predictors whose logistic function gives probability;
#'              otherwise return probability itself
#' @return Matrix of predictors/posterior proabilities
#' @export
setMethod('predict', 'bbl', function(object, newdata=NULL, logit=TRUE,
                                     L=NULL, useC=TRUE, verbose=1, naive=FALSE,
                                     progress.bar=FALSE){

  if(verbose<=0) progress.bar <- FALSE
# browser()
  iy <- which(object@y==colnames(object@data))
  y <- object@data[,iy]   # y is from training data
  xi0 <- object@data[,-iy]
  
  if(is.null(newdata)) 
    data <- xi0           # self-prediction
  else{
    if(object@y %in% colnames(newdata))
      newdata <- newdata[,colnames(newdata)!=object@y]
    if(NCOL(newdata)!= NCOL(xi0)) 
      stop('Test data predictors do not match model')
    if(!is.null(colnames(newdata)))
      xi <- newdata[,match(colnames(newdata),colnames(xi0))]
  }
  
  Ly <- length(object@groups)
  m <- NCOL(object@data) - 1
  if(is.null(L)){
    for(i in seq_len(m)){
      if(object@type=='numeric') 
        L <- c(L, max(object@data[,i])+1)
      else
        L <- c(L, length(object@predictors[[i]]))
    }
  } else L <- rep(L, m)
  
  h <- object@h
  J <- object@J
  nsample <- NROW(xi)
  lz <- py <- rep(0, Ly)
  
  if(object@type=='numeric'){
    if(!is.numeric(xi[1,1])) stop('Numeric model requires numeric data')
    xid <- as.matrix(xi)
  }else{
    xid <- matrix(0, nrow=nsample, ncol=m)
    for(i in seq_len(m)){
      if(sum(!levels(factor(xi[,i])) %in% object@predictors[[i]])>0)
        stop('Levels in test data not in trained model')
      xid[,i] <- match(xi[,i],object@predictors[[i]]) - 1
    }
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
        E[iy] <- ham(x, h[[iy]], J[[iy]], numeric=object@type=='numeric',
                     naive=naive) - lz[iy] + log(py[iy])
    }else
      E <- predict_class(x, c(Ly), h, J, c(object@type=='numeric'), lz, py,
                         c(naive))
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

ham <- function(x, h, J, numeric, naive){

  m <- length(h)
  
  e <- 0
  for(i in seq_len(m)){
    if(x[i]==0) next()
    if(numeric) e <- e + h[[i]][1]*x[i]
    else if(length(h[[i]])<x[i]) next()
    else e <- e + h[[i]][x[i]]
    if(naive) next()
    for(j in seq_len(m)){
      if(j==i | x[j]==0) next()
      if(numeric) e <- e + J[[i]][[j]]*x[i]*x[j]/2
      else if(NROW(J[[i]][[j]])<x[i] | 
              NCOL(J[[i]][[j]])<x[j]) next()
      else e <- e + J[[i]][[j]][x[i],x[j]]/2
    }
  }
  return(e)  
}

#' @export
setMethod('[', 'bbl', function(x,i,j){
  
  if(!missing(i)){
    i <- as.vector(i)
    x@data <- x@data[i,]
  }
  return(x)
})
