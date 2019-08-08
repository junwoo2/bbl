#' Predict class from new data using bbm model
#' 
#' @param newdata List of names \code{xi} and \code{y}; new data for 
#'                which prediction is made
#' @param logit Return predictors whose logistic function gives probability;
#'              otherwise return probability itself
#' @return Matrix of predictors/posterior proabilities
#' @export
setMethod('predict', 'bbm', function(object, newdata=NULL, logit=TRUE){
  
# browser()
  if(is.null(newdata)) data <- object@data # self-prediction
  else data <- newdata
  iy <- which(colnames(data)=='y')
  xi <- data[,-iy]
  y <- data$y
  
  Ly <- length(object@groups)
  L <- length(object@predictors)
  h <- object@h
  J <- object@J
  nsample <- NROW(xi)
  lz <- py <- rep(0, Ly)
  
  xid <- matrix(0, nrow=nsample, ncol=m)
  for(i in seq_len(m))
    xid[,i] <- match(xi[,i],object@predictors) - 1
  
  for(iy in seq_len(Ly)){
    lz[iy] <- Zeff(L, m, h[[iy]], J[[iy]], xid)  # log partition function
    py[iy] <- sum(y==object@groups[iy])        # marginal distribution P(y)
  }
  py <- py/nsample
  
  ay <- matrix(0, nrow=nsample, ncol=Ly)
  for(k in seq_len(nsample)){
    x <- xid[k,]
    for(iy in seq_len(Ly-1)){
      w0 <- ham(x, h[[iy]], J[[iy]]) - lz[iy] + log(py[iy])
      e <- 0
      for(iyp in seq_len(Ly)){
        if(iyp==iy) next()
        w <- ham(x, h[[iyp]], J[[iyp]]) - lz[iyp] + log(py[iyp])
        e <- e + exp(w - w0)
      }
      ay[k, iy] <- -log(e)
    }
    pry <- 1-sum(1/(1+exp(-ay[k,-Ly])))  # sum rule for sum_y {Pr(y|x)}=1
    ay[k, Ly] <- log(pry/(1-pry))
  }
  
  if(!logit) ay <- 1/(1+exp(-ay))  # posterior probability
  rownames(ay) <- seq_len(nsample)
  colnames(ay) <- object@groups
  return(ay)
})

ham <- function(x, h, J){

  m <- NROW(h)
  
  e <- 0
  for(i in seq_len(m)){
    if(x[i]==0) next()
    e <- e + h[i,x[i]]
    for(j in seq_len(m)){
      if(j==i | x[j]==0) next()
      e <- e + J[[i]][[j]][x[i],x[j]]/2
    }
  }
  return(e)  
}

# Pseudo partition function using data xi
Zeff <- function(L, m, h, J, xi){
  
  nsample <- NCOL(xi)
  lz <- 0
  for(k in seq_len(nsample)) for(i in seq_len(m)){
    z <- 1
    for(a in seq(1, L-1)){
      e <- h[i,a]
      for(j in seq(1,m)){
        if(j==i) next()
        xj <- xi[k,j]
        if(xj==0) next()
        e <- e + J[[i]][[j]][a,xj]
      }
      z <- z + exp(e)
    }
    lz <- lz + log(z)
  }
  lz <- lz/(nsample*m)
  
  return(lz)
}

#' @export
setMethod('[', 'bbm', function(x,i,j){
  
  if(!missing(i)){
    i <- as.vector(i)
    x@data <- x@data[i,]
  }
  return(x)
})
