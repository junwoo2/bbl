#' Predict class from new data using bbm model
#' 
#' @param newdata List of names \code{xi} and \code{y}; new data for 
#'                which prediction is made
#' @param logit Return predictors whose logistic function gives probability;
#'              otherwise return probability itself
#' @return Matrix of predictors/posterior proabilities
#' @export
setMethod('predict', 'bbm', function(object, newdata=NULL, logit=TRUE,
                                     computeZ=FALSE, mf=FALSE, useC=TRUE){

  if(is.null(newdata)) data <- object@data # self-prediction
  else data <- newdata
  iy <- which(object@y==colnames(data))
  xi <- data[,-iy]
  y <- data[,iy]
  
  Ly <- length(object@groups)
  m <- NCOL(object@data) - 1
  L <- NULL
  for(i in seq_len(m)){
    if(object@type=='numeric') 
      L <- c(L, max(object@data[,i]))
    else
      L <- c(L, length(object@predictors[[i]]))
  }
  
  h <- object@h
  J <- object@J
  nsample <- NROW(xi)
  lz <- py <- rep(0, Ly)
  
  if(object@type=='numeric'){
    if(!is.numeric(xi[1,1])) stop('Numeric model requires numeric data')
    xid <- as.matrix(xi)
  }else{
    xid <- matrix(0, nrow=nsample, ncol=m)
    for(i in seq_len(m))
      xid[,i] <- match(xi[,i],object@predictors[[i]]) - 1
  }
  
  for(iy in seq_len(Ly)){
    if(computeZ)
      lz[iy] <- Zeff(L, m, h[[iy]], J[[iy]], xid, object@type=='numeric', mf=mf)  # log partition function
    else
      lz[iy] <- object@lz[iy]
    py[iy] <- sum(y==object@groups[iy])        # marginal distribution P(y)
  }
  py <- py/nsample
  
  if(!useC){
    ay <- matrix(0, nrow=nsample, ncol=Ly)
    for(k in seq_len(nsample)){
      x <- xid[k,]
      E <- rep(0, Ly)
      for(iy in seq_len(Ly))
        E[iy] <- ham(x, h[[iy]], J[[iy]], 
                     numeric=object@type=='numeric') - lz[iy] + log(py[iy])
      for(iy in seq_len(Ly))
        ay[k,iy] <- -log(sum(exp(E[-iy]-E[iy])))
    }
  }else
    ay <- predict_class(xid, c(Ly), h, J, c(object@type=='numeric'), lz, py)

  if(!logit) ay <- 1/(1+exp(-ay))  # posterior probability
  rownames(ay) <- seq_len(nsample)
  colnames(ay) <- object@groups
  return(ay)
})

ham <- function(x, h, J, numeric){

  m <- length(h)
  
  e <- 0
  for(i in seq_len(m)){
    if(x[i]==0) next()
    if(numeric) e <- e + h[[i]][1]*x[i]
    else e <- e + h[[i]][x[i]]
    for(j in seq_len(m)){
      if(j==i | x[j]==0) next()
      if(numeric) e <- e + J[[i]][[j]]*x[i]*x[j]/2
      else e <- e + J[[i]][[j]][x[i],x[j]]/2
    }
  }
  return(e)  
}

# Pseudo partition function using data xi
Zeff <- function(L, m, h, J, xi, numeric, mf=mf){
  
  nsample <- NROW(xi)
  if(numeric) Lp <- rep(1, m)
  else Lp <- L-1
  if(mf){
    fi <- vector('list',m)
    f0 <- rep(0, m)
    for(i in seq_len(m)){
      fi[[i]] <- rep(0, L[i]-1)
      for(l in seq_len(L[i]-1))
        fi[[i]][l] <- mean(xi[,i]==l)
      f0[i] <- 1 - sum(fi[[i]])
    }
    lz <- - sum(log(f0))
    for(i in seq(1,m-1)) for(j in seq(i+1,m))
      for(l0 in seq_len(Lp[i])) for(l1 in seq_len(Lp[j])){
        dz <- J[[i]][[j]][l0,l1]*fi[[i]][l0]*fi[[j]][l1]
        if(numeric) dz <- dz * l0*l1;
        lz <- lz - dz
      }
    return(lz)
  }
    
  lz <- 0
  for(k in seq_len(nsample)) for(i in seq_len(m)){
    z <- 1
    for(a in seq(1, L[i]-1)){
      if(numeric) e <- h[[i]][1]*a
      else e <- h[[i]][a]
      for(j in seq(1,m)){
        if(j==i) next()
        xj <- xi[k,j]
        if(xj==0) next()
        if(numeric)
          e <- e + J[[i]][[j]]*a*xj/2
        else
          e <- e + J[[i]][[j]][a,xj]/2
      }
      z <- z + exp(e)
    }
    lz <- lz + log(z)
  }
  lz <- lz/nsample
  
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
