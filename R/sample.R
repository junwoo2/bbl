eh <- function(si, h, J, numeric=FALSE){

  N <- length(si)
  e <- 0

  for(i in seq_len(N)){
    if(si[i]==0) next
    if(numeric) e <- e + h[[i]][1]*si[i]
    else e <- e + h[[i]][si[i]]
    if(i < N){
      for(j in seq(i+1,N)){
        if(si[j]==0) next
        if(numeric) e <- e + J[[i]][[j]][1]*si[i]*si[j]
        else e <- e + J[[i]][[j]][si[i],si[j]]
      }
    }
  }
  return(exp(e))
}

enum <- function(si, L, i, h, J, e=NULL, numeric=FALSE){

  for(s in seq(0,L[i]-1)){
    si[i] <- s
    if(i>1) e <- enum(si, L, i-1, h, J, e, numeric=numeric)
    else e <- rbind(e, cbind(t(si),eh(si, h, J, numeric=numeric)))
  }
  N <- length(si)
  return(e)
}

#' Generate random samples of observed/hidden states
#' @param nsample Sample size
#' @param L Number of hidden states per site 
#' @param n Total number of sites
#' @param nrepl Number of replicates of chains to stitch
#' @param h field parameter
#' @param J Coupling parameters. List of length equal to maximum distance. 
#'          d'th element is a square matrix of dimension d.
#' @param numeric Return numeric code of factors minus 1 (0,...,L-1)
#' @export
sample_xi <- function(nsample=1, numeric=FALSE, L=NULL, predictors=NULL, 
                      h, J, code_out=FALSE){

  if(numeric){ 
    if(is.null(L)) stop('L must be given for numeric model')
    nvar <- length(L)
  } else{
    L <- NULL
    for(p in predictors) L <- c(L, length(p))
    if(length(predictors)!=length(h) | length(predictors)!=length(J))
    stop("Predictors and h,J sizes don't match")
    nvar <- length(predictors)
  }

  e <- enum(si=rep(0,nvar), L=L, i=nvar, h, J, numeric=numeric)
  sid <- sample(NROW(e), size=nsample, replace=TRUE, prob=e[,nvar+1])
  si <- e[sid,seq_len(nvar)]
  
  if(!code_out & !numeric){
    for(i in seq_len(nvar)){
      x <- data.frame(x=factor(predictors[[i]][si[,i]+1], 
                                          levels=predictors[[i]]))
      if(i==1) fsi <- x
      else fsi <- cbind(fsi,x)
    }
  } else fsi <- as.data.frame(si)
  rownames(fsi) <- seq_len(nsample)
  colnames(fsi) <- seq_len(nvar)
  return(fsi)
}

#' @export
randompar <- function(predictors, h0=0, dh=1, J0=0, dJ=1, distr='unif'){
  
  m <- length(predictors)
  L <- NULL
  for(p in predictors) L <- c(L, length(p))
  
  h <- J <- vector('list',m)
  for(i in seq_len(m)) 
    J[[i]] <- vector('list',m)
  
  for(i in seq_len(m)){
    if(distr=='unif')
      h[[i]] <- runif(n=L[i]-1,min=h0-dh,max=h0+dh)
    else
      h[[i]] <- rnorm(n=L[i]-1,mean=h0, sd=dh)
    names(h[[i]]) <- predictors[[i]][-1]
    for(j in seq(i,m)){
      if(i==j) x <- 0
      else{ 
        if(distr=='unif') x <- runif(n=(L[i]-1)*(L[j]-1), min=J0-dJ, max=J0+dJ)
        else x <- rnorm(n=(L[i]-1)*(L[j]-1), mean=J0, sd=dJ)
      }
      x <- matrix(x, nrow=L[i]-1, ncol=L[j]-1)
      rownames(x) <- predictors[[i]][-1]
      colnames(x) <- predictors[[j]][-1]
      J[[i]][[j]] <- x
      if(i!=j) J[[j]][[i]] <- t(x)
    }
  }
  
  return(list(h=h, J=J))
}