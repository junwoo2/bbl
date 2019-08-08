eh <- function(si, L, h, J){

  N <- length(si)
  e <- 0

  for(i in 1:N){
    if(si[i]==0) next
    e <- e + h[i,si[i]]
    if(i < N){
      for(j in (i+1):N){
        if(si[j]==0) next
        e <- e + J[[i]][[j]][si[i],si[j]]
      }
    }
  }
  return(exp(e))
}

enum <- function(si, L, i, h, J, e=NULL){

  for(s in seq(0,L-1)){
    si[i] <- s
    if(i>1) e <- enum(si, L, i-1, h, J, e)
    else e <- rbind(e, cbind(t(si),eh(si, L, h, J)))
  }
  N <- length(si)
  if(nrow(e)==L^N) e[,N+1] <- e[,N+1]/sum(e[,N+1])
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
sample_si <- function(nsample=1, predictors=c('0','1'), 
                      nrepl=1, nsite, h, J, mc=FALSE,
                      nstep=1000, progress.bar=TRUE, numeric=FALSE){

  if(!is.character(predictors)) 
    predictors <- as.character(predictors)
  predictors <- unique(predictors)
  L <- length(predictors)
  
  if(!all(predictors[-1]==colnames(h))) 
    stop("Predictors and h column names don't match")
  
  if(!mc){
    if(NROW(h)!=nsite | NCOL(h)!=L-1 | length(J)!=nsite |
       !all(dim(J[[1]][[2]])==c(L-1,L-1)))
         stop('Incorrect dimension of parameters')
    e <- enum(si=rep(0,nsite), L=L, i=nsite, h, J)
    si <- NULL
    for(k in 1:nsample){
      sid <- sample(1:nrow(e), size=nrepl, replace=TRUE, 
                    prob=e[,nsite+1])
      si0 <- e[sid,1:nsite]
      si0 <- c(t(si0))
      si <- rbind(si, si0)
    }
  }
  else{
    mc <- mc.sample(si=rep(1,nsite),L=L,h,J,nstep=nstep, 
                    progress.bar=progress.bar)
    sid <- sample(nstep,size=nsample, replace=FALSE)
    si <- mc[sid,]
  }
  
  if(!numeric)
    si <- as.data.frame(matrix(predictors[si+1], nrow=nsample, ncol=nsite))
  rownames(si) <- seq_len(nsample)
  colnames(si) <- seq_len(nsite)
  return(si)
}

mc.sample <- function(si, L=L, h, J, nstep=1000, progress.bar=TRUE){
  
  N <- length(si)
  dmax <- length(J)
  sa <- NULL
  pb <- txtProgressBar(style=3)
  for(istep in seq_len(nstep)){
    for(i in seq_len(N)){
      sii <- si[i]
      if(sii==1) e0 <- 0
      else e0 <- energy(i, si, h, J)
      sp <- seq_len(L)
      sp <- sp[sp!=sii]
      s2 <- sample(sp, size=1)
      si2 <- si
      si2[i] <- s2
      if(s2==1 ) e1 <- 0
      else e1 <- energy(i, si2, h, J)
      de <- e1-e0
      move <- TRUE
      if(de>0) if(runif(n=1)<exp(-de)) move <- FALSE
      if(move) si <- si2
    }
    sa <- rbind(sa,si)
    setTxtProgressBar(pb,value=istep/nstep)
  }
  close(pb)
  return(sa)
}

energy <- function(i, si, h, J){
  
  N <- length(si)
  e <- h[si[i]-1]
  for(j in seq_len(N)){
    if(si[j]>1) 
      e <- e + J[i,j,2*(si[i]-1)+si[j]-1]
  }
  return(e)
}

#' @export
GenRandomPar <- function(m, predictors=NULL, L=NULL, h0=0, dh=1, J0=0, dJ=1, 
                         distr='unif'){
  
  if(!is.null(predictors)){
    predictors <- unique(predictors)
    L <- length(predictors)
  } else if(is.null(L)) stop('Either predictors or L must be given')
  
  if(distr=='unif')
    h <- matrix(runif(n=m*(L-1),min=h0-dh,max=h0+dh), nrow=m, ncol=L-1)
  else
    h <- matrix(rnorm(n=m*(L-1),mean=h0, sd=dh), nrow=m, ncol=L-1)
  if(!is.null(predictors))
    colnames(h) <- predictors[-1]
  J <- vector('list',m)
  for(i in seq_len(m)) J[[i]] <- vector('list',m)

  for(i in seq_len(m)) for(j in seq(i,m)){
    if(i==j) x <- 0
    else{ 
      if(distr=='unif') x <- runif(n=(L-1)^2, min=J0-dJ, max=J0+dJ)
      else x <- rnorm(n=(L-1)^2, mean=J0, sd=dJ)
    }
    x <- matrix(x, nrow=L-1, ncol=L-1)
    if(!is.null(predictors)) rownames(x) <- colnames(x) <- predictors[-1]
    J[[i]][[j]] <- x
    if(i!=j) J[[j]][[i]] <- t(x)
  }
  
  return(list(h=h, J=J))
}