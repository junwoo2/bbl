#' Maximum likelihood estimate
#' 
#' Perform inference of (h,J) for a given discrete target y
#' 
#' Given an input hidden state configuration, either pseudo-likelihood
#' of mean-field theory is used to find the maximum likelihood estimate
#' of single-site \code{h} and coupling \code{J} parameters.
#' 
#' @param xi Hidden state configuration. A matrix of dimension
#'           n by m, where n is the number of samples and m is 
#'           the number of sites. Each element ranges in value from
#'           1 to \code{L}.
#' @param L Number of possible levels for each hidden variable.
#' @param dmax Maximum distance between sites with non-zero interaction.
#' @param method Method to be used for maximum likelihood estimation. Either
#'            pseudo-likelihood (\code{'pseudo'}) or mean-field theory 
#'            (\code{'mft'}).
#' @param lambda Penalizer for interaction \code{J} with \code{'pseudo'}.
#' @param Lh  Penalizer for single-site parameter \code{h} with
#'            \code{'pseudo'}.
#' @param eps Regularization parameter for \code{'mft'}; allowed value
#'            is between 0 (no interaction) and 1 (full interaction).
#' @return List of inferred parameters \code{h} and \code{J}.
#' @export

mlestimate <- function(xi, L=2, method='pseudo', numeric=FALSE, lambda=0, 
                       symmetrize=TRUE, eps=1, nprint=100, itmax=10000,
                       tolerance=1e-5, verbose=1, naive=FALSE,
                       nullcount=1){
  
  if(method=='pseudo'){
    Lambda <- c(lambda)
    Nprint <- c(nprint)
    Itmax <- c(itmax)
    Tol <- c(tolerance)
    Verbose <- c(verbose)
    if(!is.numeric(xi)) stop('Input data to mlestimate must be numeric')
    if(!naive)
      theta <- pseudo_mle(xi, L, numeric, Lambda, Nprint, Itmax, Tol, Verbose)
    else
      theta <- naive_bayes(xi=xi, L=L, numeric=numeric, verbose=verbose)
    
    m <- NCOL(xi)
    if(numeric) Lp <- 1
    else Lp <- L-1
    h <- matrix(0, nrow=m, ncol=Lp)
    for(i in seq_len(m)) h[i,] <- theta$h[[i]]
    J <- vector('list',m)
    for(i in seq_len(m)) J[[i]] <- vector('list',m)
    for(i in seq(1,m)) for(j in seq(i,m)){
      x <- matrix(theta$J[[i]][[j]], nrow=Lp, ncol=Lp, byrow=TRUE)
      xt <- matrix(theta$J[[j]][[i]], nrow=Lp, ncol=Lp, byrow=TRUE)
      if(i<j & symmetrize){ 
        x <- (x + t(xt))/2
        xt <- t(x)
      }
      J[[i]][[j]] <- x
      J[[j]][[i]] <- xt
    }
    return(list(h=h, J=J, mle=theta$lkl, lz=theta$lz))
  }
  else if(method=='mf'){
    theta <- meanfield(xi=xi, L=L, eps=eps, numeric=numeric,
                       nullcount=nullcount)
    return(list(h=theta$h, J=theta$J))
  } else stop('unknown method in mlestimate')

}

# Naive Bayes (no interaction)
naive_bayes <- function(xi, L, numeric=FALSE, verbose){
  
  nsample <- NROW(xi)
  m <- NCOL(xi)
  
  if(numeric) Lp <- 1
  else Lp <- L-1
  
  h <- matrix(0, nrow=m, ncol=Lp)
  
  f0 <- colMeans(xi==0)
  for(l in seq_len(L-1)){
    if(numeric)
      h[,1] <- log(colMeans(xi==l)/f0)/l
    else
      h[,l] <- log(colMeans(xi==l)/f0)
  }
  
  J <- vector('list',m)
  for(i in seq_len(m)){ 
    J[[i]] <- vector('list',m)
    for(j in seq_len(m)) J[[i]][[j]] <- 0
  }
  
  return(list(h=h, J=J, lkl=NA))
}

#' @export
meanfield <- function(xi, L, eps=1, numeric=FALSE, nullcount=0){

  nsite <- ncol(xi)
  nsample <- nrow(xi)
  if(numeric) Lp <- 1
  else Lp <- L-1
  
  csi <- matrix(0, nrow=nsample, ncol=nsite*Lp)
  for(k in seq_len(nsample)) for(i in seq_len(nsite)){
    if(xi[k,i]==0) next()
    if(numeric)
      csi[k, (i-1)*Lp+1] <- xi[k,i]
    else
      csi[k, (i-1)*Lp+xi[k,i]] <- 1
  }
  if(sum(colSums(csi)==0)>0) nullcount <- 1
  mi <- (1/(Lp+1) + colSums(csi))/(nsample+nullcount)
  mi <- matrix(mi, nrow=nsite, ncol=Lp, byrow=TRUE)
# cij <- cov(csi)
  cij <- (t(csi) %*% csi + 1/(Lp+1)^2) / (nsample-1+nullcount)
  cijb <- eps*cij + (1-eps)*mean(diag(cij))*diag(1,nrow=nsite*Lp)

  if(!numeric){
    m0 <- 1-rowSums(mi)
    h <- log(mi/m0)
  } else
    h <- invav(mi,L)

  Jij <- solve(a=cijb, b=diag(-1, nrow=nsite*Lp))
  diag(Jij) <- 0
  for(i in seq_len(nsite)) for(j in seq_len(nsite)){
    if(i==j) next
    for(l in seq_len(Lp)) for(q in seq_len(Lp))
      h[i,l] <- h[i,l] - Jij[(i-1)*Lp+l,(j-1)*Lp+q]*mi[j,q]
  }
  J <- vector('list',m)
  for(i in seq_len(m)) J[[i]] <- vector('list',m)
  for(i in seq(1,nsite)) for(j in seq(i,nsite)){
    if(i==j){ 
      J[[i]][[i]] <- matrix(0,nrow=Lp,ncol=Lp)
      next()
    }
    w <- matrix(0, nrow=Lp, ncol=Lp)
    for(l in seq_len(Lp)) for(q in seq_len(Lp))
      w[l,q] <- Jij[(i-1)*Lp+l,(j-1)*Lp+q]
    J[[i]][[j]] <- w
    J[[j]][[i]] <- t(w)
  }
  
  return(list(h=h, J=J))
}

# Inverse function of mean single-site average mi

invav <- function(mi, L){
  
  f <- function(h, mi){
    up <- down <- 0
    eh <- exp(h)
    for(x in seq(0,L-1)){
      ehx <- eh^x
      up <- up + x*ehx
      down <- down + ehx
    }
    m-up/down
  }
  
  z <- NULL
  for(m in mi){
    if(m==0) hi <- -Inf
    else if(m==L-1) hi <- Inf
    else hi <- uniroot(f, interval=c(-10,10), mi=m)$root
    z <- c(z, hi)
  }
  
  z <- matrix(z, nrow=length(mi),ncol=1)
  return(z)
}