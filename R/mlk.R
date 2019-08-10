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
                       tolerance=1e-5, verbose=1, naive=FALSE){
  
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
    
  } else if(method=='mft'){
    theta <- mle_mft(xi=xi, L=L, dmax=dmax, eps=eps)
  } else stop('unknown method in mlestimate')

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
mle_mft <- function(xi, L, dmax=1, eps=1){

  nsite <- ncol(xi)
  nsample <- nrow(xi)
  csi <- matrix(0, nrow=nsample, ncol=nsite*(L-1))
  for(k in 1:nsample){ 
    for(i in 1:N){
      if(xi[k,i]==1) next
      j <- (i-1)*(L-1)
      csi[k, j+xi[k,i]-1] <- 1
    }
  }
  cij <- cov(csi)
  cijb <- eps*cij + (1-eps)*mean(diag(cij))*diag(1,nrow=nsite*(L-1))
  mi <- matrix(colMeans(csi), nrow=nsite, ncol=L-1, byrow=TRUE)
  m0 <- 1-rowSums(mi)
  Jij <- solve(a=cijb, b=diag(-1, nrow=nsite*(L-1)))
  diag(Jij) <- 0
  hi <- log(mi/m0)
  for(i in 1:N) for(j in 1:N){
    if(i==j) next
    for(l in 1:(L-1)) for(q in 1:(L-1))
      hi[i,l] <- hi[i,l] - Jij[(i-1)*(L-1)+l,(j-1)*(L-1)+q]*mi[j,q]
  }
  h <- colMeans(hi)
  J <- vector('list',dmax)
  for(d in 1:dmax) J[[d]] <- matrix(0,nrow=L-1,ncol=L-1)
  s <- rep(0, dmax)  # counter
  for(i in 1:(N-1)) for(j in (i+1):N){
    d <- j-i
    if(d > dmax) next
    s[d] <- s[d] + 1
    for(l in 1:(L-1)) for(q in 1:(L-1))
      J[[d]][l,q] <- J[[d]][l,q] + Jij[(i-1)*(L-1)+l,(j-1)*(L-1)+q]
  }
  for(d in 1:dmax) J[[d]] <- J[[d]]/s[d]
  
  return(list(h=h, J=J))
}
