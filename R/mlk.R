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

mlestimate <- function(xi, L=2, method='pseudo', lambda=0, 
                       symmetrize=TRUE, eps=1){
  
  if(method=='pseudo'){
    Lambda <- c(lambda)
    theta <- pseudo_mle(xi, L, Lambda)
  } else if(method=='mft'){
    theta <- mle_mft(xi=xi, L=L, dmax=dmax, eps=eps)
  } else stop('unknown method in mlestimate')

  m <- NCOL(xi)
  h <- matrix(0, nrow=m, ncol=L-1)
  for(i in seq_len(m)) h[i,] <- theta$h[[i]]
  J <- vector('list',m)
  for(i in seq_len(m)) J[[i]] <- vector('list',m)
  for(i in seq(1,m)) for(j in seq(i,m)){
    x <- matrix(theta$J[[i]][[j]], nrow=L-1, ncol=L-1, byrow=TRUE)
    xt <- matrix(theta$J[[j]][[i]], nrow=L-1, ncol=L-1, byrow=TRUE)
    if(i<j & symmetrize){ 
      x <- (x + t(xt))/2
      xt <- t(x)
    }
    J[[i]][[j]] <- x
    J[[j]][[i]] <- xt
  }
  
  return(list(h=h, J=J, mle=theta$lkl))
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
