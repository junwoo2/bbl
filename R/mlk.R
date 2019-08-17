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

mlestimate <- function(xi, method='pseudo', L=NULL, numeric=FALSE, lambda=0, 
                       symmetrize=TRUE, eps=1, nprint=100, itmax=10000,
                       tolerance=1e-5, verbose=1, prior.count=TRUE){
  
  m <- NCOL(xi)
  if(is.null(L))
    L <- apply(xi, 2, max)
  else L <- rep(L-1, m)
  if(numeric) Lp <- rep(1, NCOL(xi))
  else Lp <- L
  
  if(method=='pseudo'){
    Lambda <- c(lambda)
    Nprint <- c(nprint)
    Itmax <- c(itmax)
    Tol <- c(tolerance)
    Verbose <- c(verbose)
    if(!is.numeric(xi[1,1])) stop('Input data to mlestimate must be numeric')
    xi <- as.matrix(xi)
    theta <- pseudo_mle(xi, numeric, Lambda, Nprint, Itmax, Tol, Verbose)
    
    h <- theta$h
    J <- vector('list',m)
    for(i in seq_len(m)) J[[i]] <- vector('list',m)
    for(i in seq(1,m)){ 
      for(j in seq(i,m)){
        x <- matrix(theta$J[[i]][[j]], nrow=Lp[i], ncol=Lp[j], byrow=TRUE)
        xt <- matrix(theta$J[[j]][[i]], nrow=Lp[j], ncol=Lp[i], byrow=TRUE)
        if(i<j & symmetrize){ 
          x <- (x + t(xt))/2
          xt <- t(x)
        }
        J[[i]][[j]] <- x
        J[[j]][[i]] <- xt
      }
    }
    return(list(h=h, J=J, mle=theta$lkl, lz=theta$lz))
  }
  else if(method %in% c('mf','nb')){
    if(method=='nb') eps <- 0  # no interaction
    theta <- meanfield(xi=xi, L=L, Lp=Lp, eps=eps, numeric=numeric,
                       prior.count=prior.count)
    return(list(h=theta$h, J=theta$J))
  } else stop('unknown method in mlestimate')

}

meanfield <- function(xi, L, Lp, eps=1, numeric=FALSE, prior.count=TRUE){

  nsite <- ncol(xi)
  nsample <- nrow(xi)
  
  bad <- colSums(xi)==0  # predictors without variations
  id <- which(!bad)      # index in original vector of those that are good
  nvar <- nsite - sum(bad)
  xir <- xi[,!bad]
  
  csi <- matrix(0, nrow=nsample, ncol=sum(Lp[!bad]))
  for(k in seq_len(nsample)){ 
    Ls <- 0
    for(i in seq_len(nvar)){
      if(xir[k,i]>0){
        if(numeric) csi[k, Ls+1] <- xir[k,i]
        else csi[k, Ls+xir[k,i]] <- 1
      }
      Ls <- Ls + Lp[id[i]]
    }
  }

  if(prior.count){
    count0 <- NULL
    for(i in seq_len(nvar)){
      w <- Lp[id[i]]
      count0 <- c(count0, rep(1/(w+1),w))
    }
    m0 <- (count0 + colSums(csi))/(nsample+1)
  } else
    m0 <- colSums(csi) /nsample
  mi <- vector('list',nvar)
  Ls <- 0
  for(i in seq_len(nvar)){
    mi[[i]] <- m0[seq(Ls+1,Ls+Lp[id[i]])]
    Ls <- Ls + Lp[id[i]]
  }
  csi <- csi - matrix(m0, nrow=nsample, ncol=Ls, byrow=TRUE)
  if(prior.count){
    csi <- csi + matrix(count0/nsample, nrow=nsample, ncol=NCOL(csi), byrow=TRUE)
    cij <- t(csi) %*% csi / (nsample + 1)
  }
  else
    cis <- t(csi) %*% csi / nsample
  
  cijb <- eps*cij + (1-eps)*mean(diag(cij))*diag(1,nrow=NROW(cij))

  if(!numeric){
    h <- vector('list',nsite)
    for(i in seq_len(nsite)){
      if(bad[i]) h[[i]] <- 0
      else{
        mx <- mi[[which(id==i)]]
        h[[i]] <- log(mx/(1-sum(mx)))
      }
    }
  } else
    h <- invav(mi,L, bad)

  Jij <- solve(a=cijb, b=diag(-1, nrow=sum(Lp[!bad])))
  diag(Jij) <- 0
  Lsi <- 0
  for(i in seq_len(nsite)){
    if(bad[i]) next()
    for(l in seq_len(Lp[i])){
      Lsj <- 0
      for(j in seq_len(nsite)){
        if(i!=j){
          if(bad[j]) next()
          mx <- mi[[which(id==j)]]
          for(q in seq_len(Lp[j]))
            h[[i]][l] <- h[[i]][l] - Jij[Lsi+l, Lsj+q]*mx[q]
        }
        Lsj <- Lsj + Lp[j]
      }
    }
    Lsi <- Lsi + Lp[i]
  }

  J <- vector('list',nsite)
  Lsi <- 0
  for(i in seq_len(nsite)){
    J[[i]] <- vector('list',nsite)
    if(bad[i]){
      for(j in seq(nsite)) J[[i]][[j]] <- matrix(0)
      next()
    }
    Lsj <- 0
    for(j in seq_len(nsite)){
      if(bad[j]){ 
        J[[i]][[j]] <- matrix(0)
        next()
      }
      w <- matrix(0, nrow=Lp[i], ncol=Lp[j])
      if(i==j) J[[i]][[i]] <- w
      else{ 
        for(l in seq_len(Lp[i])) for(q in seq_len(Lp[j]))
          w[l,q] <- Jij[Lsi+l,Lsj+q]
        J[[i]][[j]] <- w
        J[[j]][[i]] <- t(w)
      }
      Lsj <- Lsj + Lp[j]
    }
    Lsi <- Lsi + Lp[i]
  }
  
  return(list(h=h, J=J))
}

# Inverse function of mean single-site average mi

invav <- function(mi, L, bad){
  
  f <- function(h, m, L){
    up <- down <- 0
    eh <- exp(h)
    for(x in seq(0,L)){
      ehx <- eh^x
      up <- up + x*ehx
      down <- down + ehx
    }
    m-up/down
  }
  
  z <- list()
  id <- which(!bad)
  for(i in seq_along(mi)){ 
    if(bad[i]) zm <- rep(0, length(mi[[i]]))
    else{
      zm <- NULL
      for(m in mi[[id==i]]){
        if(m==0) hi <- -Inf
        else if(m==L[i]) hi <- Inf
        else hi <- uniroot(f, interval=c(-10,10), m=m, L=L[i])$root
        zm <- c(zm, hi)
      }
    }
    z[[i]] <- zm
  }
  return(z)
}