#' Maximum likelihood estimate
#' 
#' Perform inference of bias and interaction parameters for a single response group 
#' 
#' Given numeric data matrix, either pseudo-likelihood
#' of mean-field theory is used to find the maximum likelihood estimate
#' of bias \code{h} and interaction \code{J} parameters.
#' 
#' @param xi Data matrix. In contrast to \code{data} slot in \code{\link{bbl}}
#'        object, this matrix is always numeric with elements ranging from 
#'        zero to positive integral upper bound. \code{\link{train}} does
#'        the transformation from factors to this numeric matrix for
#'        \code{bbl} objects of \code{type = 'factors'}. 
#' @param method \code{c('pseudo','mf')} for pseudo-likelihood maximization or
#'        mean field inference.
#' @param numeric Numeric model in which there is a single bias parameter
#'        for each predictor and a single interaction parameter for each
#'        distinct predictor pair. Elements of \code{xi} multipled by 
#'        these parameters determine the probability. If \code{FALSE} (the
#'        \code{'factors'} model), the numbers of parameters are one less than
#'        the number of factor levels and their products for bias and interaction,
#'        respectively. 
#' @param lambda Vector of L2 regularization parameters for 
#'        \code{method = 'pseudo'}. Inference will be repeated for each values 
#'        of \code{lambda}.
#' @param symmetrize Enforce the symmetry of interaction parameters by
#'        taking mean values of the matrix and its trace:
#'        \eqn{J_{ij}^{(y)}(x_1,x_2)=J_{ji}^{(y)}(x_2,x_1)}.
#' @param eps Vector of regularization parameters for \code{mf}. Must be
#'        within the range of \eqn{epsilon \in [0,1]}. Inferce will be 
#'        repeated for each value of \code{eps}.
#' @param nprint Frequency of printing iteration progress under \code{'pseudo'}.
#' @param itmax Maximum number of iterations for \code{'pseudo'}.
#' @param tolerance Upper bound for fractional changes in pseduo-likelihood
#'        values before termiating iteration in \code{'pseudo'}.
#' @param verbose Verbosity level.
#' @param prior.count Use prior count for \code{method = 'mf'} to reduce
#'        numerical instability.
#' @param naive Naive Bayes inference. Equivalent to \code{method = 'mf'} together
#'        with \code{eps = 0}.
#' @return List of inferred parameters \code{h} and \code{J}. See 
#'        \code{link{bbl}} for parameter structures.
#' @examples
#' set.seed(535)
#' predictors <- list()
#' for(i in 1:5) predictors[[i]] <- c('a','c','g','t')
#' par <- randompar(predictors)
#' par
#' xi <- sample_xi(nsample=5000, predictors=predictors, h=par$h, J=par$J,
#'                 code_out=TRUE)
#' head(xi)
#' ps <- mlestimate(xi=xi, method='pseudo')
#' ps$h
#' ps$J[[1]]
#' mf <- mlestimate(xi=xi, method='mf', eps=0.9)
#' plot(x=unlist(par$h), y=unlist(ps$h), xlab='True', ylab='Inferred')
#' segments(x0=-2, x1=2, y0=-2, y1=2, lty=2)
#' points(x=unlist(par$J), y=unlist(ps$J), col='red')
#' points(x=unlist(par$h), y=unlist(mf$h), col='blue')
#' points(x=unlist(par$J), y=unlist(mf$J), col='green')
#' @export

mlestimate <- function(xi, method='pseudo', numeric=FALSE, lambda=0, 
                       symmetrize=TRUE, eps=1, nprint=100, itmax=10000,
                       tolerance=1e-5, verbose=1, prior.count=TRUE,
                       naive=FALSE){
  
  m <- NCOL(xi)
  L <- apply(xi, 2, max)
  if(numeric) Lp <- rep(1, NCOL(xi))
  else Lp <- L
  
  if(method=='pseudo'){
    Lambda <- c(lambda)
    Nprint <- c(nprint)
    Itmax <- c(itmax)
    Tol <- c(tolerance)
    Verbose <- c(verbose)
    Naive <- c(naive)
    if(!is.numeric(xi[1,1])
       | min(xi)<0) stop('Input data to mlestimate must be numeric and non-negative')
    xi <- as.matrix(xi)
    theta <- pseudo_mle(xi, numeric, Lambda, Nprint, Itmax, Tol, Naive, Verbose)
    
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
  else if(method=='mf'){
    theta <- meanfield(xi=xi, L=L, Lp=Lp, eps=eps, numeric=numeric,
                       prior.count=prior.count)
    return(list(h=theta$h, J=theta$J, lz=theta$lz))
  } else stop('unknown method in mlestimate')

}

# Mean field inference
meanfield <- function(xi, L, Lp, eps=1, numeric=FALSE, prior.count=TRUE){

  nsite <- ncol(xi)
  nsample <- nrow(xi)
  
  bad <- apply(xi, 2, function(x){min(x)==max(x)})  # non-varying column
  id <- which(!bad)      # index in original vector of those that are good
  nvar <- nsite - sum(bad)
  xir <- xi[,!bad, drop=FALSE]

  if(eps>0){
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
        m0 <- colSums(csi) / nsample
    
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
    Jij <- solve(a=cijb, b=diag(-1, nrow=sum(Lp[!bad])))
  } else{
    mi <- vector('list',nvar)
    for(i in seq_len(nvar)){
      tmp <- NULL
      for(l in seq(Lp[id[i]])){
        mx <- sum(xir[,i]==l)
        if(prior.count)
          mx <- (mx + 1/(Lp[id[i]]+1))/(nsample+1)
        else mx <- mx/nsample
        tmp <- c(tmp,mx)
      }
      mi[[i]] <- tmp
    }
  }

  h <- vector('list',nsite)
  for(i in seq_len(nsite)){
    if(bad[i]) h[[i]] <- 0
    else{
      mx <- mi[[which(id==i)]]
      if(numeric)
        h[[i]] <- invav(mx,L[i])
      else
        h[[i]] <- log(mx/(1-sum(mx)))
    }
  }

  if(eps>0){  
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
  }

  J <- vector('list',nsite)
  Lsi <- 0
  for(i in seq_len(nsite))
    J[[i]] <- vector('list',nsite)
  for(i in seq_len(nsite)){
    if(bad[i] | eps==0){
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
  
  lz <- Zeff(L=L, bad=bad, m=nsite, h=h, J=J, xi=xi, numeric=numeric, mf=TRUE, 
             naive=(eps==0))  # log partition function
  
  return(list(h=h, J=J, lz=lz))
}

# Inverse function of mean single-site average mi

invav <- function(mx, L){
  
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
  
  zm <- NULL
  for(m in mx){
    if(m==0) hi <- -Inf
    else if(m==L) hi <- Inf
    else hi <- stats::uniroot(f, interval=c(-10,10), m=m, L=L)$root
    zm <- c(zm, hi)
  }
  return(zm)
}

# Pseudo partition function using data xi
Zeff <- function(L, bad, m, h, J, xi, numeric, mf=mf, naive){
  
  nsample <- NROW(xi)
  if(numeric) Lp <- rep(1, m)
  else Lp <- L
  if(mf){
    fi <- vector('list',m)
    f0 <- rep(1, m)
    for(i in seq_len(m)){
      if(bad[i]) next()
      if(L[i]<1){
        fi[[i]] <- 0
        f0[i] <- 1
      } else{
        fi[[i]] <- rep(1/(L[i]+1), L[i])
        for(l in seq_len(L[i]))
          fi[[i]][l] <- fi[[i]][l] + sum(xi[,i]==l)
        fi[[i]] <- fi[[i]]/(nsample+1)
        f0[i] <- 1 - sum(fi[[i]])
      }
    }
    lz <- - sum(log(f0))
    if(naive) return(lz)
    for(i in seq(1,m-1)){
      if(bad[i]) next()
      for(j in seq(i+1,m)){
        if(bad[j]) next()
        for(l0 in seq_len(Lp[i])) for(l1 in seq_len(Lp[j])){
          if(NROW(J[[i]][[j]])<l0 |
             NCOL(J[[i]][[j]])<l1) next()
          dz <- J[[i]][[j]][l0,l1]*fi[[i]][l0]*fi[[j]][l1]
          if(numeric) dz <- dz * l0*l1;
          lz <- lz - dz
        }
      }
    }
    return(lz)
  }
  
  lz <- 0
  for(k in seq_len(nsample)) for(i in seq_len(m)){
    if(bad[i]) next()
    z <- 1
    for(a in seq(1, L[i])){
      if(numeric) e <- h[[i]][1]*a
      else if(length(h[[i]]) < a) next()
      else e <- h[[i]][a]
      if(naive) next()
      for(j in seq(1,m)){
        if(j==i) next()
        xj <- xi[k,j]
        if(xj==0) next()
        if(numeric)
          e <- e + J[[i]][[j]]*a*xj/2
        else if(NROW(J[[i]][[j]])<a | NCOL(J[[i]][[j]])<xj) next()
        else e <- e + J[[i]][[j]][a,xj]/2
      }
      z <- z + exp(e)
    }
    lz <- lz + log(z)
  }
  lz <- lz/nsample
  
  return(lz)
}
