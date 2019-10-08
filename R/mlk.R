#' Maximum likelihood estimate
#' 
#' Perform inference of bias and interaction parameters for a single response group 
#' 
#' Given numeric data matrix, either pseudo-likelihood
#' of mean-field theory is used to find the maximum likelihood estimate
#' of bias \code{h} and interaction \code{J} parameters. Normally
#' called by \code{\link{train}} rather than directly.
#' 
#' @param xi Data matrix. In contrast to \code{data} slot in \code{\link{bbl}}
#'        object, this matrix is always numeric with elements ranging from 
#'        zero to positive integral upper bound \code{L-1}. 
#'        \code{\link{train}} does
#'        the transformation from factors to this numeric matrix for
#'        \code{bbl} objects. 
#' @param method \code{c('pseudo','mf')} for pseudo-likelihood maximization or
#'        mean field inference.
#' @param L Vector of number of factor levels in each predictor. If
#'        \code{NULL}, will be inferred from \code{xi}.
#' @param lambda Vector of L2 regularization parameters for 
#'        \code{method = 'pseudo'}. Applies to interaction parameters \code{J}.
#' @param lambdah L2 parameters for \code{h} in \code{'pseudo'}.
#'         If \code{NULL}, it is set equal to \code{lambda}.
#'         \code{lambdah = 0} will free \code{h} from penalization.
#' @param symmetrize Enforce the symmetry of interaction parameters by
#'        taking mean values of the matrix and its trace:
#'        \eqn{J_{ij}^{(y)}(x_1,x_2)=J_{ji}^{(y)}(x_2,x_1)}.
#' @param eps Vector of regularization parameters for \code{mf}. Must be
#'        within the range of \eqn{\epsilon \in [0,1]}.
#' @param nprint Frequency of printing iteration progress under \code{'pseudo'}.
#' @param itmax Maximum number of iterations for \code{'pseudo'}.
#' @param tolerance Upper bound for fractional changes in pseduo-likelihood
#'        values before termiating iteration in \code{'pseudo'}.
#' @param verbose Verbosity level.
#' @param prior.count Use prior count for \code{method = 'mf'} to reduce
#'        numerical instability.
#' @param naive Naive Bayes inference. Equivalent to \code{method = 'mf'} together
#'        with \code{eps = 0}.
#' @param lz.half Divide interaction term in approximation to \eqn{\ln Z_{iy}}
#'        in \code{'pseudo'}.
#' @return List of inferred parameters \code{h} and \code{J}. See 
#'        \code{\link{bbl}} for parameter structures.
#' @examples
#' set.seed(535)
#' predictors <- list()
#' for(i in 1:5) predictors[[i]] <- c('a','c','g','t')
#' par <- randompar(predictors)
#' par
#' xi <- sample_xi(nsample=5000, predictors=predictors, h=par$h, J=par$J,
#'                 code_out=TRUE)
#' head(xi)
#' ps <- mlestimate(xi=xi, method='pseudo', lambda=0)
#' ps$h
#' ps$J[[1]]
#' mf <- mlestimate(xi=xi, method='mf', eps=0.9)
#' plot(x=unlist(par$h), y=unlist(ps$h), xlab='True', ylab='Inferred')
#' segments(x0=-2, x1=2, y0=-2, y1=2, lty=2)
#' points(x=unlist(par$J), y=unlist(ps$J), col='red')
#' points(x=unlist(par$h), y=unlist(mf$h), col='blue')
#' points(x=unlist(par$J), y=unlist(mf$J), col='green')
#' @export

mlestimate <- function(xi, method='pseudo', L=NULL, lambda=1e-5,
                       lambdah=0, symmetrize=TRUE, eps=0.9, 
                       nprint=100, itmax=10000,
                       tolerance=1e-5, verbose=1, prior.count=TRUE,
                       naive=FALSE, lz.half=FALSE){
  
  if(is.null(lambdah))
    lambdah <- lambda
  
  m <- NCOL(xi)
  La <- apply(xi, 2, max)
  if(is.null(L))
    L <- La
  else{ 
    if(!all(L >= La)) 
      stop('Data provided have predictor levels exceeding L')
    L <- L-1
  }
  
  if(naive){
    method <- 'mf'
    eps <- 0
  }
  
  xi <- as.matrix(xi)
  if(!is.numeric(xi[1,1]) | min(xi)<0) 
    stop('Input data to mlestimate must be numeric and non-negative')
  
  if(method=='pseudo'){
    Lambda <- c(lambda, lambdah)
    Nprint <- c(nprint)
    Itmax <- c(itmax)
    Tol <- c(tolerance)
    Verbose <- c(verbose)
    Lzhalf <- c(lz.half)
    Naive <- c(naive)
    theta <- pseudo_mle(xi, L, Lambda, Nprint, Itmax, Tol, Naive, 
                        Verbose, Lzhalf)
    L <- theta$L
  }
  else if(method=='mf'){
    Eps <- c(eps)
    theta <- mfwrapper(xi, L, Eps)
  }
  else stop('unknown method in mlestimate')

  h <- theta$h
  J <- vector('list',m)
  for(i in seq_len(m)) J[[i]] <- vector('list',m)
  for(i in seq(1,m)){ 
    for(j in seq(i,m)){
      x <- matrix(theta$J[[i]][[j]], nrow=L[i], ncol=L[j], byrow=TRUE)
      xt <- matrix(theta$J[[j]][[i]], nrow=L[j], ncol=L[i], byrow=TRUE)
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