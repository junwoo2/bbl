#' Main S3 function for bbl inference
#' @export
bbl <- function(formula, data, weights=NULL, verbose=1, ...){

  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  term <- terms(formula, data=data)
  xlevels <- .getXlevels(term, m=data)
  m <- length(xlevels)
  idy <- attributes(term)$response
  resp <- all.vars(cl)[idy]
  y <- data[,resp]
  x <- data[,names(xlevels)]

  label <- attr(term,'term.labels')
  ilabel <- label[vapply(label,FUN=function(x){grepl(':',x)}, logical(1))]
  ijlabel <- strsplit(ilabel,split=':')
  qJ <- matrix(FALSE, nrow=m, ncol=m)
  rownames(qJ) <- colnames(qJ) <- names(xlevels)
  for(k in seq_along(ijlabel))
    qJ[ijlabel[[k]][1],ijlabel[[k]][2]] <- TRUE
  qJ <- qJ | t(qJ)    # TRUE for all interacting pairs of predictors
  
  b <- bbl.fit(x=x, y=y, qJ=qJ, weights=weights, xlevels=xlevels,
               verbose=verbose-1, ...)
  bb <- list()
  class(bb) <- 'bbl'
  bb$coefficients <- list(h=b$h, J=b$J)
  bb$xlevels <- xlevels
  bb$terms <- term
  bb$groups <- levels(factor(y))
  bb$groupname <- resp
  bb$qJ <- qJ
  bb$model <- data[,c(resp, names(xlevels))]
  bb$lkh <- b$lkh
  bb$lz <- b$lz
  bb$weights <- weights
  bb$call <- cl
  
  return(bb)
}


#' @export
bbl.fit <- function(x, y, qJ, weights=NULL, xlevels=NULL, verbose=1, ...){
  
  le <- levels(factor(y))
  Ly <- length(le)
  h <- J <- list()
  m <- NCOL(x)
  if(is.null(xlevels)){
    xlevels <- list()
    for(i in seq_len(m))
      xlevels[[colnames(x)[m]]] <- levels(factor(x[,i]))
  }
  lkh <- 0
  lz <- rep(0, Ly)
  for(iy in seq_len(Ly)){
    da <- x[y==le[iy],]
    ny <- sum(y==le[iy])
    if(!is.null(weights))
      we <- weights[y==le[iy]]
    else we <- rep(1, sum(y==le[iy]))
    if(verbose > 1) cat(' Inference for y = ',le[iy],':\n')
    xda <- data.matrix(da) - 1
    b <- mlestimate(xi=xda, weights=we, qJ=qJ, verbose=verbose-1, ...)
    for(i in seq_len(m)){
      ni <- xlevels[[i]][-1][seq_along(b$h[[i]])]
      names(b$h[[i]]) <- ni
      for(j in seq_len(m)){
        if(NROW(b$J[[i]][[j]])>0)
          rownames(b$J[[i]][[j]]) <- ni[seq_len(NROW(b$J[[i]][[j]]))]
        if(NCOL(b$J[[i]][[j]])>0)
          colnames(b$J[[i]][[j]]) <- 
            xlevels[[j]][-1][seq_len(NCOL(b$J[[i]][[j]]))]
      }
    }
    h[[iy]] <- b$h
    J[[iy]] <- b$J
    lkh <- lkh + b$lkh*ny
    lz[iy] <- b$lz
  }
  names(lz) <- le
  
  return(list(h=h, J=J, lkh=lkh, lz=lz))
}