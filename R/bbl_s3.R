#' Main S3 function for bbl inference
#' @export
bbl <- function(formula, data, freq=NULL, xlevels=NULL,
                verbose=1, method='pseudo', 
                novarOk=FALSE, testNull=TRUE, ...){

  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  if(!is.null(freq)){
    if(length(freq)!=NROW(data))
      stop('Length of freq does not match data')
    zero <- freq==0
    data <- data[!zero,]
    freq <- freq[!zero]   # remove rows with zero freq
  }
  
  term <- terms(formula, data=data)
  idy <- attributes(term)$response
  vars <- as.character(attributes(term)$variables)
  resp <- vars[[idy+1]]
  vars <- vars[!(vars %in% c('list',resp))]
  # xlevels <- .getXlevels(term, m=data)
  if(is.null(xlevels))
    xlevels <- getxlevels(vars, data=data)
  
  formula <- formula(term)
  
  if(!novarOk){
    for(i in seq_along(xlevels)){
      if(length(xlevels[[i]])==1)
        stop(paste0('Predictor ',names(xlevels)[i],' has one level'))
    }
  }
  m <- length(xlevels)
  idy <- attributes(term)$response
  resp <- as.character(attributes(term)$variables[[idy+1]])
  y <- data[,resp]
  x <- data[,names(xlevels)]
  
  Ly <- length(levels(factor(y)))
  if(Ly==1) warning('Only one response group in data')
  
  label <- attr(term,'term.labels')
  ilabel <- label[vapply(label,FUN=function(x){grepl(':',x)}, logical(1))]
  ilabel <- gsub('`','',ilabel)
  ijlabel <- strsplit(ilabel,split=':')
  qJ <- matrix(FALSE, nrow=m, ncol=m)
  rownames(qJ) <- colnames(qJ) <- names(xlevels)
  for(k in seq_along(ijlabel))
    qJ[ijlabel[[k]][1],ijlabel[[k]][2]] <- TRUE
  qJ <- qJ | t(qJ)    # TRUE for all interacting pairs of predictors
  naive <- sum(qJ)==0
  
  b <- bbl.fit(x=x, y=y, qJ=qJ, freq=freq, xlevels=xlevels,
               verbose=verbose-1, method=method, ...) # alternative
  if(testNull)
    b0 <- bbl.fit(x=x, y=rep('pooled',length(y)), qJ=qJ, freq=freq, 
                xlevels=xlevels, verbose=verbose-1, method=method, ...) # null
  else b0 <- NULL
  bb <- list()
  class(bb) <- 'bbl'
  bb$coefficients <- list(h=b$h, J=b$J, h0=b0$h[[1]], J0=b0$J[[1]])
  bb$xlevels <- xlevels
  bb$terms <- term
  bb$groups <- levels(factor(y))
  bb$groupname <- resp
  bb$qJ <- qJ
  bb$model <- data[,c(resp, names(xlevels))]
  bb$lkh <- b$lkh
  bb$lz <- b$lz
  bb$freq <- freq
  bb$call <- cl
  if(naive) bb$method <- 'mf'
  else bb$method <- method
  df <- 0
  for(i in seq_len(m)){
    ni <- length(xlevels[[i]])-1
    df <- df + ni
    if(i==m) next()
    for(j in seq(i+1,m))
      if(qJ[i,j]) df <- df + ni*(length(xlevels[[j]])-1)
  }
  bb$df <- df
  
  return(bb)
}

#' @export
bbl.fit <- function(x, y, qJ, freq=NULL, xlevels=NULL, verbose=1, 
                    method='pseudo', ...){
  
  le <- levels(factor(y))
  Ly <- length(le)
  h <- J <- list()
  m <- NCOL(x)
  if(is.null(xlevels)){
    xlevels <- list()
    for(i in seq_len(m)){
      v <- x[,i]
      if(is.factor(v) | is.character(v))
        xlevels[[colnames(x)[m]]] <- levels(factor(v))
      else{
        v <- unique(v)
        xlevels[[colnames(x)[m]]] <- v[order(v)]
      }
    }
  }
  numericx <- unlist(lapply(X=xlevels, FUN=is.numeric))
      # true for numeric predictors
  for(i in seq_along(numericx))
    if(numericx[i]) if(min(x[,i])!=0) 
      stop(paste0("Numeric '",names(xlevels)[i],"' minimum is nonzero"))
                       
  lkh <- 0
  lz <- rep(0, Ly)
  for(iy in seq_len(Ly)){
    ny <- sum(y==le[iy])
    if(ny==0) 
      stop(paste0('No instance of "',groups[iy],'" in training data'))
    da <- x[y==le[iy],]
    if(!is.null(freq))
      frq <- freq[y==le[iy]]
    else frq <- rep(1, sum(y==le[iy]))
    if(verbose > 1) cat(' Inference for y = ',le[iy],':\n',sep='')
#   xda <- data.matrix(da) - 1
    xda <- matrix(0, nrow=NROW(da), ncol=NCOL(da))
    L <- NULL
    for(i in seq_len(NCOL(da))){
      xda[,i] <- match(da[,i],xlevels[[i]])-1
      L[i] <- length(xlevels[[i]])
    }
    b <- mlestimate(xi=xda, freq=frq, qJ=qJ, numericx, verbose=verbose-1, 
                    method=method, L=L, ...)
    names(b$h) <- names(b$J) <- names(xlevels)
    for(i in seq_len(m)){
      if(numericx[i]) ni <- 'Numeric'
      else ni <- xlevels[[i]][-1][seq_along(b$h[[i]])]
      names(b$h[[i]]) <- ni
      names(b$J[[i]]) <- names(xlevels)
      for(j in seq_len(m)){
        if(NROW(b$J[[i]][[j]])>0){
          if(numericx[i])
            rownames(b$J[[i]][[j]]) <- 'Numeric'
          else
            rownames(b$J[[i]][[j]]) <- ni[seq_len(NROW(b$J[[i]][[j]]))]
        }
        if(NCOL(b$J[[i]][[j]])>0){
          if(numericx[j])
            colnames(b$J[[i]][[j]]) <- 'Numeric'
          else
            colnames(b$J[[i]][[j]]) <- 
              xlevels[[j]][-1][seq_len(NCOL(b$J[[i]][[j]]))]
        }
      }
    }
    h[[iy]] <- b$h
    J[[iy]] <- b$J
    lkh <- lkh + b$lkh*ny
    lz[iy] <- b$lz
  }
  names(lz) <- names(h) <- names(J) <- le

  return(list(h=h, J=J, lkh=lkh, lz=lz))
}