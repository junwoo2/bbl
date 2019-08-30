#' Train Boltzmann Bayes model
#' 
#' Perform inference on data provided and store parameters
#' 
#' Use \code{data} stored in \code{bbl} object as the training data, 
#' performs either \code{pseudo} or \code{mf} inference for each response
#' group, and stores the model parameters in slots \code{h} and \code{J}. 
#' This function is a driver of \code{\link{mlestimate}} for \code{bbl} objects.  
#' 
#' @param object Object of class \code{bbl} containing training data.
#' @param method \code{c('pseudo','mf')} for pseudo-likelihood maximization or
#'        mean field.
#' @param naive Naive Bayes. Equivalent to using \code{method = 'mf'} and
#'        \code{eps = 0} in \code{\link{mlestimate}}.
#' @param verbose Verbosity level. Relays it to \code{\link{mlestimate}} 
#'        with one level lower.
#' @param fixL Use predictor level sizes \code{L} from 
#'        \code{object@predictors}. If \code{FALSE}, \code{L} is inferred 
#'        from \code{object@data}. Intended for cases where this function
#'        is called from \code{\link{crossval}}.
#' @param ... Other parameters for \code{\link{mlestimate}}.
#' @return \code{object} with slots \code{h}, \code{J}, and \code{lz} filled.
#' @examples
#' titanic <- freq2raw(as.data.frame(Titanic), Freq='Freq')
#' model <- bbl(data=titanic, y='Survived')
#' model <- train(model, method='pseudo')
#' model@h
#' model@J[[1]]
#' @export
train <- function(object, method='pseudo', naive=FALSE, verbose=1, 
                  fixL=FALSE, ...){

  predictors <- object@predictors
  groups <- object@groups
  Ly <- length(groups)
  if(Ly < 2) stop('No. of groups < 2')
  
  if(naive)
    method <- 'mf'

  data <- object@data
  xi <- data[,colnames(data)!=object@y]
  y <- data[,colnames(data)==object@y]
  
  nrl <- 0
  object@h <- object@J <- vector('list',Ly)
  lz <- c()
  for(iy in seq_len(Ly)){
    id <- which(y==groups[iy])
    if(length(id)==0) stop(paste0('No instance of "',
                                  groups[iy],'" in training data'))
    xid <- xi[id,]
    if(verbose>0) cat('  Inference for class "',object@y,'" = ',
                      groups[iy],':\n',sep='')
    xidi <- matrix(0, nrow=NROW(xid), ncol=NCOL(xid))
    if(fixL) L <- c(0, NCOL(xid))
    else L <- NULL
    for(i in seq_len(NCOL(xidi))){
      xidi[,i] <- match(xid[,i],predictors[[i]])-1   
                                      # xidi = 0, ..., L-1
      if(fixL) L[i] <- length(predictors[[i]])
    }
    mle <- mlestimate(xi=xidi, method=method, L=L, naive=naive, 
                      verbose=verbose-1, ...)
    if(verbose>0 & method=='pseudo') 
      cat('  Maximum pseudo-likelihood = ',mle$mle,'\n\n',sep='')
    
    m <- length(object@predictors)
    
    for(i in seq_len(m)){
      names(mle$h[[i]]) <- 
        object@predictors[[i]][-1][seq_along(mle$h[[i]])]
      if(naive) next()
      for(j in seq_len(m)){
        if(NROW(mle$J[[i]][[j]])>0)
          rownames(mle$J[[i]][[j]]) <- 
            object@predictors[[i]][-1][seq_len(NROW(mle$J[[i]][[j]]))]
        if(NCOL(mle$J[[i]][[j]])>0)
          colnames(mle$J[[i]][[j]]) <- 
            object@predictors[[j]][-1][seq_len(NCOL(mle$J[[i]][[j]]))]
      }
    }
    object@h[[iy]] <- mle$h
    object@J[[iy]] <- mle$J
    
    lz <- c(lz,mle$lz)
  }
  
  object@lz <- lz
  object@data <- data
  return(object)
}
