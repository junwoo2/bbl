#' Train cRfm model 
#' @param lambda Penalizer
#' @param object \code{bbm} object for model
#' @param data Data frame of training data. The column named \code{y} is
#'             used as group variables
#' @param naive Naive Bayes (no interaction; \code{J}=0)
#' @param verbose Output verbosity
#' @export
train <- function(object, method='pseudo', L=NULL, lambda=0, eps=1, data=NULL, 
                  verbose=1){
  
  predictors <- object@predictors
  groups <- object@groups
  Ly <- length(groups)
  if(Ly < 2) stop('No. of groups < 2')
  
  if(is.null(data))
    data <- object@data
  else if(!object@y %in% colnames(data)) 
    stop(paste0('Data must contain column ',object@y,'\n'))
  
  xi <- data[,colnames(data)!=object@y]
  y <- data[,colnames(data)==object@y]
  
  nrl <- 0
  object@h <- object@J <- vector('list',Ly)
  lz <- c()
  for(iy in seq_len(Ly)){
    id <- which(y==groups[iy])
    xid <- xi[id,]
    if(verbose>0) cat('  Inference for class "',object@y,'" = ',
                      groups[iy],':\n',sep='')
    if(object@type=='numeric'){ 
      if(!is.numeric(xid[1,1])) stop('Numeric model requires numeric data')
      xidi <- as.matrix(xid)
    } else{
      xidi <- matrix(0, nrow=NROW(xid), ncol=NCOL(xid))
      for(i in seq_len(NCOL(xidi)))
        xidi[,i] <- match(xid[,i],predictors[[i]])-1   # xidi = 0, ..., L-1
    }
    mle <- mlestimate(xi=xidi, L=L, numeric=object@type=='numeric', 
                      lambda=lambda, method=method, eps=eps, verbose=verbose-1)
    if(verbose>0 & method=='pseudo') 
      cat('  Maximum pseudo-likelihood = ',mle$mle,'\n\n',sep='')
    
    m <- length(object@predictors)
    
    for(i in seq(m)){
      names(mle$h[[i]]) <- object@predictors[[i]][-1][seq_along(mle$h[[i]])]
      for(j in seq(m)){
        if(length(mle$J[[i]][[j]]==0)) next()
        rownames(mle$J[[i]][[j]]) <- object@predictors[[i]][-1][seq_along(NROW(mle$J[[i]][[j]]))]
        colnames(mle$J[[i]][[j]]) <- object@predictors[[j]][-1][seq_along(NCOL(mle$J[[i]][[j]]))]
      }
    }
    object@h[[iy]] <- mle$h
    object@J[[iy]] <- mle$J
    
    if(method=='pseudo') lz <- c(lz,mle$lz)
  }
  
  if(method=='pseudo') object@lz <- lz
  object@data <- data
  return(object)
} 
