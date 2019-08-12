#' Train cRfm model 
#' @param lambda Penalizer
#' @param object \code{bbm} object for model
#' @param data Data frame of training data. The column named \code{y} is
#'             used as group variables
#' @param naive Naive Bayes (no interaction; \code{J}=0)
#' @param verbose Output verbosity
#' @export
train <- function(object, method='pseudo', lambda=0, L=NULL, eps=1,
                  data=NULL, naive=FALSE, verbose=1, nullcount=0){
  
  predictors <- object@predictors
  groups <- object@groups
  if(identical(predictors, 'numeric')){
    numeric <- TRUE
    if(is.null(L)) stop('Numeric model requires L input')
  }else{
    numeric <- FALSE
    L <- length(predictors)
  }
  Ly <- length(groups)
  nsite <- object@nsite
  
  if(is.null(data))
    data <- object@data
  else if(!'y' %in% colnames(data)) stop("Data must contain column 'y'")
  
  xi <- data[,colnames(data)!='y']
  if(numeric) if(max(xi)!=L-1) 
    stop('Numeric model data have maximum must have range of [0,L-1]')
  y <- data$y
  
  if(nsite!=NCOL(xi)) stop('No. of columns in data does not match nsite')
  
  object@h <- object@J <- vector('list',Ly)
  lz <- c()
  for(iy in seq_len(Ly)){
    id <- which(y==groups[iy])
    xid <- xi[id,]
    if(verbose>0) cat('  Inference for class y = ',groups[iy],':\n',sep='')
    if(numeric){ 
      if(!is.numeric(xid[1,1])) stop('Numeric model requires numeric data')
      xidi <- as.matrix(xid)
    } else{
      xidi <- matrix(0, nrow=NROW(xid), ncol=NCOL(xid))
      for(i in seq_len(NCOL(xidi)))
        xidi[,i] <- match(xid[,i],predictors) - 1   # xidi = 0, ..., L-1
    }
    mle <- mlestimate(xi=xidi, L=L, numeric=numeric, lambda=lambda, 
                      method=method, eps=eps, verbose=verbose-1, naive=naive,
                      nullcount=nullcount)
    if(verbose>0 & method=='pseudo') 
      cat('  Maximum pseudo-likelihood = ',mle$mle,'\n\n',sep='')
    object@h[[iy]] <- mle$h
    object@J[[iy]] <- mle$J
    if(method=='pseudo') lz <- c(lz,mle$lz)
  }
  
  if(method=='pseudo') object@lz <- lz
  object@data <- data
  return(object)
} 
