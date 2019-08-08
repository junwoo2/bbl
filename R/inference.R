#' Train cRfm model 
#' @param data Data frame of training data. The column named \code{y} is
#'             used as group variables
#' @export
train <- function(object, method='pseudo', lambda=0, data=NULL, 
                  verbose=1){
  
  predictors <- object@predictors
  groups <- object@groups
  L <- length(predictors)
  Ly <- length(groups)
  nsite <- object@nsite
  
  if(is.null(data))
    data <- object@data
  else if(!'y' %in% colnames(data)) stop("Data must contain column 'y'")
  
  xi <- data[,colnames(data)!='y']
  y <- data$y
  
  if(nsite!=NCOL(xi)) stop('No. of columns in data does not match nsite')
  
  object@h <- object@J <- vector('list',Ly)
  
  for(iy in seq_len(Ly)){
    id <- which(y==groups[iy])
    xid <- xi[id,]
    if(verbose>0) cat('Inference for class y = ',groups[iy],':\n',sep='')
    xidi <- matrix(0, nrow=NROW(xid), ncol=NCOL(xid))
    for(i in seq_len(NCOL(xidi)))
      xidi[,i] <- match(xid[,i],predictors) - 1   # xidi = 0, ..., L-1
    mle <- mlestimate(xi=xidi, L=L, lambda=lambda, method=method, verbose=verbose)
    if(verbose>0) cat('Maximum pseudo-likelihood = ',mle$mle,'\n\n',sep='')
    object@h[[iy]] <- mle$h
    object@J[[iy]] <- mle$J
  }
  
  object@data <- data
  
  return(object)
  
} 