#' Train cRfm model 
#' @export
train <- function(object, method='pseudo', lambda=0, xi=NULL, y=NULL, 
                  verbose=1){
  
  L <- object@L
  Ly <- object@Ly
  nsite <- object@nsite
  
  if(is.null(xi)) xi <- object@xi
  if(is.null(y)) y <- object@y
  
  if(nsite!=NCOL(xi)) stop('No. of columns in data does not match nsite')
  
  object@h <- object@J <- vector('list',Ly)
  
  for(iy in seq_len(Ly)){
    id <- which(y==iy-1)
    xid <- xi[id,]
    if(verbose>0) cat('Inference for class y = ',iy-1,':\n',sep='')
    mle <- mlestimate(xi=xid, L=L, lambda=lambda, method=method)
    if(verbose>0) cat('Maximum pseudo-likelihood = ',mle$mle,'\n\n',sep='')
    object@h[[iy]] <- mle$h
    object@J[[iy]] <- mle$J
  }
  
  object@xi <- xi
  object@y <- y
  
  return(object)
  
} 