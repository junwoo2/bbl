#' Predict class from new data using bbm model
#' 
#' @param newdata List of names \code{xi} and \code{y}; new data for 
#'                which prediction is made
#' @param logit Return predictors whose logistic function gives probability;
#'              otherwise return probability itself
#' @return Matrix of predictors/posterior proabilities
#' @export
setMethod('predict', 'bbm', function(object, newdata=NULL, logit=TRUE){
  
#  browser()
  if(is.null(newdata)){
    xi <- object@xi
    y <- object@y
  } else{
    xi <- newdata[[1]]
    y <- newdata[[2]]
  }
  
  Ly <- object@Ly
  h <- object@h
  J <- object@J
  nsample <- NROW(xi)
  lz <- py <- rep(0, Ly)
  
  for(iy in seq_len(Ly)){
    lz[iy] <- Zeff(L, m, h[[iy]], J[[iy]], xi)  # log partition function
    py[iy] <- sum(y==iy-1)            # marginal distribution P(y)
  }
  py <- py/nsample
  
  ay <- matrix(0, nrow=nsample, ncol=Ly)
  for(k in seq_len(nsample)){
    x <- xi[k,]
    for(iy in seq_len(Ly-1)){
      w0 <- ham(x, h[[iy]], J[[iy]]) - lz[iy] + log(py[iy])
      e <- 0
      for(iyp in seq_len(Ly)){
        if(iyp==iy) next()
        w <- ham(x, h[[iyp]], J[[iyp]]) - lz[iyp] + log(py[iyp])
        e <- e + exp(w - w0)
      }
      ay[k, iy] <- -log(e)
    }
    pry <- 1-sum(1/(1+exp(-ay[k,-Ly])))  # sum rule for sum_y {Pr(y|x)}=1
    ay[k, Ly] <- log(pry/(1-pry))
  }
  
  if(!logit) ay <- 1/(1+exp(-ay))  # posterior probability
  
  return(ay)
})

ham <- function(x, h, J){

  m <- NROW(h)
  
  e <- 0
  for(i in seq_len(m)){
    if(x[i]==0) next()
    e <- e + h[i,x[i]]
    for(j in seq_len(m)){
      if(j==i | x[j]==0) next()
      e <- e + J[[i]][[j]][x[i],x[j]]
    }
  }
  return(e)  
}

# Pseudo partition function using data xi
Zeff <- function(L, m, h, J, xi){
  
  nsample <- NCOL(xi)
  lz <- 0
  for(k in seq_len(nsample)) for(i in seq_len(m)){
    z <- 1
    for(a in seq(1, L-1)){
      e <- h[i,a]
      for(j in seq(1,m)){
        if(j==i) next()
        xj <- xi[k,j]
        if(xj==0) next()
        e <- e + J[[i]][[j]][a,xj]
      }
      z <- z + exp(e)
    }
    lz <- lz + log(z)
  }
  lz <- lz/(nsample*m)
  
  return(lz)
}
#peff.resize(L);
#double z=1;
#for(int a=0;a<L;a++){
#  double e=h1[a];
#  for(int j=0;j< nsnp;j++){
#    if(j==i0) continue;
#    int b=ci[j];
#    if(b==0) continue;
#    e+=J1[j][L*a+b-1];
#  }
#  peff[a]=exp(e);
#  z+=peff[a];
#}
#for(int a=0;a<L;a++)
#  peff[a]/=z;

