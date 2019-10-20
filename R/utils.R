#' Convert Frequency Table into Raw Data
#' 
#' Data with unique rows and a frequency column is converted into
#' data with duplicate rows.
#' 
#' The ouput data frame can be used as input to \code{\link{bbl}}.
#' 
#' @param fdata Data frame with factors in columns 
#' @param Freq Vector of frequency of each row in \code{fdata}
#' @return Raw data frame with one row per instances
#' @examples
#' Titanic
#' x <- as.data.frame(Titanic)
#' head(x)
#' titanic <- freq2raw(fdata=x[,1:3], freq=x$Freq)
#' head(titanic)
#' @export
freq2raw <- function(fdata,freq){
  
  if(length(freq)!=NROW(fdata)) 
     stop('Frequency length does not match fdata')
  n <- nrow(fdata)
  dat <- NULL
  for(i in 1:n){
    w <- fdata[rep(i,freq[i]),]
    dat <- rbind(dat, w)
  }
  rownames(dat) <- seq_len(NROW(dat))
  return(dat)
}


#' Read FASTA file
#' 
#' Read nucleotide sequence files in FASTA format
#' 
#' Sequence data in FASTA files are converted into data frame
#' suitable as input to \code{\link{bbl}}. If sequence lengths are different,
#' instances longer than those already read will be truncated. Empty sequences
#' are skipped.
#' 
#' @param file File name of FASTA input.
#' @param rownames Use the sequence annotation line in file (starts with
#'        \code{'>'}) as the row names. Will fail if there are duplicate items.
#' @return Data frame of each sequence in rows.
#' @examples
#' file <- tempfile('data')
#' write('>seq1', file)
#' write('atgcc', file, append=TRUE)
#' write('>seq2', file, append=TRUE)
#' write('gccaa', file, append=TRUE)
#' system(paste0('cat ',file))
#' x <- read.fasta(file)
#' x
#' @export
readFasta <- function(file, rownames=FALSE){
  
  if(!file.exists(file)) stop(paste0(file,' does not exist'))
  fl <- readLines(file)
  if(length(fl)==0) stop(paste0(file, ' is empty'))
  label <- dat <- NULL
  
  for(i in seq_along(fl)){
    if(substr(fl[[i]],start=1,stop=1)=='>'){
      x <- strsplit(fl[[i]],split='\t')[[1]]
      x <- paste(x,collapse=' ')
      label <- c(label, x)
    }
    else{
      seqs <- strsplit(fl[[i]],split='')[[1]]
      if(NROW(dat)>0){ 
        if(length(seqs)<NCOL(dat)) next()
        if(length(seqs)>NCOL(dat)) seqs <- seqs[seq_len(NCOL(dat))]
      }
      dat <- rbind(dat, as.data.frame(matrix(seqs,nrow=1)))
    }
  }
  if(rownames){
    if(length(label)!=NROW(dat)) 
      stop('Annotation lines do not match sequences.')
    rownames(dat) <- label
  }
  colnames(dat) <- seq_len(NCOL(dat))
  
  return(dat)
}

#' Remove non-varying predictors
#' @export
removeConst <- function(x){
  
  bad <- lapply(x, function(x){length(levels(factor(x)))==1})
  bad <- unlist(bad)
  return(x[,!bad])
}

#' Newton-Raphson for f(h)=xav or h=f^-1(xav)
#' 
fu <- function(x, L, xav){
  
  L/(1-exp(-L*x)) - 1/(1-exp(-x)) - xav
  
}
#' 1st derivative
dfu <- function(x, L, xav){
  
  -L^2*exp(-L*x)/(1-exp(-L*x))^2 + exp(-x)/(1-exp(-x))^2
  
}

#' Root finding function
nr <- function(fu, dfu, xinit=0.1, tol=1e-5, maxit=100, ...){
  
  x <- xinit
  for(i in seq_len(maxit)){
    xp <- x - fu(x, ...)/dfu(x, ...)
    df <- abs((xp/x)^2-1)
    if(df < tol) break()
    x <- xp
  }
  if(i >= maxit) warning('Maximim iteration limit reached')
  return(x)
}