#' Convert Frequency Table into Raw Data
#' 
#' Data with unique rows and a frequency column is converted into
#' data with duplicate rows.
#' 
#' The ouput data frame can be used as input to \code{\link{bbl}}.
#' 
#' @param fdata Data frame with factors in columns and a frequency column
#' @param Freq Integer column index or column name of frequency
#' @return Raw data frame with one row per instances
#' @examples
#' Titanic
#' x <- as.data.frame(Titanic)
#' head(x)
#' titanic <- freq2raw(x, Freq='Freq')
#' head(titanic)
#' @export
freq2raw <- function(fdata,Freq){
  
  if(is.character(Freq)) 
    Freq <- which(colnames(fdata)==Freq)
  if(length(Freq)==0) stop('No frequency column found')
  n <- nrow(fdata)
  dat <- NULL
  for(i in 1:n){
    w <- fdata[rep(i,fdata[i,Freq]),]
    dat <- rbind(dat, w)
  }
  rownames(dat) <- seq_len(NROW(dat))
  return(dat[,-Freq])
}

#' Compute Prediction Accuracy
#' 
#' Accuracy of predicted response probability is computed.
#' 
#' An option is provided for computing group-balanced accuracy, where
#' prediction score is calculated for each group separately and averaged.
#' @param object Object of class \code{bbl} with test data in \code{data} slot.
#' @param prediction Data frame of predicted response group probability from
#'        \code{\link{predict}}.
#' @param balanced Compute balanced accuracy. If \code{TRUE},
#'        \deqn{s = \frac{1}{K}\sum_y \frac{1}{n_y} \sum_{k\in y}
#'        \delta\left({\hat y}_k = y\right).}
#'        If \code{FALSE},
#'        \deqn{s = \frac{1}{n}\sum_{k}
#'        \delta\left({\hat y}_k = y_k\right).}
#' @return List of \code{acc} (accuracy score) and \code{yhat} (predicted
#'        response group). 
#' @examples
#' titanic <- freq2raw(as.data.frame(Titanic), Freq='Freq')
#' nsample <- NROW(titanic)
#' mod <- bbl(data=titanic, y='Survived')
#' mod <- mod[sample(nsample),]
#' mtrain <- mod[seq(nsample/2),]
#' mtest <- mod[seq(nsample/2,nsample),]
#' mtrain <- train(mtrain, method='mf')
#' pred <- predict(mtrain, newdata=mtest@data)
#' score <- accuracy(mtest, prediction=pred, balanced=TRUE) 
#' @export
accuracy <- function(object, prediction, balanced=FALSE){
  
  groups <- object@groups
  
  y <- object@data[,colnames(object@data)==object@y]
  yhat <- colnames(prediction)[apply(prediction, 1, which.max)]
  
  if(balanced){
    acav <- NULL
    for(k in seq_along(groups)){
      gr <- groups[k]
      sub <- y==gr
      acc <- mean(yhat[sub]==gr)
      acav <- c(acav,acc)
    }
    acav <- mean(acav)
  }
  else
    acav <- mean(y==yhat)
  
  return(list(acc=acav, yhat=yhat))
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
read.fasta <- function(file, rownames=FALSE){
  
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