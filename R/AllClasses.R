#' @useDynLib BBM
#' @importFrom Rcpp evalCpp
#' @export bbm
bbm <- setClass('bbm',
               slots=c(predictors='character',  # predictor factor levels
                       groups='character',  # response factor levels
                       nsite='numeric', # no of sites
                       data='data.frame',
                       h='list',      # field
                       J='list'      # coupling
))
#' @export
setMethod('initialize', signature=('bbm'),
          definition=function(.Object, predictors, groups, nsite, 
                              xi=NULL, y=NULL, ...){
            if(!is.character(predictors)) 
              predictors <- as.character(predictors)
            predictors <- unique(predictors)
            if(!is.character(groups))
              groups <- as.character(groups)
            groups <- unique(groups)
            .Object <- callNextMethod(.Object, predictors=predictors, 
                                      groups=groups, nsite=nsite, ...)
            if(!is.null(xi)) .Object@xi <- xi
            if(!is.null(y)) .Object@y <- y
            return(.Object)
          })
#' @export
setMethod('show', signature='bbm',
          definition=function(object){
            cat('An object of class ', class(object),'\n', sep='')
            cat(' Predictor state c(',sep='')
            L <- length(object@predictors)
            Ly <- length(object@groups)
            for(i in seq_len(L)){ 
              cat(object@predictors[i],sep='')
              if(i < L) cat(',',sep='')
            }
            cat(') on ',object@nsite, ' sites\n', sep='')
            cat(' Targe state c(',sep='')
            for(i in seq_len(Ly)){ 
              cat(object@groups[i],sep='')
              if(i < Ly) cat(',',sep='')
            }
            cat(')\n',sep='')
          })
