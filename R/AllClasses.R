#' Class \code{bbm} Bayesian Boltzmann model)
#' @slot predictors If a vector of characters, finite number of factors; 
#'                  if \code{'numeric'}, multiplicative model
#' @export bbm
#' @useDynLib bbm
#' @importFrom Rcpp evalCpp
bbm <- setClass('bbm',
               slots=c(predictors='character',  # predictor factor levels
                       groups='character',  # response factor levels
                       nsite='numeric', # no of sites
                       data='data.frame',
                       h='list',      # field
                       J='list'       # coupling
#                      lz='numeric'   # log partition function
))
#' @export
setMethod('initialize', signature=('bbm'),
          definition=function(.Object, predictors, groups, nsite, data=data.frame(), ...){
            if(!is.character(predictors)) 
              predictors <- as.character(predictors)
            predictors <- unique(predictors)
            if(!is.character(groups))
              groups <- as.character(groups)
            groups <- unique(groups)
            .Object <- callNextMethod(.Object, predictors=predictors, 
                                      groups=groups, nsite=nsite, data=data, ...)
            return(.Object)
          })
#' @export
setMethod('show', signature='bbm',
          definition=function(object){
            cat('An object of class ', class(object),'\n', sep='')
            if(length(predictors)!=1 | predictors[1]!='numeric'){
              cat(' Predictor state c(',sep='')
              L <- length(object@predictors)
              for(i in seq_len(L)){ 
                cat(object@predictors[i],sep='')
                if(i < L) cat(',',sep='')
              }
              cat(') on ',sep='')
            } else cat('Numeric predictors on ',sep='')
            cat(object@nsite, ' sites\n', sep='')
            cat(' Targe state c(',sep='')
            Ly <- length(object@groups)
            for(i in seq_len(Ly)){ 
              cat(object@groups[i],sep='')
              if(i < Ly) cat(',',sep='')
            }
            cat(')\n',sep='')
            if(!is.null(object@data))
              cat(' Sample size ',NROW(object@data),'\n',sep='')
          })
