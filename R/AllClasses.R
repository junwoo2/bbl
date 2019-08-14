#' Class \code{bbm} Bayesian Boltzmann model)
#' @slot predictors If a vector of characters, finite number of factors; 
#'                  if \code{'numeric'}, multiplicative model
#' @export bbm
#' @useDynLib bbm
#' @importFrom Rcpp evalCpp
bbm <- setClass('bbm',
               slots=c(type='character',   # factor or numeric
                       predictors='list',  # predictor factor levels
                       groups='character', # response factor levels
                       data='data.frame',
                       y='character',      # target variable
                       h='list',        # field
                       J='list',        # coupling
                       lz='numeric'     # log partition function
))
#' @export
setMethod('initialize', signature=('bbm'),
          definition=function(.Object, type='factors',predictors=NULL, groups=NULL,
                              data, y='y', ...){
            if(!is.data.frame(data))
              data <- as.data.frame(data)
            if(!y %in% colnames(data)) stop(paste0(y,' not in data'))
            iy <- which(colnames(data)==y)
            if(is.null(groups)) groups <- levels(factor(data[,iy]))
            xi <- data[,-iy]
            nvar <- ncol(xi)
            if(type=='numeric') predictors <- list('numeric')
            else predictors <- list()
            cnt <- NULL
            for(i in seq_len(nvar)){
              fac <- levels(factor(xi[,i]))
              if(length(fac)==1) next()
              cnt <- c(cnt, i)
              if(type=='factors') predictors[[length(cnt)]] <- fac
            }
            if(length(cnt)<nvar){
              cat(' ',nvar-length(cnt),' variables with one level removed\n',
                    sep='')
              nvar <- length(cnt)
              data <- cbind(xi[,cnt], data.frame(y=data[,iy]))
            }
            .Object <- callNextMethod(.Object, type=type, predictors=predictors, 
                                      groups=groups, data=data, y=y, ...)
            return(.Object)
          })
#' @export
setMethod('show', signature='bbm',
          definition=function(object){
            cat('An object of class ', class(object),'\n', sep='')
            iy <- which(colnames(object@data)==object@y)
            xi <- object@data[,-iy]
            nvar <- ncol(xi)
            if(object@type=='factors'){
              predictors <- object@predictors
              cat(' ',nvar, ' predictor states:\n',sep='')
              for(i in seq_len(min(3,length(predictors))))
                cat('  ',colnames(xi)[i],'=',predictors[[i]],'\n',sep=' ')
              if(nvar > 3) cat('   ...\n',sep='')
            }
            else{ 
              cat(' ',nvar,' numeric predictors:\n',sep='')
              for(i in seq_len(min(3, nvar)))
                cat('   [',min(xi[,i]),',',max(xi[,i]),']\n',sep='')
              if(nvar > 3) cat('   ...\n',sep='')
            }
            cat(' Responses:\n  ', object@y,'=',object@groups, '\n',sep=' ')
            cat(' Sample size: ',NROW(object@data),'\n',sep='')
          })
