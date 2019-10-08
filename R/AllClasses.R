#' Class \code{bbl} for Boltzmann Bayes Learner
#' @slot predictors List of vectors of characters, each corresponding
#'   to predictor factor levels.
#' @slot groups Vector of characters for response factor levels.
#' @slot data Data frame for training data. Expected to contain both
#' predictors and response.
#' @slot y Character column name of response data in \code{data}.
#' @slot h Bias parameters stored after training. List of length
#'   equal to no. of response groups, each of which is a list of
#'   length equal to no. of predictors, containing vectors of 
#'   length equal to each predictor factor levels:
#'   \eqn{h_i^{(y)}(x)} represented by
#'   \code{h[[y]][[i]][x]}.
#' @slot J Interaction parameters stored after training. List of
#'   length equal to no. of response groups, each of which is a 
#'   list of lists of dimension \eqn{m \times m}, where \eqn{m}
#'   is the number of predictors. Each element is a matrix of
#'   dimension \eqn{L_i \times L_j}, where \eqn{L_i} and \eqn{L_j}
#'   are numbers of factor levels in predictor \code{i} and \code{j}:
#'   \eqn{J_{ij}^{(y)}(x_1,x_2)} represented by
#'   \code{J[[y]][[i]][[j]][x1,x2]}.
#' @slot lz Slot used to store log partition function of response
#'   groups.
#' @name bbl-class
#' @export bbl
#' @useDynLib bbl
#' @importFrom Rcpp evalCpp
bbl <- setClass('bbl',
               slots=c(predictors='list',  # predictor factor levels
                       groups='character', # response factor levels
                       data='data.frame',
                       y='character',      # target variable
                       h='list',        # field
                       J='list',        # coupling
                       lz='numeric'     # log partition function
))

#' Instantiate \code{bbl} object
#' 
#' Creates \code{bbl} object based on arguments provided.
#' 
#' @param .Object Object of class \code{bbl} to be created.
#' @param data Data frame of training data. Expected to contain both 
#'        predictor and response variables in columns and instances 
#'        in rows.
#' @param predictors List of predictor factor levels. If not provided, 
#'        will be inferred from \code{data} columns with names not equal 
#'        to \code{y}.
#' @param groups Vector of characters for response factor levels. If not
#'        provided, will be inferred from \code{data} response column
#'        named \code{y}.
#' @param y Column name in \code{data} containing response data. Can also
#'        be column number.
#' @param ... Other parameters (not used).
#' @return  Object of class \code{bbl}.
#' @details Note that with argument \code{data} only, \code{predictor}
#'        factor levels will be ordered by default. Also some of predictors 
#'        may have values in \code{data} less than expected, which would 
#'        result in fewer levels of \code{predictor} levels. By explicitly 
#'        providing \code{predictors}, this default behavior can be overriden. 
#'        Predictors in \code{data} with a single factor level will be removed.
#'        Same applies to \code{groups} for response levels. 
#' @export
#' @examples
#' ## Factors model
#' titanic <- freq2raw(as.data.frame(Titanic), Freq='Freq')
#' summary(titanic)
#' model <- bbl(data=titanic, y='Survived')
#' model
#' @importFrom methods initialize new callNextMethod is
setMethod('initialize', signature=('bbl'),
          definition=function(.Object, data, predictors=NULL, groups=NULL, 
                              y='y', ...){
            if(!is.data.frame(data))
              data <- as.data.frame(data)
            if(!y %in% colnames(data)) stop(paste0(y,' not in data'))
            iy <- which(colnames(data)==y)
            if(is.null(groups)) 
              groups <- levels(factor(data[,iy]))
            if(length(groups)==1) stop('Only one response group present')
            xi <- data[,-iy]
            nvar <- ncol(xi)
            if(is.null(predictors)){
              predictors <- list()
              cnt <- NULL
              for(i in seq_len(nvar)){
                fac <- levels(factor(xi[,i]))
                fac <- fac[order(fac)]
                if(length(fac)==1) next()
                cnt <- c(cnt, i)
                predictors[[length(cnt)]] <- fac
              }
            } else{
              cnt <- NULL
              for(i in seq_len(nvar)){
                fac <- levels(factor(xi[,i]))
                if(sum(!fac%in%predictors[[i]])>0)
                  stop('Input data contains factor levels not in predictor')
                cnt <- c(cnt, i)
              } 
            }
            if(length(cnt)==0) stop('No variable with multiple levels')
            data <- cbind(data.frame(y=data[,iy]), xi[,cnt])
            colnames(data)[1] <- y
            if(length(cnt)<nvar)
              warning(paste0(' ',nvar-length(cnt),
                      ' variables with one level removed\n'))
            .Object <- callNextMethod(.Object, predictors=predictors,
                                      groups=groups, data=data, y=y, ...)
            return(.Object)
})

#' @describeIn bbl
#' Display \code{bbl} object content
#' 
#' Show predictor levels, response groups, and sample size.
#' @param object Object of class \code{bbl} to display
#' @examples
#' data <- data.frame(y=c('g1','g2','g2','g1','g1','g2'),
#'   x1=c(0,0,1,2,0,1), x2=c('a','a','b','a','b','b'))
#' x <- bbl(data)
#' x
#' x[1:3,]
#' @importFrom methods show
#' @export
setMethod('show', signature='bbl',
          definition=function(object){
            cat('An object of class ', class(object),'\n', sep='')
            iy <- which(colnames(object@data)==object@y)
            xi <- object@data[,-iy]
            nvar <- ncol(xi)
            predictors <- object@predictors
            cat(' ',nvar, ' predictor states:\n',sep='')
            for(i in seq_len(min(3,length(predictors))))
              cat('  ',colnames(xi)[i],'=',predictors[[i]],'\n',sep=' ')
            if(nvar > 3) cat('   ...\n',sep='')
            cat(' Responses:\n  ', object@y,'=',object@groups, '\n',sep=' ')
            cat(' Sample size: ',NROW(object@data),'\n',sep='')
            return(invisible(object))
          })
