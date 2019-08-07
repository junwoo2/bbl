#' @useDynLib BBM
#' @importFrom Rcpp evalCpp
#' @export bbm
bbm <- setClass('bbm',
               slots=c(L='numeric',  # number of factor levels in X
                       Ly='numeric',  # no. of levels in target var. Y
                       nsite='numeric', # no of sites
                       xi='matrix',
                       y='vector',
                       h='list',      # field
                       J='list'      # coupling
))
#' @export
setMethod('initialize', signature=('bbm'),
          definition=function(.Object, L=2, Ly=2, nsite, xi=NULL, y=NULL,
                              ...){
            .Object <- callNextMethod(.Object, L=L, Ly=Ly, nsite=nsite, ...)
            if(!is.null(xi)) .Object@xi <- xi
            if(!is.null(y)) .Object@y <- y
            return(.Object)
          })
#' @export
setMethod('show', signature='bbm',
          definition=function(object){
            cat('An object of class ', class(object),'\n', sep='')
            cat(' ', object@L, ' states on ',object@nsite, ' sites\n', sep='')
            cat(' ', object@Ly,' target states\n',sep='')
          })
