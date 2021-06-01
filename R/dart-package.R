#' dart: a package containing tools for simulating persistent homology dynamically.
#' @docType package
#' @name dart
#' @import methods Rcpp
#' @useDynLib dart, .registration = TRUE
#' @keywords internal
"_PACKAGE"


#' simplex_to_str
#' @export
simplex_to_str <- function(x){
	if (is.null(x) || length(x) == 0){ return(NULL) }
	to_str <- function(simplex){ sprintf("(%s)",paste0(simplex, collapse=",")) }
	if (!is.null(dim(x))){
		return(simplex_to_str_rcpp(split(x, slice.index(x,2))))
	} else if (is.numeric(x) && is.vector(x)){
		return(sprintf("(%s)",paste0(x, collapse=",")))
	} else if (is.list(x)){
		return(simplex_to_str_rcpp(x))
	} else { stop("Invalid format given.") }
}

#' str_to_simplex
#' @param x a vector, list, or other iterable of strings to give to sapply. 
#' @param ... passed to sapply. 
#' @export
str_to_simplex <- function(x, ...){
	if (is.null(x) || length(x) == 0){ return(NULL) }
	to_simplex <- function(s){ as.integer(strsplit(substring(s, first = 2L, last = nchar(s)-1L), split=",")[[1]]) }
	sapply(x, to_simplex, ...)
}


Rcpp::loadModule("PspBoolMatrix", TRUE)
Rcpp::loadModule("implicit_filtration_module", TRUE)
# Rcpp::loadModule("A_printable", TRUE)
# loadModule("PspIntMatrix", TRUE)
