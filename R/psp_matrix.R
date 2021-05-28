#' PspMatrix 
#' @param dims dimension of matrix to construct
#' @param i row indices of non-zero entries.
#' @param j column indices of non-zero entries.
#' @param x non-zero entries corresponding to positions (\code{i}, \code{j}) in the matrix.
#' @param type storage type for the non-zero entries. Deduced from \code{x} if missing. See details. 
#' @import R6 Matrix
#' @export
PspMatrix <- R6::R6Class("PspMatrix", list(
  matrix = NULL, 
  initialize = function(i=NULL, j=NULL, x=NULL, dims=NULL, type=c("bool", "integer", "numeric")){
		ijx_supplied <- all(!c(missing(i), missing(j), missing(x)))
		matrix_supplied <- FALSE
  	if (!missing(x) && !is.null(dim(x)) && missing(i) && missing(j)){
			nz_idx <- which(x != 0, arr.ind = TRUE)
			dims <- dim(x)
			i <- nz_idx[,1]
			j <- nz_idx[,2]
			x <- x[nz_idx]
			matrix_supplied <- TRUE
		}
  	m <- new(dart:::PspBoolMatrix, dims[1], dims[2])
		if (ijx_supplied || matrix_supplied){
			stopifnot(is.vector(i), is.vector(j), is.numeric(x))
			m$construct(i-1L,j-1L,x)
		} 
		self$matrix <- m
  }, 
  print = function(){
  	if (self$matrix$nnz == 0){ cat(sprintf("< empty %d x %d PSP matrix >\n", self$matrix$n_rows, self$matrix$n_cols)) }
	  else {
	    cat(sprintf("%d x %d Permutable Sparse Matrix with %d non-zero entries", self$matrix$n_rows, self$matrix$n_cols, self$matrix$nnz))
	  }
  }
))

PspMatrix$set("public", "as.Matrix", function(type="CSC") {
  return(self$matrix$as.Matrix())
})


#' @export
`[.PspMatrix` <- function(x, i=NULL, j=NULL) { 
	if (missing(i) && missing(j)){ return(x$as.Matrix()) }
	if (missing(i) && !missing(j)){
		stopifnot(is.vector(j), is.numeric(j))
		j <- as.integer(j)
		return(x$matrix$submatrix(0L, x$matrix$n_rows-1L, max(j)-1L, max(j)-1L))
	}
	if (!missing(i) && missing(j)){
		stopifnot(is.vector(i), is.numeric(i))
		i <- as.integer(i)
		return(x$matrix$submatrix(min(i)-1L, min(i)-1L, 0L, x$matrix$n_cols-1L))
	}
	stopifnot(is.vector(i), is.numeric(i))
	i <- as.integer(i)
	j <- as.integer(j)
	return(x$matrix$submatrix(min(i)-1L, max(i)-1L, min(j)-1L, max(j)-1L))
}


#' psp_matrix 
#' @description User friendly constructed of a Permutable SParse matrix (psp_matrix). \cr
#' \cr
#' A psp_matrix is a sparse matrix representation designed to allow efficient permutation. Like any traditional 
#' sparse matrix representation (CSC, CSR, etc.), the storage complexity of a psp_matrix is proportional to the 
#' number of non-zero entries. However, swapping any two rows or columns takes just constant time. With other 
#' representations, performing a row/column swap in a matrix A can require up to O(nnz(A)) where 
#' nnz(A) = number of non-zero entries in A.
#' @import Matrix 
#' @export
psp_matrix <- function(i=NULL, j=NULL, x=NULL, dims=NULL, type=c("bool", "integer", "numeric")){
	m <- do.call(PspMatrix$new, as.list(match.call())[-1])
	return(m)
}

methods::setAs(from = "PspMatrix", to = "sparseMatrix", def = function(from){
	return(from$as.Matrix())
})
methods::setAs(from = "PspMatrix", to = "dtCMatrix", def = function(from){
	m <- from$as.Matrix()
	tri <- Matrix::sparseMatrix(i=m@i, p = m@p, x = as.integer(m@x), dims = m@Dim, triangular = TRUE, index1 = FALSE)
	return(tri)
})

# coerce.PspMatrix <- function(from, to){
# 	
# }



# setClass("Rcpp_PspBoolMatrix")
# .print_psp_matrix <- setMethod("show", "Rcpp_PspBoolMatrix", function (object) {
#   if (object$nnz == 0){ cat(sprintf("< empty %d x %d PSP matrix >\n", object$n_rows, object$n_cols)) }
#   else {
#     cat(sprintf("%d x %d Permutable Sparse Matrix with %d non-zero entries", object$n_rows, object$n_cols, object$nnz))
#   }
# })


format.Rcpp_PspBoolMatrix <- function(x, ...){
	print("hello world")
}