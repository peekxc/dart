## reduction.R

#' Reduction algorithm 
#' @param D either filtration boundary matrix or an object of class "fbm" from e.g. \link[dart]{boundary_matrix}.
#' @param field Field of coefficients to use. Defaults to modulo 2 coefficients. See details.
#' @param output Type of reduction output requested. See details. 
#' @param indices Integer indices of which columns to reduce. By default all columns are reduced. 
#' @param options Options and optimizations which affect how the reduction proceeds. 
#' @param show_progress whether to print a progress bar as the columns are reduced. 
#' @description The standard PH algorithm (and its variants) 
#' @details This function applies the standard algorithm to decompose a filtration boundary matrix representation
#' into a variety of outputs. 
#' @import RcppProgress
#' @export
reduce <- function(D, field=c("mod2", "reals"), output = c("RV", "RU", "R", "V", "U"), 
									 indices = NULL,
									 options = c(clearing=TRUE), 
									 validate = FALSE, show_progress=FALSE) {
	if (missing(field) || field == "mod2"){
		coerce_field <- function(x){ x %% 2L }
	} else {
		coerce_field <- identity
	}
	use_clearing <- options[["clearing"]]
	if ("fbm" %in% class(D)){
		if (is.list(D$matrix)){
			k <- length(D$matrix)
			stopifnot(k >= 2)
			R_psp <- lapply(seq(k), function(i){ dart::psp_matrix(x = D$matrix[[i]]) })
			V_psp <- lapply(seq(k), function(i){ dart::psp_matrix(x = Matrix::Diagonal(ncol(D$matrix[[i]]))) })
			
			## Perform the reductions on each of the matrices
			nc <- 0L
			for (j in seq(k, 2L)){
				i <- j - 1L
				nc <- nc + reduce_local_pspbool(
					D1_ptr = R_psp[[i]]$matrix$as_XPtr(), V1_ptr = V_psp[[i]]$matrix$as_XPtr(), 
					D2_ptr = R_psp[[j]]$matrix$as_XPtr(), V2_ptr = V_psp[[j]]$matrix$as_XPtr(), 
					clearing = use_clearing, show_progress = show_progress
				)
			}
			res <- list(R = R_psp, V = V_psp, nc = nc)
			# D1_psp <- dart::psp_matrix(x = D$matrix[[1]])
			# V1_psp <- dart::psp_matrix(x = Matrix::Diagonal(ncol(D$matrix[[1]])))
			# D2_psp <- dart::psp_matrix(x = D$matrix[[2]])
			# V2_psp <- dart::psp_matrix(x = Matrix::Diagonal(ncol(D$matrix[[2]])))
			# 
			# 
			# 
			# 
			# 
			# ## Perform the reductions on each of the matrices
			# nc <- reduce_local_pspbool(
			# 	D1_ptr = D1_psp$matrix$as_XPtr(), V1_ptr = V1_psp$matrix$as_XPtr(), 
			# 	D2_ptr = D2_psp$matrix$as_XPtr(), V2_ptr = V2_psp$matrix$as_XPtr(), 
			# 	clearing = use_clearing, show_progress = show_progress
			# )
			# 
			# ## Clean entries up 
			# D1_psp$matrix$clean(0)
			# D2_psp$matrix$clean(0)
			# V1_psp$matrix$clean(0)
			# V2_psp$matrix$clean(0)
			
			# if (missing(output) || output == "RV"){
			# 	res <- list(R=list(D1_psp, D2_psp), V=list(V1_psp, V2_psp))
			# } else if (out == "RU" || "U"){
			# 	U1 <- coerce_field(solve(as(V1_psp, "sparseMatrix")))
			# 	U2 <- coerce_field(solve(as(V2_psp, "sparseMatrix")))
			# 	res <- ifelse(out == "RU", 
			# 				 list(R=list(D1_psp, D2_psp), U=list(psp_matrix(x=U1), psp_matrix(x=U2))), 
			# 				 list(U=list(psp_matrix(x=U1), psp_matrix(x=U2))))
			# } else if (out == "R"){
			# 	res <- list(R=list(D1_psp, D2_psp))
			# } else if (out == "V"){
			# 	res <- list(V=list(V1_psp, V2_psp))
			# }
			# res$nc <- nc
		} else {
			stopifnot(!is.null(dim(D$matrix)))
			
			## Perform the reduction on the whole matrix 
			D_psp <- dart::psp_matrix(x = D$matrix)
			V_psp <- dart::psp_matrix(x = Matrix::Diagonal(ncol(D$matrix)))
			nc <- reduce_pspbool(D_psp$matrix$as_XPtr(), V_psp$matrix$as_XPtr(), show_progress=show_progress)
			
			## Determine the output type
			if (missing(output) || output == "RV"){
				res <- list(R=D_psp, V=V_psp)
			} else if (out == "RU" || "U"){
				U <- solve(as(V_psp, "sparseMatrix"))
				res <- ifelse(out == "RU", list(R=D_psp, U=psp_matrix(x=U)), list(U=psp_matrix(x=U)))
			} else if (out == "R"){
				res <- list(R=D_psp)
			} else if (out == "V"){
				res <- list(V=V_psp)
			}
			res$nc <- nc
		}
	} else {
		stopifnot(!is.null(dim(D)))
		## Perform the reduction on the whole matrix 
		D_psp <- dart::psp_matrix(x = D)
		V_psp <- dart::psp_matrix(x = Matrix::Diagonal(ncol(D)))
		nc <- reduce_pspbool(D_ptr = D_psp$matrix$as_XPtr(), V_ptr = V_psp$matrix$as_XPtr(), 
									 show_progress = show_progress)
		
		## Determine the output type
		if (missing(output) || output == "RV"){
			res <- list(R=D_psp, V=V_psp)
		} else if (out == "RU" || "U"){
			U <- solve(as(V_psp, "sparseMatrix"))
			res <- ifelse(out == "RU", list(R=D_psp, U=psp_matrix(x=U)), list(U=psp_matrix(x=U)))
		} else if (out == "R"){
			res <- list(R=D_psp)
		} else if (out == "V"){
			res <- list(V=V_psp)
		}
		res$nc <- nc
	}
	class(res) <- "fbm_decomp"
	if (validate){
		if ("fbm" %in% class(D) && is.list(D$matrix)){
			valid <- validate_decomp(res, D)
			if (any(!valid)){
				stop(sprintf("Invalid: R reduced = %s, V triangular = %s, D := RV = %s", 
										 valid["reduced"], valid["triangular"], valid["decomposition"]))
			}
		}
	}
	attr(res, "hom_dim") <- attr(D, "hom_dim")
	# output_type <- match(output, c("RV", "RU", "R", "V", "U"))
	return(res)
}

#' @method print fbm_decomp
#' @export print.fbm_decomp 
#' @export
print.fbm_decomp <- function(x, ...){
	stopifnot("fbm_decomp" %in% class(x))
	# if (is.list(x$matrix)){
	# 	num_bm <- length(x$matrix)
	# 	bm_str <- do.call(sprintf, append(list(paste(rep("D%d", num_bm), collapse = ", ")), attr(x, "hom_dim")))
	# 	bm_str <- paste0(bm_str, sprintf("-filtration boundary %s", ifelse(num_bm == 1, "matrix", "matrices")))
	# 	sz_str <- paste(sapply(x$matrix, function(bm) do.call(sprintf, append(list("%d x %d"), dim(bm)))), collapse = ", ")
	# 	sz_str <- paste(ifelse(num_bm == 1, "Size:", "Sizes:"), sz_str)
	# } else {
	# 	bm_str <- sprintf("Filtration boundary matrix (dim=%s)", paste0(as.character(attr(x, "hom_dim")), collapse = ","))
	# 	sz_str <- sprintf("Size: %d x %d", nrow(x$matrix), ncol(x$matrix))
	# }
	# writeLines(c(bm_str, sz_str))
	cat("R = DV Decomposition")
}

#' @export
is_reduced <- function(x){
	P <- pivots(x)
	is_reduced <- !as.logical(anyDuplicated(P[,1]))
	return(is_reduced)
}


#' Pivot entries 
#' @description Returns the pivot entries in a matrix.
#' @param x a matrix type. 
#' @details The function returns the (row, column) index pairs of the lowest row entries associated 
#' to each non-zero column of \code{x}, with optimizations used if \code{x} is a sparse matrix. 
#' @export
pivots <- function(x){
	if (inherits(x, "PspMatrix")){
		return(t(x$matrix$pivots()+1L))
		# ci <- seq(x$matrix$n_cols)
		# ri <- sapply(ci-1L, function(j){ x$matrix$low_entry(j)+1 })
		# return(cbind(row=ri[ri != 0], col=ci[ri != 0]))
	} else if (as.logical(length(grep("*Matrix", class(x))))){
		if (all(c("i", "p") %in% slotNames(x))){
			cs <- diff(cumsum(x@p))
			ri <- x@i[unique(cs)]+1L
			ci <- which(!duplicated(cs))
			return(cbind(row=ri,col=ci))
		} else if ("ddiMatrix" %in% class(x)){
			stopifnot(ncol(x) == nrow(x))
			return(cbind(row=seq(nrow(x)), col=seq(nrow(x))))
		} else {
			stop("unknown matrix type supplied")
		}
	} else if (is.matrix(x)){
		ci <- seq(ncol(x))
		ri <- apply(x, 2, dart:::low_entry)
		return(cbind(row=ri[ri != 0], col=ci[ri != 0]))
	} else {
		stop("Unknown matrix type given")
	}
}

#' Validates reduction decomposition
#' @export
validate_decomp <- function(x, D = NULL, field = c("mod2", "reals")){
	stopifnot("fbm_decomp" %in% class(x))
	stopifnot(all(c("R", "V") %in% names(x)))
	if (!missing(D)){ stopifnot("fbm" %in% class(D)) }
	if (missing(field) || field == "mod2"){
		coerce_field <- function(x){ x %% 2L }
	} else {
		coerce_field <- identity
	}
	is_red <- ifelse(is.list(x$R), all(sapply(x$R, is_reduced)), is_reduced(x$R))
	is_tri <- function(v) { as.logical(Matrix::isTriangular(v)) }
	if (is.list(x$V)){
		is_upt <- all(sapply(x$V, function(v){ is_tri(as(v, "sparseMatrix")) }))
	} else {
		is_upt <- is_tri(as(x$V, "sparseMatrix"))
	}
	if (!missing(D)){
		if ("fbm" %in% class(D) && is.list(D$matrix)){
			stopifnot(length(D$matrix) == length(x$V), length(D$matrix) == length(x$R))
			# is_psp <- all(sapply(x$V, function(m) "PspMatrix" %in% class(m)))
			diff <- lapply(seq_along(D), function(i){ 
				DV <- coerce_field(as(D$matrix[[i]], "sparseMatrix") %*% as(x$V[[i]], "sparseMatrix")) 
				R <- coerce_field(as(x$R[[i]], "sparseMatrix"))
				return(coerce_field(R - DV))
			})
			is_dec <- sapply(diff, function(d) all(d == 0))
		} else if ("fbm" %in% class(D)){
			R <- coerce_field(as(x$R, "sparseMatrix"))
			B <- coerce_field(as(D$matrix, "sparseMatrix"))
			V <- coerce_field(as(x$V, "sparseMatrix"))
			diff <- coerce_field(R-(B%*%V))
			is_dec <- all(diff == 0)
		} else {
			stopifnot(!is.null(dim(D)))
			R <- coerce_field(as(x$R, "sparseMatrix"))
			B <- coerce_field(as(D, "sparseMatrix"))
			V <- coerce_field(as(x$V, "sparseMatrix"))
			diff <- coerce_field(R-(B%*%V))
			is_dec <- all(diff == 0)
		}
	}
	if (missing(D)){
		return(c(reduced=is_red, triangular=is_upt))
	} 
	return(c(reduced=is_red, triangular=is_upt, decomposition=is_dec))
}

#' Persistent homology computation 
#' @description This function offers a generic persistent homology computation using the reduction algorithm. 
#' While certain optimizations are available (see \link{reduce}), this implementation is not guaranteed to be 
#' necessarily efficient for all filtrations, though it is relatively efficient. One can specify whether to 
#' compute the persistence of multiple dimensions simultaneously or alternatively specific dimensions of homology 
#' can be chosen by setting \code{dim} accordingly. 
#' @param x point cloud matrix, 'dist' object, or Filtration.
#' @export
ph <- function(x, dims = c(0L, 1L), type = c("rips", "cech", "delaunay")){
	contiguous <- all(dims == seq(0, max(dims)))
	if (contiguous){
		fi <- dart::rips_filtration(x, dim = max(dims)+1L)
		D <- dart::boundary_matrix(fi, dims = "all")
		rv <- dart::reduce(D, validate = FALSE)
		colnames(rv$R) <- simplex_to_str(fi$simplices)
		bc <- dart::read_barcodes(rv, essential = FALSE, dims = dims, values = fi$grading, collapse = TRUE)
	} else {
		fi <- dart::rips_filtration(x, dim = max(dims)+1L)
		D <- dart::boundary_matrix(fi, dims = c(dims, max(dims)+1L))
		rv <- dart::reduce(D, validate = FALSE)
		values <- lapply(attr(D, "hom_dim"), function(di){ fi$k_grades(di) })
		bc <- dart::read_barcodes(rv, essential = FALSE, dims = dims, values = values, collapse = TRUE)
	}
	return(bc)
}
	
	
	# for (n in 3:15){
	# 	for (i in 1:100){
	# 		set.seed(i)
	# 		x1 <- replicate(2, rnorm(n))
	# 		x2 <- x1 + replicate(2, runif(nrow(x1)))
	# 		r1 <- rips_filtration(x1, diameter = Inf, dim = 2)
	# 		r2 <- rips_filtration(x2, diameter = Inf, dim = 2)
	# 		S <- vineyards_schedule(r1$simplices, r2$simplices, w0 = r1$grading, w1 = r2$grading, schedule_order = "linear")
	# 	}
	# }



#' Read bardcodes
#' @param x a reduced boundary matrix with column labels, or the output object of \link{reduce}. See details. 
#' @param dims the homology class dimensions to include in the output.
#' @param values optional function values associated with the filtration. If excluded the index persistence is reported. See details.  
#' @param collapse whether to remove persistent pairs that have 0 persistence. Defaults to true. 
#' @param essential whether to keep essential classes which do not die. Defaults to false. 
#' @description Extracts the persistence diagram from reduced matrices.
#' @details Given a reduced boundary matrix \code{x}, this function extracts the persistent pairs
#' representing points on the persistence diagram. If \code{x} is a matrix, it is expected to have simplices (as strings)
#' in its \code{colnames} to extract dimension information, otherwise if \code{x} is a list encoding each reduction matrix 
#' separately the dimension information is inferred from the supplied \code{dims}.
#' @export 
read_barcodes <- function(x, dims = c(0L, 1L), values=NULL, collapse=TRUE, essential=FALSE){
	.matrix_like <- function(m){ !is.null(dim(m)) }
	stopifnot(inherits(x, "fbm_decomp") || .matrix_like(x))
	stopifnot(is.vector(dims) && all(is.numeric(dims)))
	
	## If a decomposition was given but the elements are just matrices, 
	## reduce to the full matrix case
	if (inherits(x, "fbm_decomp")){
		R <- x$R
	} else {
		R <- x
	}
	
	## If values supplied, then check that: 
	##	1) they are a numeric vectors matching the number columns of R, if R is a matrix, or
	##	2) they are a list of vectors, each matching the number of columns of each submatrix R
	if (!missing(values)){ 
		if (.matrix_like(R)){
			stopifnot(is.numeric(values) && length(values) == ncol(R)) 
		} else {
			stopifnot(is.list(values))
			stopifnot(all(sapply(values, length) == sapply(R, ncol)))
		}
	}
	
	## If R is a matrix, treat it like the full reduced boundary matrix 
	if (.matrix_like(R)){
		{ nr <- nrow(R); nc <- ncol(R) }
		s_dims <- sapply(str_to_simplex(colnames(R)), length) - 1L
		P <- pivots(R)
		
		## The pivots define the non-essential persistence pairs
		barcodes <- cbind(s_dims[P[,1]], P)
		colnames(barcodes) <- c("dimension", "birth", "death")
		
		## If requested, record essential classes which never die 
		if (essential){
			E <- setdiff(seq.int(nc), as.vector(P))
			barcodes <- rbind(barcodes, cbind(s_dims[E], E, Inf))
		}
		
		## Fix order to the lexicographical 
		barcodes <- barcodes[order(barcodes[,1], barcodes[,2], barcodes[,3]),]
		
		## If weights supplied, use report them
		if (!missing(values) && !is.null(values)){
			barcodes[,2] <- values[barcodes[,2]]
			barcodes[barcodes[,3] != Inf, 3] <- values[barcodes[barcodes[,3] != Inf, 3]]	
		}
		
		## If requested, remove 0-persistence pairs 
		if (collapse){
			persistence <- abs(barcodes[,2] - barcodes[,3])
			barcodes <- barcodes[persistence > .Machine$double.eps,,drop=FALSE]
		}
	} else if (is.list(R) && inherits(x, "fbm_decomp")) {
		s_dims <- attr(x, "hom_dim")
		stopifnot(is.numeric(s_dims), length(s_dims) == length(R), length(R) >= 2L)
		if (!missing(values)){ stopifnot(is.list(values), all(sapply(values, length) == sapply(R, ncol))) }
		k <- length(R)
		barcodes <- matrix(0.0, nrow = 0, ncol = 3, dimnames = list(NULL, c("dimension", "birth", "death")))
		for (i in seq(1L, k)){
			{ nr <- nrow(R[[i]]); nc <- ncol(R[[i]]) }
			P <- pivots(R[[i]])
			if (length(P) > 0){
				barcodes <- rbind(barcodes, cbind(s_dims[i]-1L, P))
			}
			
			## If requested, record essential classes which never die 
			if (essential){
				E <- setdiff(seq.int(nc), as.vector(P[,2]))
				if (i < k){
					E <- setdiff(E, pivots(R[[i+1]])[,1])
				}
				if (length(E) > 0){
					barcodes <- rbind(barcodes, cbind(s_dims[i], E, Inf))
				}
			}
		}
		
		## If weights supplied, use report them
		if (!missing(values) && !is.null(values)){
			for (di in s_dims){
				di_intervals <- barcodes[,1] == di
				birth_vals <- values[[match(di, s_dims)]]   
				death_vals <- values[[match(di+1L, s_dims)]]
				barcodes[di_intervals, 2] <- birth_vals[barcodes[di_intervals,2]]
				if (!is.null(death_vals)){
					non_essential <- di_intervals & barcodes[,3] != Inf
					barcodes[non_essential, 3] <- death_vals[barcodes[non_essential,3]]
				} else {
					barcodes[di_intervals, 3] <- Inf
				}
			}
		}
		
		## If requested, remove 0-persistence pairs 
		if (collapse){
			persistence <- abs(barcodes[,2] - barcodes[,3])
			barcodes <- barcodes[persistence > .Machine$double.eps,,drop=FALSE]
		}
	} else {
		stop("Invalid input to 'x'--must be matrix or decomposition.")
	}
	
	## Fix order to the lexicographical 
	barcodes <- barcodes[order(barcodes[,1], barcodes[,2], barcodes[,3]),,drop=FALSE]
	
	## Only grab the dimensions requested (todo: make this more efficient)
	barcodes <- barcodes[barcodes[,1] %in% dims,,drop=FALSE]
	return(barcodes)
}

# permute_decomp <- function(x, i, j){
# 	x$R
# }


# balls_to_bins <- function(){
# 	
# }
# 
# bins_to_balls <- function(p, ){
# 	
# }


#' Reduction <- R6::R6Class("Reduction", list(
#'   filtration = list(),
#'   V = Matrix::rsparsematrix(nrow = 0, ncol = 0, density = 0.0), 
#'   field = "mod2"
#' ))
#' 
#' # Reduction$set("active", "R", function(value){
#' # 	if (!missing(value)){ stop("R is a read-only property. ") }
#' # 	
#' # })
#' 
#' 
#' Reduction$set("public", "pivots", function(R){
#' 	row_pivots <- apply(B1, 2, function(x) {
#' 		idx <- which(as.vector(x) != 0)
#' 		ifelse(length(idx) == 0, 0, tail(idx, 1L))
#' 	})
#' 	
#' })
#' 
#' Reduction$set("public", "reduce", function(D){
#' 	V <- Matrix::Matrix(diag(ncol(D)), sparse = TRUE)
#' 	J <- seq(ncol(D))
#' 	P <- vector(mode = "integer", length = 0L)
#' 	E <- vector(mode = "integer", length = 0L)
#' 	.pivot_row <- function(x) { idx <- which(as.vector(x) != 0); ifelse(length(idx) == 0, 0, tail(idx, 1L)) }
#' 	.search_p <- function(pivot_j, j){
#' 		if (pivot_j == 0 || j == 1 || length(P) == 0){ return(0) }
#' 		potential_ranks <- which(P < rank_szudsick(matrix(c(pivot_j, j))))
#' 		if (length(potential_ranks) == 0){ return(0) }
#' 		potential_pivots <- unrank_szudsick(P[potential_ranks])
#' 		true_pivots <- which(apply(potential_pivots, 2, function(pivot){ pivot[1] == pivot_j }))
#' 		if (length(true_pivots) == 0){ return(0) }
#' 		return(rank_szudsick(matrix(potential_pivots[,head(true_pivots, 1L)])))
#' 	}
#' 	for (j in J){
#' 		R_j <- D[,j,drop=FALSE]
#' 		pivot_j <- .pivot_row(R_j)
#' 		sk <- .search_p(pivot_j, j)
#' 		while(sk != 0){
#' 			pivot_k <- unrank_szudsick(sk)[1]
#' 			k <- unrank_szudsick(sk)[2]
#' 			R_k <- D %*% V[,k,drop=FALSE]
#' 			lambda <- R_j[pivot_j]/R_k[pivot_k]
#' 			R_j <- R_j - lambda * R_k 
#' 			V[,j,drop=FALSE] <- V[,j,drop=FALSE] - lambda * V[,k,drop=FALSE]
#' 			pivot_j <- .pivot_row(R_j)
#' 			sk <- .search_p(pivot_j, j)
#' 		} 
#' 		if (pivot_j != 0){
#' 			P <- sort(c(P, rank_szudsick(matrix(c(pivot_j, j)))))
#' 		} else { E <- c(E, j) }
#' 	}
#' })
#' 
#' # restore_l <- function(R, V, I, state = c(FALSE, FALSE, FALSE, FALSE)){
#' # 	.pivot_row <- function(x) { idx <- which(as.vector(x) != 0); ifelse(length(idx) == 0, 0, tail(idx, 1L)) }
#' # 	{ d_low <- .pivot_row(R[,I[1]]); d_colR <- R[,I[1],drop=FALSE]; d_colV <- V[,I[1],drop=FALSE] }
#' # 	if (length(I) == 1){ return(list(R=R, V=V, d_colR=d_colR, d_colV=d_colV)) }
#' # 	for (k in seq(2, length(I))){
#' # 		{ new_d_low <- .pivot_row(R[,I[k]]); new_d_colR <- R[,I[k],drop=FALSE]; new_d_colV <- V[,I[k],drop=FALSE] }
#' # 		if (!state[3]){
#' # 			R[,I[1]] <- (R[,I[1]] + R[,I[k]]) %% 2
#' # 			V[,I[1]] <- (V[,I[1]] + V[,I[k]]) %% 2
#' # 		} else {
#' # 			R[,I[k]] <- (R[,I[k]] + d_colR) %% 2
#' # 			V[,I[k]] <- (V[,I[k]] + d_colV) %% 2
#' # 		}
#' # 		test <- ifelse(state[2], new_d_low > d_low, d_low > new_d_low)
#' # 		if (test){ # && new_d_low != 0
#' # 			{ d_low <- new_d_low; d_colR <- new_d_colR; d_colV <- new_d_colV }
#' # 		}
#' # 	}
#' # 	return(list(R=R, V=V, d_colR=d_colR, d_colV=d_colV))
#' # }
#' 

## permutes 'M' such that column/row i is moved to the jth position
#' permute_move
#' @description Applies a move permutation to a matrix.
#' @details Applies the permutation given by moving column or row 'i' (or both) to column
#' or row 'j' (or both).
#' @export
permute_move <- function(M, i, j, dims = c("both", "rows", "cols")){
	if (i == j){ return(M) }
	.move_permutation <- function(i, j, n){
		if (i > j){
			p <- c(i, setdiff(seq(j, n), i))
			if (j != 1){ p <- c(seq(j-1), p) }
		} else {
			p <- c(setdiff(seq(j), i), i)
			if (length(p) < n){ p <- c(p, seq(j+1, n)) }
		}
		return(p)
	}
	## Apply the permutation to a matrix if the dimensions aren't null
	if (!is.null(dim(M))){
		stopifnot(!((missing(dims) || dims == "both") && max(c(i,j)) > min(dim(M))))
		stopifnot(!((dims == "rows") && max(c(i,j)) > nrow(M)))
		stopifnot(!((dims == "cols") && max(c(i,j)) > ncol(M)))
		if (missing(dims) || dims == "both"){
			pr <- .move_permutation(i,j,nrow(M))
			pc <- Matrix::invPerm(.move_permutation(i,j,ncol(M)))
			return(pbgrad::permutation_matrix(pr) %*% M %*% pbgrad::permutation_matrix(pc))
		} else if (dims == "rows"){
			pr <- .move_permutation(i,j,nrow(M))
			return(pbgrad::permutation_matrix(pr) %*% M)
		} else if (dims == "cols"){
			pc <- Matrix::invPerm(.move_permutation(i,j,ncol(M)))
			return(M %*% pbgrad::permutation_matrix(pc))
		} else { stop("Invalid dims given.") }
	} else if (is.vector(M)){
		p <- .move_permutation(i,j,length(M))
		return(M[p])
	} else {
		stop("Input must be a matrix or a vector.")
	}
}


# 
# # ## Moves the column at position i to position j
# # permute_col <- function(M, i, j){
# # 	p <- c(setdiff(seq(j), i), i)
# # 	if (j < ncol(M)){ p <- c(p, seq(j+1, ncol(M))) }
# # 	pc <- pbgrad::permutation_matrix(Matrix::invPerm(p))
# # 	return(M %*% pc)
# # }
# #
# # permute_row <- function(M, i, j){
# # 	p <- c(setdiff(seq(j), i), i)
# # 	if (j < nrow(M)){ p <- c(p, seq(j+1, nrow(M))) }
# # 	pr <- pbgrad::permutation_matrix(p)
# # 	return(pr %*% M)
# # }
# #
# # # Moves the i'th row/column to the j'th row/column
# # permute_ij <- function(M, i, j){
# # 	permute_row(permute_col(M, i, j), i, j)
# # }
# #
# # move_left <- function(R, V, i, j){
# # 	stopifnot(i >= j, ncol(R) == ncol(V), j < ncol(R), j <= ncol(R))
# # 	if (i == j){ return(list(R=R,V=V)) }
# # 	.pivot_row <- function(x) { idx <- which(as.vector(x) != 0); ifelse(length(idx) == 0, 0, tail(idx, 1L)) }
# # 	I <- seq(j,i)[which(V[j,seq(j,i)] != 0)]
# # 	# I <- seq(j,i)[which(V[seq(j,i),i] != 0)]
# # 	stopifnot(length(I) > 0)
# #
# # 	## Solution 1 for obtaining J
# # 	J <- which(R[i,] != 0)
# # 	J <- intersect(J, with(list(low_j=apply(R[,J,drop=FALSE], 2, .pivot_row)), { J[low_j >= j & low_j <= i] }))
# #
# # 	## Restore R, V
# # 	# i_d <- which.min(apply(R[,I,drop=FALSE], 2, .pivot_row))
# # 	result1 <- restore(R, V, I)
# # 	R_new <- result1$R
# # 	V_new <- result1$V
# #
# # 	## Restore R_1, V_1
# # 	if (length(J) > 1){
# # 		result2 <- restore(R_new, V_new, J)
# # 		R_new <- result2$R
# # 		V_new <- result2$V
# # 	}
# #
# # 	## Generate the permutation
# # 	Rh <- permute_move(R_new, i, j)
# # 	Vh <- permute_move(V_new, i, j)
# #
# # 	## Replace newly moved columns with the donors and return
# # 	Rh[,i] <- permute_move(result1$d_colR, i, j, dims = "rows")
# # 	Vh[,i] <- permute_move(result1$d_colV, i, j, dims = "rows")
# # 	return(list(R=Rh, V=Vh))
# # }
# 
# # move_left <- function(R, V, i, j, state=c(FALSE, FALSE, FALSE, FALSE)){
# # 	# stopifnot(j < i, ncol(R) == ncol(V), i < ncol(R), j <= ncol(R))
# # 	if (i == j){ return(list(R=R,V=V)) }
# # 	.pivot_row <- function(x) { idx <- which(as.vector(x) != 0); ifelse(length(idx) == 0, 0, tail(idx, 1L)) }
# # 	if (state[4]){
# # 		I <- seq(j,i)[which(V[i,seq(j,i)] != 0)]
# # 	} else {
# # 		I <- seq(j,i)[which(V[seq(j,i), i] != 0)]
# # 	}
# # 	stopifnot(length(I) > 0)
# #
# # 	## Solution 1 for obtaining J (should be 9,10)
# # 	J <- which(R[i,] != 0)
# # 	J <- intersect(J, with(list(low_j=apply(R[,J,drop=FALSE], 2, .pivot_row)), { J[low_j >= i & low_j <= j] }))
# # 	# print(length(J))
# #
# # 	## Restore R, V
# # 	# i_d <- which.min(apply(R[,I,drop=FALSE], 2, .pivot_row))
# # 	if (!state[1]){ I <- rev(I) }
# # 	result1 <- restore_l(R, V, I)
# # 	R_new <- result1$R
# # 	V_new <- result1$V
# #
# # 	## Restore R_1, V_1
# # 	if (length(J) > 1){
# # 		print("J")
# # 		result2 <- restore(R_new, V_new, J, state)
# # 		R_new <- result2$R
# # 		V_new <- result2$V
# # 	}
# #
# # 	## Generate the permutation
# # 	Rh <- permute_move(R_new, i, j)
# # 	Vh <- permute_move(V_new, i, j)
# #
# # 	## Replace newly moved columns with the donors and return
# # 	if (state[5]){
# # 		Rh[,j] <- permute_move(result1$d_colR, i, j, "rows")
# # 		Vh[,j] <- permute_move(result1$d_colV, i, j, "rows")
# # 	} else {
# # 		Rh[,i] <- permute_move(result1$d_colR, i, j, "rows")
# # 		Vh[,i] <- permute_move(result1$d_colV, i, j, "rows")
# # 	}
# # 	return(list(R=Rh, V=Vh))
# # }
# 
# 
# 
# # restore(reduction$R, reduction$V, I = )
# 
# # wrapper to generate a boundary matrix
# Reduction$set("public", "boundary_matrix", function(dim = "all"){
# 	stopifnot(missing(dim) || dim == "all" || is.numeric(dim))
# 	if (missing(dim) || dim == "all"){
# 		S <- sapply(self$filtration, length)
# 		ns <- sum(S)
# 		M <- Matrix::Matrix(0, nrow = ns, ncol = ns, sparse = TRUE)
# 		if (length(S) <= 1){ return(M) } # otherwise length(S) >= 2
# 		cs <- cumsum(S)
# 		for (j in seq(2, length(self$filtration))){
# 			col_idx <- setdiff(seq(cs[j]), seq(cs[j-1]))
# 			row_idx <- setdiff(seq(cs[j-1]), seq(cs[j-2]))
# 			M[row_idx, col_idx] <- self$boundary_matrix(j-1)
# 		}
# 		nv <- S[1]
# 		sn <- unlist(sapply(seq(length(self$filtration)), function(i){ simplex_to_str(unrank_lex(self$filtration[[i]], n = nv, k = i)) }))
# 		dimnames(M) <- list(sn, sn)
# 		return(M)
# 	} else if (is.numeric(dim)){
# 		d <- dim + 1L
# 		if (length(self$filtration) < d){ return(matrix(numeric(length = 0L), nrow = 0L, ncol = 0L))}
# 		domain_simplices <- self$filtration[d]
# 		if (d == 1){
# 			nv <- length(self$filtration[[1]])
# 			M <- Matrix::Matrix(matrix(numeric(length = 0L), nrow = 0L, ncol = nv, dimnames = list(NULL, self$filtration[[1]])))
# 			return(M)
# 		} else {
# 			nv <- length(self$filtration[[1]])
# 			{ C <- self$filtration[[d]]; R <- self$filtration[[d-1]] }
# 			{ nc <- length(C); nr <- length(R) }
# 			face_indices <- lapply(C, function(c){ rank_lex(combn(unrank_lex(x = c, n = nv, k = d), m = d-1), n = nv) })
# 			j <- rep(seq(length(face_indices)), each = d)
# 			i <- unlist(lapply(face_indices, function(lex_idx){ match(lex_idx, R) }))
# 			x <- rep((-1)^seq(0, d-1), nc)
# 			cn <- simplex_to_str(unrank_lex(C, n = nv, k = d))
# 			rn <- simplex_to_str(unrank_lex(R, n = nv, k = d-1))
# 			D <- Matrix::sparseMatrix(i, j, x = x, dims = c(nr, nc), dimnames = list(rn,cn))
# 			return(D)
# 		}
# 	}
# })

