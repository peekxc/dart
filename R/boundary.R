#' @title Boundary matrix
#' @description Constructs filtration boundary matrices
#' @param x list of simplices, simplex tree, or a filtration.
#' @param dim dimension of boundary matrix desired
#' @param labeled whether to add simplex labels to the corresponding matrices.
#' @param field intended 
#' @details This function computes and returns the elementary chains associated with a filtered simplicial complex.
#' The filtration is expected to be given as list of a simplices.  
#' @import Matrix
#' @export
boundary_matrix <- function(x, dims="all", labeled=c("none", "columns", "both"), field=c("mod2", "integer", "real")){
	if (is.list(x)){ x <- filtration(x) }
	show_labels  <- !missing(labeled) || (labeled == TRUE)
	labeled_true <- show_labels && (labeled == TRUE)
	if (inherits(x, "Rcpp_SimplexTree")){
		if (missing(dims) || dims == "all"){
			D <- dart:::boundary_matrix_st_full(x$as_XPtr())
			if (show_labels && (labeled_true || labeled %in% c("columns", "both"))){
				colnames(D) <- simplex_to_str(as.list(simplextree::level_order(x)))	
				if (labeled == "both" || labeled_true){
					rownames(D) <- colnames(D)
				}
			}
		} else {
			D <- lapply(dims, function(k){ 
				m <- dart:::boundary_matrix_st(x$as_XPtr(), k) 
				if (show_labels && (labeled_true || labeled %in% c("columns", "both"))){
					colnames(m) <- simplex_to_str(as.list(simplextree::k_simplices(x, k)))
					if (labeled == "both" || labeled_true){
						rownames(D) <- colnames(D)
					}
				}
				return(m)
			})
		}
	} else if (inherits(x, "Filtration")){
		if (missing(dims) || dims == "all"){
			p <- Matrix::invPerm(x$shortlex_order)
			D <- dart:::boundary_matrix_st_full(x$complex$as_XPtr())
			D <- D[p,p,drop=FALSE]
			if (show_labels && (labeled_true || labeled %in% c("columns", "both"))){
				colnames(D) <- simplex_to_str(x$simplices)
				if (labeled == "both" || labeled_true){
					rownames(D) <- colnames(D)
				}
			}
		} else {
			D <- lapply(dims, function(k){ 
				m <- dart:::boundary_matrix_st(x$complex$as_XPtr(), k) 
				p <- Matrix::invPerm(order(x$ranks[x$dimensions == k]))
				q <- Matrix::invPerm(order(x$ranks[x$dimensions == k-1L]))
				m <- m[q,p,drop=FALSE]
				if (show_labels && labeled_true || labeled %in% c("columns", "both")){
					colnames(m) <- simplex_to_str(x$k_simplices(k))
					if (labeled == "both" || labeled_true){
						rownames(D) <- colnames(D)
					}
				}
				return(m)
			})
		}
	}
	out <- structure(list(matrix=D), hom_dim=dims, class="fbm")
	return(out)
}

boundary_matrix_backup <- function(){
		# fi_ptr <- x$.__enclos_env__[["private"]]$.filtration$as_XPtr()
		# if (missing(dim) || dim == "all"){
		# 	D <- dart:::boundary_matrix_fi_full(fi_ptr)
		# 	if (labeled){
		# 		colnames(D) <- simplex_to_str(x$simplices)	
		# 		rownames(D) <- colnames(D)
		# 	}
		# } else {
		# 	D <- lapply(dim, function(k){ 
		# 		m <- dart:::boundary_matrix_fi(fi_ptr, k) 
		# 		if (labeled){
		# 			colnames(m) <- simplex_to_str(x$k_simplices(k))	
		# 		}
		# 		return(m)
		# 	})
		# }
		if (!missing(dim) && is.numeric(dim)){
		di <- sapply(x, length)
		nv <- sum(di == 1L)
		D <- vector(mode = "list", length = length(dim))
		for (d in dim){
			f <- do.call(cbind, x[di == d])
			s <- do.call(cbind, x[di == d+1L])
			if (is.null(f) && is.null(s)){
				return(Matrix::Matrix(matrix(0, ncol = 0, nrow = 0), sparse = TRUE))
			}
			f_idx <- rankr::rank_comb(x = f, n = nv)
			m_idx <- do.call(rbind, lapply(seq(ncol(s)), function(i){
				face_idx <- match(rankr::rank_comb(combn(x = s[,i], d), n = nv), f_idx)
				cbind(sort(face_idx), i)
			}))
			ns <- ifelse(is.null(ncol(s)), 0, ncol(s))
			M <- Matrix::sparseMatrix(
				i = m_idx[,1], j = m_idx[,2], 
				x = rep((-1)^(seq(d+1)-1L), ns), 
				dimnames = list(simplex_to_str(f), simplex_to_str(s)), 
				dim = c(ncol(f), ns)
			)
			D[[match(d, dim)]] <- M
		}
		if (length(dim) == 1L){ D <- D[[1L]] }
	} else if (missing(dim) || dim == "all"){
		di <- sapply(x, length)
		m_idx <- do.call(rbind, lapply(seq_along(x), function(i){
			d <- di[i]
			if (d <= 1L){ return(matrix(0, nrow = 0, ncol = 2)) }
			face_idx <- unlist(lapply(combn(x[[i]], d-1L, simplify = FALSE), function(sigma){ match(list(sigma), x) }))
			cbind(face_idx, i)
		}))
		entries <- unlist(lapply(di, function(d){ 
			if (d == 1L) return(NULL) 
			(-1)^(seq(d)-1L) 
		}))
		D <- Matrix::sparseMatrix(
			i = m_idx[,1], j = m_idx[,2], 
			x = entries, 
			dimnames = list(simplex_to_str(x), simplex_to_str(x)), 
			dims = rep(length(x), 2)
		)
		dim <- sort(unique(di))
	}
}

# boundary_matrix.simplex_tree <- function(){
# 	stopifnot(class(x) %in% c("Rcpp_SimplexTree", "Rcpp_Filtration"))
# 	
# 		k1 <- rankr::rank_comb(straverse(k_simplices(R, k = 1), identity))
# 		k2 <- straverse(k_simplices(R, k = 2), function(x){
# 			findInterval(rankr::rank_comb(combn(x, length(x)-1L)), k1)
# 		})
# 		M <- Matrix::sparseMatrix(i = , j = , x = rep((-1)^(seq(d+1)-1L), ncol(s)))
# 	}
# }



#' Prints filtration boundary matrix objects.
#' @param x a filtration boundary matrix object.
#' @param ... unused. 
#' @method print fbm
#' @export print.fbm 
#' @export
print.fbm <- function(x, ...){
	stopifnot("fbm" %in% class(x))
	if (is.list(x$matrix)){
		num_bm <- length(x$matrix)
		bm_str <- do.call(sprintf, append(list(paste(rep("D%d", num_bm), collapse = ", ")), attr(x, "hom_dim")))
		bm_str <- paste0(bm_str, sprintf("-filtration boundary %s", ifelse(num_bm == 1, "matrix", "matrices")))
		sz_str <- paste(sapply(x$matrix, function(bm) do.call(sprintf, append(list("%d x %d"), dim(bm)))), collapse = ", ")
		sz_str <- paste(ifelse(num_bm == 1, "Size:", "Sizes:"), sz_str)
	} else {
		bm_str <- sprintf("Filtration boundary matrix (dim=%s)", paste0(as.character(attr(x, "hom_dim")), collapse = ","))
		sz_str <- sprintf("Size: %d, %d", nrow(x$matrix), ncol(x$matrix))
	}
	writeLines(c(bm_str, sz_str))
}
			# bm_dim <- sapply(seq(num_bm), function(i){ length(str_to_simplex(colnames(D[[i]][,1,drop=FALSE]))) })-1L
	# tri <- rips$simplices[di == 3]
	# m_idx <- do.call(rbind, lapply(seq_along(tri), function(i){
	# 	face_idx <- match(rank_comb(combn(tri[[i]], 2), rips$n_simplices[1]), e_idx)
	# 	cbind(sort(face_idx), i)
	# }))
	# M <- Matrix::sparseMatrix(i = m_idx[,1], j = m_idx[,2], x = rep(c(1.0,-1.0,1.0), length(tri)))
	# rownames(M) <- simplex_to_str(do.call(cbind, rips$simplices[di == 2]))
	# colnames(M) <- simplex_to_str(do.call(cbind, rips$simplices[di == 3]))
	# 
	# x_names <- names(x)
	# D <- matrix(0, nrow = length(x), ncol = length(x))
	# positions <- lapply(seq(length(x)), function(i){
	# 	c_simplex <- x[[i]]
	# 	d <- length(c_simplex)
	# 	if (d == 1){ return(NA) }
	# 	combn(c_simplex, d-1, function(tau){
	# 		Position(function(sigma){ 
	# 			if (length(sigma) != length(tau)){ return(FALSE) }
	# 			return(all(sigma == tau))
	# 		}, x)
	# 	})
	# })
	# for (j in seq(length(positions))){
	# 	idx <- positions[[j]]
	# 	if (length(idx) > 1){
	# 		D[positions[[j]],j] <- 1
	# 	}
	# }
	# return(D)
