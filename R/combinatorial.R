#' Kendall distance 
#' @description Measures the kendall distance between two vectors
#' @param x a vector of comparable elements. 
#' @param y a vector of elements of the same kind as \code{x}, or NULL. See details. 
#' @param normalize whether to normalize the distance to [0,1]. 
#' @details This functions computes the Kendall distance between arbitrary vectors, which is equivalent 
#' to counting the number of inversions between the pairs of elements between the two vectors.  
#' Its expected that both \code{x} and \code{y} have comparable elements, and that 
#' the type of the elements can be ordered, e.g. \link[base]{sort} and \link[base]{sort} can be applied. 
#' See below for example inputs. \cr
#' \cr
#' If \code{normalize}=\code{TRUE}, the inversion count is normalized to [0,1].\cr
#' \cr
#' If \code{y} is not supplied, the Kendall distance from \code{x} to the identity 
#' permutation of \code{x} is reported. 
#' @examples 
#' x <- sample(seq(50))
#' 
#' ## This reports the number of inversions needed to sort x into ascending order
#' kendall_dist(x)
#' 
#' ## The distance can also computed between permutations
#' kendall_dist(x, sample(x))
#' 
#' ## Any type that follows a natural order can be used 
#' kendall_dist(letters[1:20], sample(letters[1:20]))
#' @export
kendall_dist <- function(x, y = NULL, normalize = FALSE){
	stopifnot(is.vector(x) || !is.null(dim(x)))
	if (!is.null(dim(x)) && is.matrix(x)){
		return(perm_dist_mat(x, kendall = TRUE, normalize = normalize))
	} else {
		if (missing(y) || is.null(y)){
			if (is.integer(x)){
				inv_count <- inversion_count(x)
			} else {
				inv_count <- inversion_count(match(x, sort(unique(x))))
			}
		} else {
			stopifnot(length(x) == length(y), all(x %in% y))
			inv_count <- inversion_count(match(y,x))
		}
		return(ifelse(normalize, inv_count/choose(length(x), 2L), inv_count))
	}
	stop("Invalid input detected. Please input either an integer matrix or a pair of integer vectors.")
}

#' Spearman distance
#' @description Computes the Spearman's footrule distance. 
#' @param x a vector of comparable elements. 
#' @param y a vector of elements of the same kind as \code{x}, or NULL. See details. 
#' @param normalize whether to normalize the distance to [0,1]. 
#' @details This function computes the Spearmans footrule or "displacement" distance between arbitrary vectors.
#' The Spearman distance can be thought of as a type of l_1 distance between vectors, as it effectively just 
#' measures the sum absolute difference of the ranks of the elements between vectors. 
#' Its expected that both \code{x} and \code{y} have comparable elements, and that 
#' the type of the elements can be ordered, e.g. \link[base]{sort} and \link[base]{sort} can be applied. \cr
#' \cr
#' If \code{normalize}=\code{TRUE}, the inversion count is normalized to [0,1].\cr
#' \cr
#' If \code{y} is not supplied, the Spearman distance from \code{x} to the identity 
#' permutation of \code{x} is reported. 
#' @export
spearman_dist <- function(x, y = NULL, normalize = FALSE){
	stopifnot(is.vector(x) || !is.null(dim(x)))
	if (!is.null(dim(x)) && is.matrix(x)){
		return(perm_dist_mat(x, kendall = FALSE, normalize = normalize))
	} else {
		I <- sort(unique(x))
		if (missing(y) || is.null(y)){
			displacement <- sum(abs(match(I, x) - x))
		} else {
			stopifnot(length(x) == length(y), all(x %in% y))
			displacement <- sum(abs(seq_along(x) - match(x, y)))
		}
		return(ifelse(normalize, displacement/(0.50*(length(x)^2)), displacement))
	}
	stop("Invalid input detected. Please input either an integer matrix or a pair of integer vectors.")
}


#' Rank combinations 
#' @description Ranks combinations lexicographically
#' @param x (k x m) matrix of the k-combinations to rank, or a list of combinations. 
#' @param n choice of number system to rank against. Must be at least as large as the largest value in \code{x}.
#' @details This function is a ranking function which maps k-combinations to their 
#' position in the lexicographically ordered combinatorial number system. The combinatorial number system
#' of degree \code{k} is a bijection from all (\code{n} choose \code{k}) combinations to the natural numbers. 
#' Here, the ranks returned start from 1, and the combinations are expected to not contain 0s. 
#' The choice of \code{n} affects the ranking, but the ranking will always rank combinations in 
#' lexicographical order, e.g. \code{rank_combn(c(1,2,3)) < rank_combn(c(1,2,4)) < ...}. 
#' #' For example, if \code{x}  equals \code{combn(5,2)}, then \code{all(rank_combn(x) == seq_len(choose(5,2)))} 
#' is \code{TRUE}. \cr
#' \cr
#' In general, ranking combinations is faster than unranking them, since the computation reduces
#' essentially to summing up binomial coefficient. In the case where k = 2, the ranking is 
#' particularly efficient, since a closed-form solution is known. 
#' @export
rank_combn <- function(x, n){
	if (is.list(x)){
		stopifnot(all(sapply(x, function(s){ is.numeric(s) & all(s > 0 & s <= n) })))
		return(sapply(x, function(s){ rank_combn_single_rcpp(s-1L, n)+1L }))
	} else {
		if (is.null(dim(x))){ x <- matrix(as.vector(x)) }
		stopifnot(is.numeric(x), is.matrix(x), all(x != 0 & x <= n))
		stopifnot(choose(n, nrow(x)) <= (2^64 - 1L))
		return(rank_combn_rcpp(x-1L, n)+1L)
	}
}

#' Unrank combinations 
#' @description Ranks combinations lexicographically
#' @param x vector of ranks. 
#' @param n choice of number system to rank against. Must be at least as large as the largest value in \code{x}.
#' @param k 
#' @details This function provides the inverse of the ranking function which maps k-combinations to their lexicographical 
#' position in the (\code{n} choose \code{k}) combinatorial number system. For example, if \code{x}
#' equals \code{combn(5,2)}, then \code{all(unrank_combn(seq_len(choose(5,2))) == x)} is \code{TRUE}.
#' The choice of \code{n} affects the ranking, but the ranking will always rank combinations in 
#' lexicographical order. \cr
#' \cr
#' In general, unranking can be slower than ranking. The current implementation uses an implicit O(k*log(n))
#' via a binary search to perform the unranking, except in the case where k = 2, which is particularly efficient 
#' as a closed-form solution is known. 
#' @export
unrank_combn <- function(x, n, k){
	stopifnot(is.numeric(x), length(n) == 1, all(x != 0 & x <= choose(n, max(k))))
	if (length(k) == 1L){
		return(unrank_combn_rcpp(x-1L, n, k)+1L)
	} else {
		stopifnot(length(k) == length(x))
		return(lapply(seq_len(length(x)), function(i){ unrank_combn_single_rcpp(x[i]-1L, n, k[i])+1L }))
	}
}


## ---- LCS/LIS related ----

#' perm_lcs 
#' @description Computes the Longest Common Subsequence (LCS) of a pair of permutations. 
#' @details If both \code{p} and \code{q} are supplied, the LCS of between them is computed 
#' and returned. Alternatively if only \code{p} is supplied, the LCS between p and the identity 
#' permutation is returned. In both cases, computing the LCS reduces to computing the LIS. 
#' @export
perm_lcs <- function(p, q = NULL){
	if (missing(q) || is.null(q)){
		stopifnot(all(p == as.integer(p)))
		return(lis(p))
	} else {
		stopifnot(all(p == as.integer(p)), all(q == as.integer(q)))
		return(p[longest_inc_subseq(match(p,q))])
	}
}

#' lis
#' @description Computes the Longest Increasing Subsequence (LIS) of an integer sequence.
#' @export
lis <- function(p){
	stopifnot(is.numeric(p), as.integer(p) == p)
	stopifnot(all(p[order(p)] == seq_along(p))) ## check is permutation
	return(p[longest_inc_subseq(X = p)])
}

#' all_lcs
#' @description Computes all LCS's between two integer vectors.
#' @param p integer vector.
#' @param q integer vector.
#' @return (k x m) matrix of LCS's, where k = size of LCS and m = number of LCS.
#' @export
all_lcs <- function(p, q){
	stopifnot(is.numeric(p), is.numeric(q), all(as.integer(p) == p), all(as.integer(q) == q))
	return(do.call(cbind, all_lcs(p, q)))
}

#' enumerate_lis
#' @description Enumerates all LIS's of a sequence
#' @param x permutation, as an integer vector. 
#' @export
enumerate_lis <- function(x, max_n = Inf){
	stopifnot(is.numeric(x), all(as.integer(x) == x), all(sort(x) == seq_along(x)))
	if (length(x) <= 1L){ return(x) }
	prev <- function(y, i){ ifelse(any(y < i), max(y[y < i]), NA) }
	succ <- function(y, i){ ifelse(any(y > i), min(y[y > i]), NA) }
	E <- c() # T data structure
	L <- rep(0L, length(x)) # L[i] := length of LIS ending in i
	
	## Stage 1
	for (i in seq_along(x)){
		m <- x[i]
		E <- c(E, m)
		if (!is.na(prev(E, m))){ 
			L[m] <- L[prev(E, m)] + 1L
		} else {
			L[m] <- 1L
		}
		if (!is.na(succ(E, m)) && L[succ(E, m)] == L[m]){
			E <- E[E != succ(E, m)]
		}
	}
	
	## Stage 2
	L <- L[x] ## not mentioned in the paper, but absolutely necessary
	index <- which.max(L)
	S <- vector("integer", length = L[index])
	S[length(S)] <- index
	j <- L[index] - 1L
	for (i in seq(index, 1L)){
		if (L[i] == j){
			S[j] <- i
			j <- j - 1L
		}
	}
	# current_lis <- x[S]
	
	# index <- which.max(L[x])
	# k <- L[x[index]]
	# S <- vector("integer", length = k)
	# S[k] <- index
	# j <- k - 1L
	# for (i in seq(index, 1L)){
	# 	if (L[x[i]] == j){
	# 		S[j] <- i
	# 		j <- j - 1L
	# 	}
	# }
	
	## Enumerate LIS
	# L <- L[Matrix::invPerm(x)]
	L1 <- sapply(seq_along(L), function(j){ 
    idx <- which(L == L[j])
    idx <- idx[idx < j]
    ifelse(length(idx) == 0L, NA, max(idx))
	})
	L2 <- sapply(seq_along(L), function(j){ 
    idx <- which(L == (L[j]-1L))
    idx <- idx[idx < j]
    ifelse(length(idx) == 0L, NA, max(idx))
	})
	
	## Recursively enumerate LIS's
	k <- max(L)
	results <- list()
	cc <- 1L
	enumerate <- function(z, out){
		if (cc > max_n){ return(); }
		if (L[z] == k || ((L[z]+1L) <= k) && z < out[L[z]+1L]){
			out[L[z]] <- z
		} else{ return(NULL); }
		z1 <- L2[z]
		if (is.na(z1)){
			results[[cc]] <<- out
			cc <<- cc + 1L
		} else {
			enumerate(z1, out)
		}
		while(!is.na(L1[z1])){
			enumerate(L1[z1], out)
			z1 <- L2[z1]
		}
	}
	enumerate(tail(x[S],1L), x[S])
	return(do.call(cbind, results))
}




