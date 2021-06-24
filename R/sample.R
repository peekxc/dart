#' sample_filtration
#' @description Samples a given filtration, respecting the face relation
#' @export
sample_filtration <- function(S){
	stopifnot(is.character(S))
	ds <- sapply(str_to_simplex(S, simplify = FALSE), length)
	P <- seq(length(ds))
	for (u in unique(ds)){
		P[ds == u] <- P[ds == u][sample(sum(ds == u))]
	}
	S <- S[P]
	stopifnot(!as.logical(anyDuplicated(S)))
	for (i in seq(length(S))){
		sigma <- as.vector(str_to_simplex(S[i]))
		d <- length(sigma)
		face_idx <- match(simplex_to_str(combn(sigma, d-1)), S)
		if (!is.na(face_idx) && any(face_idx > i)){
			j <- max(face_idx)
			tmp <- S[j]
			S[j] <- S[i]
			S[i] <- tmp
		}
	}
	return(S)
}

#' Random Vietoris-Rips complex 
#' @description Generates a Rips complex from sampled points on the unit square.
#' @export
r_rips_complex <- function(n, radius, coords = FALSE, ...){
	xy <- cbind(runif(n), runif(n))
	d <- parallelDist::parallelDist(xy, method = "euclidean")
	R <- simplextree::rips(d, eps = 2*radius, ...)
	if (coords){ attr(R, "coords") <- xy }
	return(R)
}