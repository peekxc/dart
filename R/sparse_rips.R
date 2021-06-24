#' Sparse rips setup procedure
#' @param x matrix of points or 'dist' object.
#' @param alpha the maximum scale used to define the neighborhood graph.
#' @description This function compute the information necessary
#' to efficiently construct sparse rips complexes/filtrations up to scale \code{alpha}
#' at varying levels of sparsification. Specifically, this function constructs the neighborhood 
#' graph representing the 1-skeleton of the rips complex up to scale \code{alpha} of \code{x}, 
#' the greedy permutation of \code{x}, and the covering radii parameters all of which are 
#' used to define sparse rips complexes and filtrations.
#' @export 
sparse_rips <- function(x, alpha = "default"){
	stopifnot(is.matrix(x) || inherits(x, "dist"))
	
	## Start with a greedy permutation
	if ("dist" %in% class(x)){
		n <- inverse.choose(length(x), 2L)
		p <- landmark::landmarks_maxmin(x, num = n)
		d <- x
	} else {
		n <- nrow(x)
		p <- landmark::landmarks_maxmin(x, num = n)
		d <- parallelDist::parallelDist(x)
	}
	
	## Get the insertion radii 
	lambda <- c(max(d), sapply(seq_along(p)[-1], function(i){
		indices <- dart::rank_combn(rbind(p[i], p[1:(i-1)]), n = n)
		min(d[indices])
	}))
	
	## Start off with the full rips complex to filter
	if (missing(alpha) || alpha == "default"){
		alpha <- simplextree::enclosing_radius(d)
	} else {
		stopifnot(is.numeric(alpha))
	}
	G <- simplextree::rips(d, alpha, dim = 1L)
	out <- list(graph=G, permutation=p, radii=lambda, alpha=alpha)
	return(out)
}	

#' Rips sparsification procedure
#' @description Converts a given rips complex into a sparsified rips complex or filtration 
#' @details This function is meant to be used on the outputs of \code{sparse_rips} function. 
#' Given a neighborhood graph \code{x} (as a simplicial complex) with \emph{n} vertices, a greedy permutation 
#' \code{p} and weight vector \code{w} each of length \emph{n}, this function produces a sparse rips 
#' complex (or filtration) up to scale \code{alpha} whose persistence diagram is a (1+\code{epsilon}) approximation 
#' the persistence diagram of \code{x}.  
#' @export 
sparsify_rips <- function(x, p, w, epsilon, alpha, dim = 1L, prune_vertices = FALSE, filtered = FALSE){
	#stopifnot(is.numeric(epsilon), epsilon <= 1.0, epsilon >= 0)
	if (inherits(x, "Filtration")){ 
		x <- x$complex
	}
	stopifnot("Rcpp_SimplexTree" %in% class(x))
	stopifnot(is.vector(p), is.vector(w), is.numeric(w), length(p) == length(w))
	stopifnot(x$n_simplices[1] == length(p))
	sparse_complex <- simplextree::simplex_tree(as.list(x$vertices))
	lambda_inorder <- w[match(seq_along(p), p)]
	weights <- dart:::sparse_complex(
		x$as_XPtr(), sparse_complex$as_XPtr(), 
		lambda_inorder, alpha, epsilon
	)
	simplextree::expand(sparse_complex, k=dim)
	# f$complex$degree(f$complex$vertices)
	if (filtered){
		return(filtration(sparse_complex, weights))
	}
	return(sparse_complex)
}

#' Sparse rips complex
#' @param x point cloud matrix, or a 'dist' object.
#' @param epsilon degree of approximation. See details. 
#' @param alpha diameter to build the filtration to. Default to the enclosing radius of the space. 
#' @param dim highest dimension of simplices to include in the filtration.  
#' @desciption This is a thin wrapper for \code{sparse_rips} and \code{sparsify_rips} which 
#' are used in tandem to construct a sparse rips complex. The resulting rips complex at scale 
#' \code{alpha} is a sparse subcomplex of the rips complex of \code{x}, which itself is a subcomplex 
#' of the sparse rips complex at scale (1+\code{epsilon} * \code{alpha}).
#' @references http://donsheehy.net/research/cavanna15geometric.pdf
#' @return a simplextree object.
#' @export
sparse_rips_complex <- function(x, epsilon=0.10, alpha = "default", dim = 1L, ...){
	g <- sparse_rips(x, alpha=alpha)
	f <- sparsify_rips(g$graph, g$permutation, g$radii, epsilon, g$alpha, dim, filtered=FALSE, ...)
	return(f)
}

#' Sparse rips filtration
#' @param x point cloud matrix, or a 'dist' object.
#' @param epsilon degree of approximation. See details. 
#' @param alpha diameter to build the filtration to. Default to the enclosing raidus of the space. 
#' @param dim highest dimension of simplices to include in the filtration. 
#' @description This is a thin wrapper for \code{sparse_rips} and \code{sparsify_rips} which 
#' are used in tandem to construct a sparse rips filtration. The resulting rips filtration at scale 
#' \code{alpha} has a persistence diagram which is a (1+\code{epsilon})-approximation of the persistence
#' diagram of the rips filtration of \code{x} at scale \code{alpha}.
#' @references http://donsheehy.net/research/cavanna15geometric.pdf
#' @return a filtration object.
#' @export
sparse_rips_filtration <- function(x, epsilon=0.10, alpha = "default", dim = 1L, ...){
	g <- sparse_rips(x, alpha=alpha)
	f <- sparsify_rips(g$graph, g$permutation, g$radii, epsilon, g$alpha, dim, filtered=TRUE, ...)
	return(f)
}


#' Geodesic distance 
#' @description Constructs geodesic distances on a neighborhood graph
#' @param x point cloud matrix 
#' @param k number of k-nearest neighbors to form the neighborhood graph with/ See details. 
#' @param r radius of balls to center at each point to from the neighborhood graph with. See details. 
#' @param warn_disconnected whether to warn the user if the neighborhood graph is not connected. 
#' @details TODO 
#' @export
geodesic_dist <- function(x, k = 15L, r = NULL, warn_disconnected=TRUE){
	stopifnot(is.matrix(x))
	W <- matrix(0L, nrow = nrow(x), ncol = nrow(x))
	if (missing(k) && !missing(r)){
		knn <- RANN::nn2(x, x, radius = r)
		eps_inf <- sqrt(.Machine$double.xmax)
		for (i in seq(nrow(x))){
			non_inf <- knn$nn.dists[i,-1] != eps_inf
			W[i,knn$nn.idx[i,-1][non_inf]] <- knn$nn.dists[i,-1][non_inf]
		}
	} else {
		stopifnot(is.numeric(k))
		knn <- RANN::nn2(x, x, k = 15L)
		for (i in seq(nrow(x))){
			W[i,knn$nn.idx[i,-1]] <- knn$nn.dists[i,-1]
		}
	}
	W <- pmax(W + t(W))
	G <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
	d <- as.dist(igraph::distances(G))
	if (any(d == Inf)){
		warning("More than one connected compoinent detected")
	}
	return(d)
	# min_eps <- min(igraph::E(igraph::mst(G))$weight)
}



