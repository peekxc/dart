## filtration.R 
## Represents a set of S3 methods and generics for representing filtered simplicial complexes
Filtration <- R6::R6Class(
  "Filtration",
  public = list(
  	complex = NULL, 
    initialize = function(simplices, grades = NULL, persistent=TRUE) {
    	stopifnot(is.list(simplices) || inherits(simplices, "Rcpp_SimplexTree"))
    	## Choose what to do based on the input
    	if (is.list(simplices)){
    		## If list, assume grades given in filtration order
    		self$complex <- simplextree::simplex_tree(simplices)
    		if (missing(grades) || is.null(grades)){ grades <- seq_len(sum(self$complex$n_simplices)) }
    		if (length(grades) == self$complex$n_simplices[1]){
    			grades <- grades[match(simplices[sapply(simplices, length) == 1], self$complex$vertices)]
    		} else if (length(grades) == self$complex$n_simplices[2]){
    			grades <- grades[match(simplices[sapply(simplices, length) == 2], as.list(k_simplices(self$complex, 1L)))]
    		} else {
    			stopifnot(length(grades) == sum(self$complex$n_simplices))
    			grades <- grades[match(simplices,as.list(level_order(self$complex)))]
    		}
    	} else {
    		## Assume grades given in shortlex order
    		self$complex <- simplices
    		if (missing(grades) || is.null(grades)){ grades <- seq_len(sum(self$complex$n_simplices)) }
    		stopifnot(is.numeric(grades), length(grades) %in% c(head(self$complex$n_simplices, 2), sum(self$complex$n_simplices)))
    	}
    	private$.is_persistent <- as.logical(persistent)
      private$.filtration <- new(dart:::ImplicitFiltration, self$complex$as_XPtr(), grades)
    }, 
    finalize = function() {
			# 	    print("finalizing filtration")
			#     	if (private$.is_persistent){
			#     		saveRDS(simplextree::serialize(self$complex), file = private$store_file)
			#     	}
	  }
  ), 
  private = list(
  	.filtration = NULL, 
  	.is_persistent = NULL, 
  	.ptr = NULL, 
  	.stored_file = NULL
  )
)

Filtration$set("public", "print", function(...){
	print("Filtration")
	invisible(self)
})

Filtration$set("public", "at", function(i){
	sx <- private$.filtration$simplex_at(i-1L)
	gr <- private$.filtration$grade_at(i-1L)
	return(list(simplex=sx, grade=gr))
})

Filtration$set("active", "simplices", function(value){
	if (missing(value)){
		return(private$.filtration$simplices)
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
	invisible(self)
})

Filtration$set("active", "grading", function(value){
	if (missing(value)){
		return(private$.filtration$grading)
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})

Filtration$set("public", "k_simplices", function(k){
	return(private$.filtration$k_simplices(k))
})
Filtration$set("public", "k_grades", function(k){
	return(private$.filtration$k_grades(k))
})

Filtration$set("active", "ranks", function(value){
	if (missing(value)){
		as.vector(unlist(private$.filtration$simplex_ranks))+1L
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})

Filtration$set("active", "shortlex_order", function(value){
	if (missing(value)){
		as.vector(private$.filtration$shortlex_perm)+1L
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})

Filtration$set("active", "dimensions", function(value){
	if (missing(value)){
		as.integer(as.vector(unlist(private$.filtration$simplex_dims)))
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})

Filtration$set("active", "vertices", function(value){
	if (missing(value)){
		as.vector(unlist(private$.filtration$k_simplices(0)))
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})

Filtration$set("active", "edges", function(value){
	if (missing(value)){
		return(do.call(cbind, private$.filtration$k_simplices(1)))
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})
Filtration$set("active", "triangles", function(value){
	if (missing(value)){
		return(do.call(cbind, private$.filtration$k_simplices(2)))
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})
Filtration$set("active", "quads", function(value){
	if (missing(value)){
		return(do.call(cbind, private$.filtration$k_simplices(3)))
	} else {
		stop("Filtrations are read-only once created.", call. = FALSE)
	}
})

#' Filtered simplicial complex
#' @description Creates a filtered simplicial complex.
#' @param simplices list of simplices in filtration order, or a simplex tree. 
#' @param grades numeric vector of grades of 'weights' of simplices, in filtration order 
#' if \code{simplices} is a list or shortlex order is a simplex tree. 
#' @details This function creates a (simplicial) filtration, i.e. a simplicial complex filtered 
#' according to some choice of grading function which induces a partial order on the simplices of the complex. 
#' Given \code{simplices} and a (possibly much smaller) vector of \code{grades} or 'weights', this function 
#' efficiently builds and implicitly stores the simplices and their grades in the filtration order. If
#' \code{grades} are not totally ordered, then a total order is imposed using the lexicographical ordering
#' of on the facet poset of \code{simplices}. \cr
#' \cr 
#' The filtration object itself acts as a light encoding on the input filtration; the simplices of the filtration 
#' are stored implicitly using their rank in the combinatorial number system (see \link{rank_combn}), and the only 
#' the minimal grades are stored. That is, the grades for the k-simplices are computed on demand from the set of 
#' minimal grades and the k-simplices for k > 0 are generated on the fly only when needed. Although the memory complexity 
#' of storing all of the ranks of \code{simplices} matches the size \code{simplices} itself, in practice the
#' the combinatorial number system provides a very cache-efficient means of storage. 
#' \cr
#' To create a filtration, minimally the \code{simplices} must be provided, and optionally the \code{grades} if those 
#' need be stored as well. If \code{simplices} is a list of simplices, both \code{simplices} and \code{grades} are assumed 
#' to be in filtration order; alternatively, if \code{simplices} is a simplex tree, \code{grades} are assumed to be 
#' shortlex order (see \link[simplextree]{level_order}). If the length of \code{grades} matches the number of vertices
#' in the filtration, a lower star filtration is built and only the vertex grades are stored; if the length of \code{grades} 
#' matches the number of edges, a flag filtration is built and only the edge grades are stored; for all other filtrations, 
#' the length of \code{grades} must match the number of simplices given, in which case the full vector is stored. \cr
#' \cr
#' In all cases, the k-simplices are first converted to their lexicographical ranks in (fast) \eqn{O(mk)} time, and then 
#' the corresponding filtration with \code{m} simplices is created in \eqn{O(m \log m)} time. The returned filtration provides 
#' efficient access to simplices and grades in the filtration order. Specifically, k-simplices can be accessed at any index 
#' \code{i} in at most O(k) time, depending on type of filtration built (O(1) access for generic filtrations). The grade 
#' associated with a k-simplex is computed in at most \eqn{O(k log m[k])} time, again depending on the type of filtration
#' built. The grading for all k-simplices is computed separately in at most \eqn{O(k m[k])} time. \cr
#' \cr
#' @return a \code{Filtration} \link[R6]{R6} object. 
#' @export
filtration <- function(simplices, grades = NULL, persistent = TRUE){
	return(Filtration$new(simplices, grades, persistent))
}

#' Rips filtration
#' @param x either a row-oriented point cloud or a \link[stats]{dist} object.
#' @param diameter diameter of balls to place at every point. Defaults to the 2 * < enclosing radius of \code{x} >.
#' @param dim maximum dimension to expand the rips complex to do. Defaults to 1, producing the neighborhood graph of \code{x}. 
#' @export
rips_filtration <- function(x, diameter = "default", dim = 1L, ...){
	if (inherits(x, "dist")){
		d <- x
	} else {
		d <- parallelDist::parallelDist(x)
	}
	if (missing(diameter) || diameter == "default"){
		diameter <- 2.0*simplextree::enclosing_radius(d)
	}
	stopifnot(is.numeric(dim), is.numeric(diameter))
	n <- inverse.choose(length(d), k = 2)
	ind_to_insert <- which(d <= diameter)
	edges <- unrank_combn(ind_to_insert, n, 2L)
	st <- simplextree::simplex_tree() %>% 
			simplextree::insert(as.list(seq(n))) %>% 
			simplextree::insert(edges) %>% 
			simplextree::expand(k = dim)
	return(filtration(simplices = st, grades = d[ind_to_insert], ...))
}



# 
# create_thing <- function(){
# 	## assumes you want to make a variable named 'st' 
# 	.st <- simplextree::simplex_tree(1:5)
# 	attr(.st, "ptr") <- .st$as_XPtr()
# 	file_to_store <- tempfile()
# 	xptr::register_xptr_finalizer(attr(.st, "ptr"), function(x){ 
# 		print("Saving to file")
# 		print(x)
# 		saveRDS(simplextree::serialize(.st), file = file_to_store)
# 	}, onexit = TRUE)
# 	
# 	delayedAssign("st", {
# 		if (xptr::is_null_xptr(attr(.st, "ptr"))){
# 			print("reloading simplex tree")
# 			.st <- readRDS(file_to_store)
# 		}
# 		print(attr(.st, "ptr"))
# 		.st
# 	})
# 	# makeActiveBinding("st", function(x){ 
# 	# 	if (missing(x)){
# 	# 		if (xptr::is_null_xptr(attr(.st, "ptr"))){
# 	# 			return(readRDS(file_to_store))
# 	# 		}
# 	# 		print(attr(.st, "ptr"))
# 	# 		return(.st)
# 	# 	} else {
# 	# 		stop("cant reassign st")
# 	# 	}
# 	# }, env = .GlobalEnv)
# 	# on.exit({
# 	# 	to_save <- simplextree::serialize(st)
# 	# 	to_save$file <- tempfile()
# 	# 	saveRDS(to_save, file = to_save$file)
# 	# }, add = TRUE)
# 	# st_ptr <- st$as_XPtr()
# 
# }
# 
# ## See https://github.com/r-lib/R6/issues/42
# ST <- R6::R6Class(
#   classname = "ST",
#   public = list(
#     ptr = NULL,
#     st = NULL, 
#     store_file = NULL, 
#     initialize = function(...) {
#     	print('init')
#       self$st <- simplextree::simplex_tree(...)
# 			self$ptr <- self$st$as_XPtr()
# 			self$store_file <- file.path("~/st.rds")
# 			# xptr::register_xptr_finalizer(attr(.st, "ptr"), function(x){ 
# 			# 	print("Saving to file")
# 			# 	print(x)
# 			# 	saveRDS(simplextree::serialize(.st), file = file_to_store)
# 			# }, onexit = TRUE)
#     }, 
#     finalize = function() {
# 	    print("finalizing simplex tree")
#     	saveRDS(simplextree::serialize(self$st), file = self$store_file)
# 	  }, 
#     check_deleted = function(){
# 	  	if (xptr::is_null_xptr(self$ptr)){
#     		self$st <- simplextree::deserialize(readRDS(self$store_file))
#     		self$store_file <- file.path("~/st.rds")
#     		self$ptr <- self$st$as_XPtr()
# 	  	}
#     	invisible(self)
#     },
#     print = function(...) {
#     	self$check_deleted()
#     	print(self$st)
#     	invisible(self)
#     }
#   ),
#   cloneable = TRUE,
#   active = list(
#     tree = function(value) {
#       if (missing(value)) {
#       	if (xptr::is_null_xptr(self$ptr)){
#       		self$st <- simplextree::deserialize(readRDS(self$store_file))
#       		self$store_file <- file.path("~/st.rds")
#       		self$ptr <- self$st$as_XPtr()
#       	}
#       }
#     	return(self$st)
#     }
#   ) 
# )
# st <- ST$new(1:5)
# 
# # Circle2$set("public", "clone", function(){
# # 	print("This is a custom clone method!") # Test print statement
# # 	Circle2$new(self$radius)
# # })
# c1 <- Circle2$new(10)
# ClassFilter <- function(x){ inherits(get(x), 'Circle2') }
# Objs <- Filter( ClassFilter, ls(all.names = TRUE) )
# # implicit_filtration <- function(){
# # 	new(ImplicitFiltration, )
# # }
# 
# # When implementing a vector class, you should implement these methods: 
# # length, [, [<-, [[, [[<-, c. (If [ is implemented rev, head, and tail should all work).
# 
# test_pair <- function(x){
# 	data.frame(value=x, name=as.character(x), row.names = "")
# }
# TP <- structure(sapply(seq(10), test_pair, simplify = FALSE), class="test_pair")
# TP2 <- structure(sapply(rev(seq(10)), test_pair, simplify = FALSE), class="test_pair")
# 
# format.test_pair <- function(x, ...){
# 	sprintf("Test pair instance of length: %d", length(x))
# }
# 
# print.test_pair <- function(x, ...){ cat(format(x, ...), "\n") }
# 																		
# 
# # length.test_pair <- function(x, ...){ length(x) }
# 
# `[.test_pair` <- function(x, i){
# 	if (i > length(x)){ return(NULL) }
# 	getElement(x, i)
# }
# 
# `[<-.test_pair` <- function(x, i){
# 	if (i > length(x)){ return(NULL) }
# 	getElement(x, i)
# }

