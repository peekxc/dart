#' bipersistence
#' @description Computes the information needed for 2D persistence
#' @param x (m x d) matrix of m points representing a point set.
#' @param f function on \code{x} to create bifiltration with. 
#' @param 
#' @export 
bipersistence <- function(x, f, H = 0L, xbin = 10, ybin = 10, max_dist = 1.0, rivet_path = "~/rivet", ...){
	stopifnot(is.matrix(x), length(f) == nrow(x))
	
	## Compute the bigraded Betti numbers
	ro <- bigraded_betti(x, f = dens_v, H = H, xbin = xbin, ybin = ybin, max_dist = max_dist, 
										   rivet_path = rivet_path, ...)
	
	## Use bigraded Betti numbers to compute critical lines -> line arrangment
	LA <- dual_line_arrangement(ro)
	
	# centroids <- t(sapply(LA, sf::st_centroid))

}

#' Line arrangement under point-line duality
#' @description Computes the cell complex of the line arrangement induced by the multigraded Betti numbers. 
#' @param x output of \code{bigraded_betti}
#' @return Polygons making up the 2-cells of the arrangement. 
#' @export
dual_line_arrangement <- function(x){
	## Calculate anchors 
	S <- rbind(x$betti0[,1:2], x$betti1[,1:2])
	A <- compute_anchors(S)
	if (any(A[,1] < 0)){ A[,1] <- -A[,1] }
	
	## Use anchors to get lines making line arrangement. 
	## Determine bounds of arrangement using intersections
	L <- cbind(m=A[,1], b= -A[,2])
	I <- combn(nrow(L),2)
	a_pts <- apply(I, 2, function(ii){
    i <- ii[1]; j <- ii[2]
    a <- L[i,1]; b <- L[i,2]
    c <- L[j,1]; d <- L[j,2]
    xc <- (b-d)/(c-a)
  	if (xc > 0){ return(c(xc, a*xc + b)) }
  	return(NULL)
	})
	if (!is.null(a_pts)){
		L_pts <- do.call(rbind, a_pts)
		L_pts <- L_pts[!(abs(L_pts[,1]) == Inf | abs(L_pts[,2]) == Inf),,drop=FALSE]
	} else {
		L_pts <- matrix(0, nrow = 0, ncol = 2)
	}

	## Compute line arrangement
	if (nrow(L_pts) == 0){
		 LA <- polygonize_la(L, xlim = c(0, max(ro$x_grades)), ylim = max(ro$y_grades)*c(-1,1))
	} else {
		 LA <- polygonize_la(L, xlim = c(0, max(L_pts[,1])), ylim = c(min(c(-A[,2], A[,2])), max(L_pts[,2])))
	}
	return(LA)
}

#' Dual graph
#' @description Computes the dual graph of a line arrangement. 
#' @param arrangement geometry collection of polygons (output from \code{dual_line_arrangement}).
#' @details Uses the relation t(D) %*% D to compute the dual graph, where \code{D} is the boundary matrix 
#' of the polygons. The dual graph has one vertex for each face, and one edge for each edge connecting adjacent
#' faces in the collection. This function is particularly useful for planar subdivisions / line arrangements.  
#' @import kdtools Matrix
#' @return edge list of the dual graph. 
#' @export
dual_graph <- function(arrangement){
	stopifnot("GEOMETRYCOLLECTION" %in% class(arrangement))
	stopifnot(all(sapply(arrangement, function(x) "POLYGON" %in% class(x))))
	
	## Faster dual graph
	all_segments <- lapply(LA, function(p){
		pts <- matrix(unlist(p), ncol=2)
		do.call(rbind, lapply(seq(nrow(pts)-1), function(i){ 
			if (pts[i,1] < pts[i+1,1]){
				c(pts[i,], pts[i+1,]) 
			} else {
				c(pts[i+1,], pts[i,]) 
			}
		}))
	})
	segments_ordered <- kdtools::kd_sort(unique(do.call(rbind, all_segments)))
	face_idx <- lapply(all_segments, function(S){
		apply(S, 1, function(s){ kdtools::kd_lower_bound(segments_ordered, s) })
	})
	jj <- rep(seq_along(LA), times = sapply(face_idx, length))
	B <- Matrix::sparseMatrix(i = unlist(face_idx), j = jj, x = rep(1, length(unlist(face_idx))))
	nz <- which((t(B) %*% B) != 0, arr.ind = TRUE)
	E <- nz[nz[,1] != nz[,2],]
	DG <- unique(cbind(pmin(E[,1], E[,2]), pmax(E[,1], E[,2])))
	
	# ## Get polygon boundary intersections
	# poly_intersects <- combn(length(arrangement), 2, function(x){ 
	# 	length(sf::st_intersects(arrangement[[x[1]]], arrangement[[x[2]]])[[1]]) > 0 
	# })
	# ## Determine if intersections are point-wise or along an edge
	# intersect_at_end <- apply(combn(length(arrangement), 2)[,poly_intersects], 2, function(x){
	# 	int_type <- class(sf::st_intersection(arrangement[[x[1]]], arrangement[[x[2]]]))
	# 	return("LINESTRING" %in% int_type)
	# })
	# ## The edges of the dual graph are the pairs of faces that share an edge 
	# edges <- combn(length(arrangement), 2)[,which(poly_intersects)[intersect_at_end]]
	return(DG)
}


#' Bigraded Betti numbers 
#' @description Returns the bigraded Betti numbers of a bifiltration. 
#' @details The function uses the RIVET software to compute the bigraded betti numbers of a bifiltration 
#' built from point cloud \code{x} and function \code{f}.
#' @param x point set.
#' @param f function values on \code{x}. 
#' @param H homology dimension to compute.
#' @param xbin number of bins to coarsen \code{f} with. 
#' @param ybin number of bins to coarsen the rips filtration values on \code{x} with. 
#' @param max_dist maximum diameter to consider for rips filtration. 
#' @param rivet_path directory where \emph{rivet_console} is located. Defaults to "~/rivet". Must be valid. 
#' @param save.path if supplied, an input file giving the parameters to run RIVET is saved here. 
#' @param xreverse unused. 
#' @param yreverse unused.
#' @export
bigraded_betti <- function(x, f, H=1L, xbin=10L, ybin=10L, max_dist=NULL, rivet_path="~/rivet",
													 save.path=NULL, xreverse=FALSE, yreverse=FALSE){
	stopifnot(is.numeric(f), is.matrix(x), length(f) == nrow(x))
	rivet_input <- list()
	rivet_input <- append(rivet_input, "--datatype points_fn")
	rivet_input <- append(rivet_input, paste0(as.character(f), collapse = ","))
	pts <- apply(x, 1, function(x) paste0(as.character(x), collapse = ","))
	rivet_input <- append(rivet_input, pts)
	if (!missing(save.path) && is.character(save.path)){
		rivet_input_fn <- file.path(save.path)
	} else {
		rivet_input_fn <- file.path(tempfile())
		on.exit(unlink(rivet_input_fn), add = TRUE)
	}
	writeLines(unlist(rivet_input), con = rivet_input_fn)
	
	rivet_opts <- sprintf("-H %d --xbins %d --ybins %d", H, xbin, ybin)
	if (!missing(max_dist) && is.numeric(max_dist)){
		rivet_opts <- paste(rivet_opts, sprintf("--maxdist %g", max_dist))
	}
	if (xreverse){ rivet_opts <- paste(rivet_opts, "--xreverse") }
	if (yreverse){ rivet_opts <- paste(rivet_opts, "--yreverse") }
	rivet_opts <- paste(rivet_opts, "--betti")
	rivet_dir <- normalizePath(rivet_path)
	rivet_output_fn <- file.path(tempfile())
	on.exit(unlink(rivet_output_fn), add = TRUE)
	# if (file.exists(rivet_output_fn)){ file.remove(rivet_output_fn)	}
	rivet_cmd <- paste(file.path(rivet_dir, "rivet_console"), rivet_input_fn, rivet_opts, ">", rivet_output_fn)
	message(rivet_cmd)
	rivet_out <- system(rivet_cmd)

	if (rivet_out == 0){
		out <- read_betti(rivet_output_fn)
		return(out)
	} else {
		stop("RIVET failed to execute properly. Check inputs.")
	}
}

#' Bigraded_betti2
#' @param x list of simplices in filtration order 
#' @param f1 first filter function 
#' @param f2 second filter function
#' @export
bigraded_betti2 <- function(x, f1, f2, H=1L, xbin=10L, ybin=10L, rivet_path="~/rivet",
													 save.path=NULL, xreverse=FALSE, yreverse=FALSE){
	
	# L <- sparse_rips(x)
	# R <- sparsify_rips(x=L$graph, p=L$permutation, w=L$radii, epsilon=0.1, alpha=L$alpha, dim=2, filtered=TRUE)
	# 
	stopifnot(is.numeric(f1), is.numeric(f2), is.list(x), length(f1) == length(x), length(f2) == length(x))
	rivet_input <- list()
	rivet_input <- append(rivet_input, "--datatype bifiltration")
	
	bifiltration <- lapply(seq_along(s), function(i){
		sprintf("%s ; %g %g", paste0(as.character(s[[i]]), collapse = " "), f1[i], f2[i])
	})
	rivet_input <- c(rivet_input, bifiltration)
	# rivet_input <- append(rivet_input, paste0(as.character(f), collapse = ","))
	# pts <- apply(x, 1, function(x) paste0(as.character(x), collapse = ","))
	# rivet_input <- append(rivet_input, pts)
	
	if (!missing(save.path) && is.character(save.path)){
		rivet_input_fn <- file.path(save.path)
	} else {
		rivet_input_fn <- file.path(tempfile())
		on.exit(unlink(rivet_input_fn), add = TRUE)
	}
	writeLines(unlist(rivet_input), con = rivet_input_fn)
	
	rivet_opts <- sprintf("-H %d --xbins %d --ybins %d", H, xbin, ybin)
	if (xreverse){ rivet_opts <- paste(rivet_opts, "--xreverse") }
	if (yreverse){ rivet_opts <- paste(rivet_opts, "--yreverse") }
	rivet_opts <- paste(rivet_opts, "--betti")
	rivet_dir <- normalizePath(rivet_path)
	rivet_output_fn <- file.path(tempfile())
	on.exit(unlink(rivet_output_fn), add = TRUE)
	# if (file.exists(rivet_output_fn)){ file.remove(rivet_output_fn)	}
	rivet_cmd <- paste(file.path(rivet_dir, "rivet_console"), rivet_input_fn, rivet_opts, ">", rivet_output_fn)
	message(rivet_cmd)
	rivet_out <- system(rivet_cmd)

	if (rivet_out == 0){
		out <- read_betti(rivet_output_fn)
		return(out)
	} else {
		stop("RIVET failed to execute properly. Check inputs.")
	}
}




# Reads in the Betti number file produced by RIVET 
# x := file location 
# value := whether to return grades associated with the betti number of the original (raw) indices 
# return 1) (0,1,2) Bigraded Betti numbers, the non-zero Hilbert (dimension) function values,
# and the grading for the module
read_betti <- function(x, value=TRUE){
	
	## Read the filename
	rivet_out <- readLines(x)
	
	## Get x-grades 
	xb <- which(rivet_out == "x-grades")+1L
	xe <- which(rivet_out == "y-grades")-1L
	xc <- Filter(function(x) nchar(x) != 0, rivet_out[xb:xe])
	xc <- unname(sapply(xc, function(num) eval(parse(text=num))))
	
	## Get y-grades 
	yb <- which(rivet_out == "y-grades")+1L
	ye <- which(rivet_out == "Dimensions > 0:")-1L
	yc <- Filter(function(x) nchar(x) != 0, rivet_out[yb:ye])
	yc <- unname(sapply(yc, function(num) eval(parse(text=num))))
	
	## Get Hilbert function values
	hb <- which(rivet_out == "Dimensions > 0:")+1L
	he <- which(rivet_out == "Betti numbers:")-1L
	hc <- Filter(function(x) nchar(x) > 0, rivet_out[hb:he])
	hx <- strsplit(gsub(pattern = "\\((\\d+), (\\d+), (\\d+)\\)", replacement = "\\1 \\2 \\3 ", x = hc), " ")
	hc <- do.call(rbind, lapply(hx, as.numeric))
	hc <- cbind(xc[hc[,1]+1L], yc[hc[,2]+1L], hc[,3])
	colnames(hc) <- c("x", "y", "dim")
	
	## Get non-zero Betti-0 numbers
	b0b <- which(rivet_out == "xi_0:")+1L
	b0e <- which(rivet_out == "xi_1:")-1L
	if (b0b != b0e){
		betti0_index <- t(str_to_simplex(rivet_out[b0b:b0e]))
		betti0_value <- cbind(x=xc[betti0_index[,1]+1], y=yc[betti0_index[,2]+1], value=betti0_index[,3])
	} else {
		betti0_index <- NULL
		betti0_value <- NULL
	}
	
	## Get non-zero Betti-1 numbers
	b1b <- which(rivet_out == "xi_1:")+1L
	b1e <- which(rivet_out == "xi_2:")-1L
	if (b1b != b1e){
		betti1_index <- t(str_to_simplex(rivet_out[b1b:b1e]))
		betti1_value <- cbind(x=xc[betti1_index[,1]+1], y=yc[betti1_index[,2]+1], value=betti1_index[,3])
	} else {
		betti1_index <- NULL
		betti1_value <- NULL
	}
	
	## Get non-zero Betti-2 numbers
	b2b <- which(rivet_out == "xi_2:")+1L
	if (is.numeric(b2b) && b2b != (length(rivet_out)+1L)){
		betti2_index <- t(str_to_simplex(rivet_out[b2b:length(rivet_out)]))
		betti2_value <- cbind(x=xc[betti2_index[,1]+1], y=yc[betti2_index[,2]+1], value=betti2_index[,3])
	} else {
		betti2_index <- NULL
		betti2_value <- NULL
	}
	
	out <- list(hf=hc, x_grades=xc, y_grades=yc)
	if (value){
		out <- modifyList(out, list(betti0=betti0_value, betti1=betti1_value, betti2=betti2_value))
	} else {
		out <- modifyList(out, list(betti0=betti0_index, betti1=betti1_index, betti2=betti2_index))
	}
	return(out)
}

compute_anchors <- function(S){
	I <- combn(nrow(S), 2)
	A <- lapply(1:ncol(I), function(i){
		x <- I[,i]
		pt1 <- S[x[1],] 
    pt2 <- S[x[2],] 
    if (all(pt1 < pt2)){ return(NULL) }
    else if (any(pt1 < pt2) || (pt1[1] == pt2[1] || pt1[2] == pt2[2])){
    	return(apply(S[x,], 2, max))
    }
    return(NULL)
	})
	A <- unique(do.call(rbind, A))
	return(A)
}

## projects a series of points in 'x' onto a line given by y = mx + b
proj_L <- function(x, m, b){
	is_above <- (x[,1]*m + b) < x[,2]
	proj <- ifelse(rep(is_above, each = 2), c((x[,2]-b)/m, x[,2]), c(x[,1], m*x[,1]+b))
	proj <- matrix(proj, ncol = 2, byrow = FALSE)
	return(proj)
	# projection <- apply(x, 1, function(pt){
	# 	if (pt[2] == m*pt[1] + b){ return(pt) }
	# 	is_above <- pt[1]*m + b < pt[2]
	# 	ifelse(rep(is_above, 2), c((pt[2]-b)/m, pt[2]), c(pt[1], m*pt[1]+b))
	# })
	# return(t(projection))
}

## polygonizes a line arrangement, returning the 2-cells 
## L := (n x 2) matrix of lines y = ax + b
polygonize_la <- function(L, xlim=c(0,1), ylim=c(0,1)){
	# X <- apply(L, 1, function(ab){ (ylim[2]-ab[2])/ab[1] })
	# Y <- apply(L, 1, function(ab){ ab[1]*xlim[2] + ab[2] })
	# XY1 <- cbind(X, ylim[2])
	# XY2 <- cbind(xlim[2], Y)
	# XY3 <- cbind(-L[,2]/L[,1], ylim[1])
	# XY4 <- cbind(xlim[1], L[,2])
	# 
	# ## Collect valid points in the box
	# m <- nrow(L)
	# t <- sapply(seq(m), function(i){ (xlim[1] <= XY1[i,1] & XY1[i,1] <= xlim[2]) })
	# r <- sapply(seq(m), function(i){ (ylim[1] <= XY2[i,2] & XY2[i,2] <= ylim[2]) })
	# b <- sapply(seq(m), function(i){ (xlim[1] <= XY3[i,1] & XY3[i,1] <= xlim[2]) })
	# l <- sapply(seq(m), function(i){ (ylim[1] <= XY4[i,2] & XY4[i,2] <= ylim[2]) })
	# 
	# pc <- list(XY1, XY2, XY3, XY4)
	# sides <- list(t,r,b,l)
	# 
	# valid_idx <- combn(4,2)
	# valid_lines <- do.call(rbind, lapply(1:ncol(valid_idx), function(i){
	# 	x <- valid_idx[,i]
	# 	valid <- sides[[x[1]]] & sides[[x[2]]] 
	#   cbind(pc[[x[1]]][valid,], pc[[x[2]]][valid,])
	# }))
	# edges <- lapply(seq(nrow(valid_lines)), function(i){ 
	# 	sf::st_linestring(rbind(valid_lines[i,1:2], valid_lines[i,3:4])) 
	# })

	## Get segments intersecting given lines with boxed region
	internal_segments <- lapply(1:nrow(L), function(i){
    el <- L[i,]
    int_pts <- matrix(c(xlim[1], el[2], xlim[2], xlim[2]*el[1] + el[2]), nrow = 2, byrow = TRUE)
    sf::st_linestring(int_pts)
	})
	
	## Make constrained box
	inner <- list(
		sf::st_linestring(matrix(c(xlim[1], xlim[2], ylim[2], ylim[2]), ncol = 2)), 
		sf::st_linestring(matrix(c(xlim[2], xlim[2], ylim[2], ylim[1]), ncol = 2)),
		sf::st_linestring(matrix(c(xlim[2], xlim[1], ylim[1], ylim[1]), ncol = 2)),
		sf::st_linestring(matrix(c(xlim[1], xlim[1], ylim[1], ylim[2]), ncol = 2))
	)
	
	## Form all the cells fo the line arrangement
	arrangement <- sf::st_multilinestring(c(inner, internal_segments))
	arrangement <- sf::st_union(arrangement, by_feature = TRUE)
	polys <- sf::st_polygonize(arrangement)
	return(polys)
}

## Form the dual graph
# dual_graph <- function(arrangement){
# 	stopifnot("GEOMETRYCOLLECTION" %in% class(arrangement))
# 	stopifnot(all(sapply(arrangement, function(x) "POLYGON" %in% class(x))))
# 	
# 	## Get polygon boundary intersections
# 	poly_intersects <- combn(length(arrangement), 2, function(x){ 
# 		length(sf::st_intersects(arrangement[[x[1]]], arrangement[[x[2]]])[[1]]) > 0 
# 	})
# 	## Determine if intersections are point-wise or along an edge
# 	intersect_at_end <- apply(combn(length(arrangement), 2)[,poly_intersects], 2, function(x){
# 		int_type <- class(sf::st_intersection(arrangement[[x[1]]], arrangement[[x[2]]]))
# 		return("LINESTRING" %in% int_type)
# 	})
# 	## The edges of the dual graph are the pairs of faces that share an edge 
# 	edges <- combn(length(arrangement), 2)[,which(poly_intersects)[intersect_at_end]]
# 	return(edges)
# }

#' @export
plot_hilbert <- function(x, betti=c(0,1,2), show_grid=TRUE, highlight=NULL, ...){
	# Plot Hilbert function
	pt_col <- rgb(0,0,0,x$hf[,3]/max(x$hf[,3]))
	if ("new" %in% names(list(...)) && list(...)[["new"]]){
		points(x$hf[,1:2], col = adjustcolor("white", alpha.f = 0), ...)
	} else {
		x_step <- diff(x$x_grades)[1]
		y_step <- diff(x$y_grades)[1]
		params <- list(...)
		if (!"xlim" %in% names(params)){ params <- modifyList(params, list(xlim=range(x$hf[,1]) + c(-1,1)*x_step)) }
		if (!"ylim" %in% names(params)){ params <- modifyList(params, list(ylim=range(x$hf[,2]) + c(-1,1)*y_step)) }
		params <- modifyList(params, list(x=x$hf[,1:2], col= adjustcolor("white", alpha.f = 0)))
		do.call(plot, params)
	}

	## Each row of hc indicates the south-west corner of the grid cell to highlight
	if (show_grid){
		points(expand.grid(x$x_grades, x$y_grades), pch = 20, col = adjustcolor("black", alpha.f = 0.10))
	}
	alpha.hf <- seq(0.20, 0.85, length.out = max(x$hf[,3])+1)
	u_hf <- sort(unique(x$hf[,3]))
	x_step <- abs(diff(x$x_grades)[1])
	y_step <- abs(diff(x$y_grades)[1])
	for (i in seq(nrow(x$hf))){
		x_idx <- match(x$hf[i,1], x$x_grades)
		y_idx <- match(x$hf[i,2], x$y_grades)
		box_col  <- adjustcolor("black", alpha.f = alpha.hf[x$hf[i,3]])
		if (!missing(highlight) && is.numeric(highlight)){
			if (x$hf[i,3] %in% highlight){
				box_col  <- adjustcolor("blue", alpha.f = alpha.hf[x$hf[i,3]])
			}
		}
		rect(xleft = x$x_grades[x_idx], xright = x$x_grades[x_idx] + x_step, 
				 ybottom = x$y_grades[y_idx], ytop = x$y_grades[y_idx] + y_step, 
				 col = box_col, lwd = 0)
	}
	if (is.logical(betti)){
		if (as.logical(betti)){
			points(x$betti0[,1:2], col = adjustcolor("green", alpha.f = 0.75), pch = 20, cex = 2.0)
			points(x$betti1[,1:2], col = adjustcolor("red", alpha.f = 0.75), pch = 20, cex = 2.0)
			points(x$betti2[,1:2], col = adjustcolor("yellow", alpha.f = 0.75), pch = 20, cex = 2.0)
		}
	} else {
		if (0 %in% betti){ points(x$betti0[,1:2], col = adjustcolor("green", alpha.f = 0.75), pch = 20, cex = 2.0) }
		if (1 %in% betti){ points(x$betti1[,1:2], col = adjustcolor("red", alpha.f = 0.75), pch = 20, cex = 2.0) }
		if (2 %in% betti){ points(x$betti2[,1:2], col = adjustcolor("yellow", alpha.f = 0.75), pch = 20, cex = 2.0) }
	}
}
