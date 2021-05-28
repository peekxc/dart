#' simulate_vineyard
#' @description Simulates a vineyard of persistence diagrams
#' @param R a square reduced matrix
#' @param V a square upper-triangular matrix
#' @param S a (2 x m) matrix giving the indices of the transpositions.
#' @param f an (optional) function to execute for each decomposition. See details. 
#' @param D an (optional) boundary matrix \code{R} was derived from, used to verify the decomposition is correct. 
#' @export
simulate_vineyard <- function(R, V, S, f = NULL, D = NULL, check_valid=TRUE, progress=FALSE){
	stopifnot(is.matrix(S), nrow(S) == 2, all(S <= ncol(R)))
	stopifnot(nrow(R) == ncol(R), nrow(V) == ncol(V), ncol(R) == ncol(V))
	
	if (progress){ pb <- txtProgressBar(min = 0, max = ncol(S), style = 3) }
	
	for (column in seq(ncol(S))){
		i <- min(S[,column])
		j <- max(S[,column])
		status <- 0L
		
		Pos <- apply(R, 2, function(x) all(x == 0))
		Low <- apply(R, 2, phtools:::low_entry)
		RLow <- rep(0L, length(Low))
		RLow[Low[Low != 0]] <- which(Low != 0)
		
		if (Pos[i] && Pos[j]){
			if (V[i,j] != 0){ V[,j] <- xor(V[,i], V[,j]) }
			k <- RLow[i]
			l <- RLow[j]
			if (k != 0 && l != 0 && R[i,l] != 0){
				if (k < l){
					R <- phtools::permute_move(R, i, j)
					V <- phtools::permute_move(V, i, j)
					R[,l] <- xor(R[,k], R[,l])
					V[,l] <- xor(V[,k], V[,l])
					status <- 1
				} else {
					R <- phtools::permute_move(R, i, j)
					V <- phtools::permute_move(V, i, j)
					R[,k] <- xor(R[,l], R[,k])
					V[,k] <- xor(V[,l], V[,k])
					status <- 2
				}
			} else {
				## Case 1.2
				R <- phtools::permute_move(R, i, j)
				V <- phtools::permute_move(V, i, j) 
				status <- 3
			}
		} else if (!Pos[i] && !Pos[j]){
			## Case 2
			if (V[i,j] != 0){
				## Case 2.1
				if (Low[i] < Low[j]){
					## Case 2.1.1
					R[,j] <- xor(R[,i], R[,j])
					R <- phtools::permute_move(R, i, j)
					V[,j] <- xor(V[,i], V[,j])
					V <- phtools::permute_move(V, i, j)
					status <- 4
				} else {
					## Case 2.1.2
					R[,j] <- xor(R[,i], R[,j])
					R <- phtools::permute_move(R, i, j)
					R[,j] <- xor(R[,i], R[,j])
					V[,j] <- xor(V[,i], V[,j])
					V <- phtools::permute_move(V, i, j)
					V[,j] <- xor(V[,i], V[,j])
					status <- 5
				}
			} else {
				## Case 2.2
				R <- phtools::permute_move(R, i, j)
				V <- phtools::permute_move(V, i, j) 
				status <- 6
			}
		} else if (!Pos[i] && Pos[j]){
			## Case 3
			if (V[i,j] != 0){
				## Case 3.1
				R[,j] <- xor(R[,i], R[,j])
				R <- phtools::permute_move(R, i, j)
				R[,j] <- xor(R[,i], R[,j])
				V[,j] <- xor(V[,i], V[,j])
				V <- phtools::permute_move(V, i, j)
				V[,j] <- xor(V[,i], V[,j])
				status <- 7
			} else {
				## Case 3.2
				R <- phtools::permute_move(R, i, j)
				V <- phtools::permute_move(V, i, j)
				status <- 8
			}
		} else {
			## Case 4
			if (V[i,j] != 0){ V[,j] <- xor(V[,i], V[,j]) }
			R <- phtools::permute_move(R, i, j)
			V <- phtools::permute_move(V, i, j) 
			status <- 9
		}
		if (!missing(D)){ D <- phtools::permute_move(D, i, j) }
		if (check_valid){
			is_reduced <- pbgrad::is_reduced(R)
			is_triangular <- Matrix::isTriangular(V %% 2)
			if (!missing(D)){
				dv <- (D %*% V) %% 2
				is_correct <- (all(((dv - R) %% 2) == 0))
			} else { is_correct <- TRUE }
			if (any(!c(is_reduced, is_triangular, is_correct))){
				browser()
				stop("Invalid decomposition detected.")
			}
		}
		if (!missing(f) && !missing(D)){ f(R, V, D, status, column) } 
		else if (!missing(f)){ f(R, V, status, column) }
		if (progress){ setTxtProgressBar(pb, column) }
	}
	if (progress){ close(pb) }
	
	return(list(R=R,V=V))
}

# Local vineyard simulation
simulate_vineyard_local <- function(R1, V1, S, R2=NULL, V2=NULL, f = NULL, D1 = NULL, D2 = NULL, check_valid=TRUE, progress=FALSE){
	if (progress){ pb <- txtProgressBar(min = 0, max = ncol(S), style = 3) }

	for (column in seq(ncol(S))){
		i <- min(S[,column])
		j <- max(S[,column])
		status <- 0L

		## Maps col-indices (1-dim) to row-indices (0-dim)
		Pos <- apply(R1, 2, function(x) all(x == 0))
		Low <- apply(R1, 2, phtools:::low_entry) ## row-index in 0-simplices
		
		## Maps row-indices (1-dim) to col-indices (2-dim) 
		Low2 <- apply(R2, 2, phtools:::low_entry) 
		RLow <- rep(0L, nrow(R2))
		RLow[Low2[Low2 != 0]] <- which(Low2 != 0)   ## col-index in 1-simplices

		if (Pos[i] && Pos[j]){
			if (V1[i,j] != 0){ V1[,j] <- xor(V1[,i], V1[,j]) }
			k <- RLow[i]
			l <- RLow[j]
			if (k != 0 && l != 0 && R2[i,l] != 0){
				if (k < l){
					R1 <- permute_move(R1, i, j, dims = "cols")
					V1 <- permute_move(V1, i, j, dims = "both")
					R2 <- permute_move(R2, i, j, dims = "rows")
					R2[,l] <- xor(R2[,k], R2[,l])
					V2[,l] <- xor(V2[,k], V2[,l])
					status <- 1
				} else {
					R1 <- permute_move(R1, i, j, dims = "cols")
					V1 <- permute_move(V1, i, j, dims = "both")
					R2 <- permute_move(R2, i, j, dims = "rows")
					R2[,k] <- xor(R2[,l], R2[,k])
					V2[,k] <- xor(V2[,l], V2[,k])
					status <- 2
				}
			} else {
				## Case 1.2
				R1 <- permute_move(R1, i, j, dims = "cols")
				R2 <- permute_move(R2, i, j, dims = "rows")
				V1 <- permute_move(V1, i, j)
				status <- 3
			}
		} else if (!Pos[i] && !Pos[j]){
			## Case 2
			if (V1[i,j] != 0){
				## Case 2.1
				if (Low[i] < Low[j]){
					## Case 2.1.1
					R1[,j] <- xor(R1[,i], R1[,j])
					R1 <- permute_move(R1, i, j, dims = "cols")
					R2 <- permute_move(R2, i, j, dims = "rows")
					V1[,j] <- xor(V1[,i], V1[,j])
					V1 <- permute_move(V1, i, j)
					status <- 4
				} else {
					## Case 2.1.2
					R1[,j] <- xor(R1[,i], R1[,j])
					R1 <- permute_move(R1, i, j, dims = "cols")
					R2 <- permute_move(R2, i, j, dims = "rows")
					R1[,j] <- xor(R1[,i], R1[,j])
					V1[,j] <- xor(V1[,i], V1[,j])
					V1 <- permute_move(V1, i, j)
					V1[,j] <- xor(V1[,i], V1[,j])
					status <- 5
				}
			} else {
				## Case 2.2
				R1 <- permute_move(R1, i, j, dims = "cols")
				R2 <- permute_move(R2, i, j, dims = "rows")
				V1 <- permute_move(V1, i, j)
				status <- 6
			}
		} else if (!Pos[i] && Pos[j]){
			## Case 3
			if (V1[i,j] != 0){
				## Case 3.1
				R1[,j] <- xor(R1[,i], R1[,j])
				R1 <- permute_move(R1, i, j, dims = "cols")
				R2 <- permute_move(R2, i, j, dims = "rows")
				R1[,j] <- xor(R1[,i], R1[,j])
				V1[,j] <- xor(V1[,i], V1[,j])
				V1 <- permute_move(V1, i, j)
				V1[,j] <- xor(V1[,i], V1[,j])
				status <- 7
			} else {
				## Case 3.2
				R1 <- permute_move(R1, i, j, dims = "cols")
				R2 <- permute_move(R2, i, j, dims = "rows")
				V1 <- permute_move(V1, i, j)
				status <- 8
			}
		} else {
			## Case 4
			if (V1[i,j] != 0){ V1[,j] <- xor(V1[,i], V1[,j]) }
			R1 <- permute_move(R1, i, j, dims = "cols")
			R2 <- permute_move(R2, i, j, dims = "rows")
			V1 <- permute_move(V1, i, j)
			status <- 9
		}
		if (!missing(D1)){ D1 <- permute_move(D1, i, j, dims = "cols") }
		if (!missing(D2)){ D2 <- permute_move(D2, i, j, dims = "rows") }
		if (check_valid){
			is_reduced_1 <- pbgrad::is_reduced(R1)
			is_reduced_2 <- pbgrad::is_reduced(R2)
			is_triangular_1 <- Matrix::isTriangular(V1 %% 2)
			is_triangular_2 <- Matrix::isTriangular(V2 %% 2)
			if (!missing(D1)){
				dv <- (D1 %*% V1) %% 2
				is_correct_1 <- (all(((dv - R1) %% 2) == 0))
			} else { is_correct_1 <- TRUE }
			if (!missing(D2)){
				dv <- (D2 %*% V2) %% 2
				is_correct_2 <- (all(((dv - R2) %% 2) == 0))
			} else { is_correct_2 <- TRUE }
			if (any(!c(is_reduced_1, is_triangular_1, is_correct_1, is_reduced_2, is_triangular_2, is_correct_2))){
				browser()
				stop("Invalid decomposition detected.")
			}
		}
		if (!missing(f) && !missing(D1)){ f(R1, R2, V1, V2, D1, D2, status, column) }
		else if (!missing(f)){ f(R1, R2, V1, V2, status, column) }
		if (progress){ setTxtProgressBar(pb, column) }
	}
	if (progress){ close(pb) }

	return(list(R1=R1,R2=R2,V1=V1,V2=V2))
}

#' vineyards_schedule
#' @description Creates a schedule of transpositions needs to transform one filtration into another. 
#' @details This function assumes K0 and K1 contain the same elements. 
#' @param K0 filtration
#' @param K1 filtation, as a list of simplices
#' @param w0 filtration grades for \code{K0}, in the same order as \code{K0}
#' @param w1 filtration grades for \code{K1}, in the same order as \code{K1}
#' @export
vineyards_schedule <- function(K0, K1, w0=NULL, w1=NULL, schedule_order="default"){
	stopifnot(is.list(K0), is.list(K1))
	if (missing(w0)){ w0 <- seq(length(K0)) }
	if (missing(w1)){ w1 <- seq(length(K1)) }
	stopifnot(length(w0) == length(K0), length(w1) == length(K1))
	stopifnot(all(K0 %in% K1), all(K1 %in% K0))
	
	## Ensure both filtrations are permutations of each other
	## If identity permutation, then return no schedule 
	k1_matching <- match(K1, K0)
	stopifnot(!any(is.na(k1_matching)))
	if (all(k1_matching == seq_along(K0))){ return(matrix(0, nrow = 2L, ncol = 0L)) }
	
	## Use asymptotically slower but in practice faster O(n^2) inversion algorithm 
	inv <- dart:::inversions(k1_matching-1L)
	stopifnot(ncol(inv) == kendall_dist(k1_matching))	
	
	## If wanted, attempt to form the homotopy, otherwise don't bother  and just report the inversions
	if (missing(schedule_order) || schedule_order %in% c("default", "insertion")){
		S <- dart:::relative_transpositions(length(K0), inv[,rev(seq(ncol(inv)))])+1L
		# is_valid_schedule <- all(apply(S, 2, diff) == 1L)
		is_valid_schedule <- all(diff(S)[seq(from = 1, to = ncol(S)-1L, by = 2)] == 1L)
		if (!is_valid_schedule){
			stop("Failed to accurately order transpositions due to rounding issues.")
		}
		return(list(schedule=S, lambda=NA))
	} else {
		int_points <- dart:::span_intersections(inv, w0, w1[match(K0, K1)], 0.0, 1.0)
		order_idx <- rev(order(int_points[1,], decreasing = TRUE))
		S <- dart:::relative_transpositions(length(K0), inv[,order_idx])+1L
		is_valid_schedule <- all(diff(S)[seq(from = 1, to = ncol(S)-1L, by = 2)] == 1L)
		if (!is_valid_schedule){
			stop("Failed to accurately order transpositions according to linear homotopy due to rounding issues.")
		}
		return(list(schedule=S, lambda=int_points[1,order_idx]))
	}
}



	## Converts transpositions by id into transpositions by position
	# relative_tr <- function(p, transpositions){
	# 	S <- matrix(0L, ncol = 0, nrow = 2)
	# 	for (i in seq(ncol(transpositions))){
	# 		tr_idx <- match(transpositions[,i], p)
	# 		S <- cbind(S, tr_idx)
	# 		p[tr_idx] <- p[rev(tr_idx)]
	# 	}
	# 	return(unname(S))
	# }
	# 
	## Magic to reverse the order in which tranpositions occur 
	
# }
	
	
	# order_idx <- order(int_points[1,])
	# inv_ordered <- inv[,order_idx]
	# dup_values <- unique(int_points[1,duplicated(int_points[1,])])
	# for (ux in dup_values){
	# 	idx <- which(int_points[1,] == ux)
	# 	inv_ordered[,idx] <- inv_ordered[,rev(idx)]
	# }

	# k1_matching <- match(K1, K0)
	# inv <- dart:::inversions(k1_matching)
	# int_points <- dart:::span_intersections(inv-1L, w0, w1, 0.0, 1.0)
	# S <- relative_tr(seq_along(K0), inv_ordered)
	# is_valid_schedule <- all(apply(relative_tr(P, inv_ordered), 2, diff) == 1L)
	# stopifnot(is_valid_schedule)
	
	## Generate the schedule
	# P <- seq_along(K0)
	# Q <- k1_matching
	# S <- generate_schedule(P, Q, inv[,order_idx], check = TRUE)

	# show_tr <- function(p, tr){
	# 	cat(p)
	# 	cat("\n")
	# 	for (i in seq(ncol(tr))){
	# 		p[tr[,i]] <- p[rev(tr[,i])]
	# 		cat(ifelse(p %in% p[tr[,i]] , crayon::black(p), crayon::red(p)))
	# 		cat("\n")
	# 	}
	# }
	# show_tr(Q, relative_tr(Q, inv)[,order_idx])
# }

	## Given two permutations 'p' and 'q' and a set of m indices to transpose between then, this 
	## function returns a (2 x m) matrix of relative positions which successively describe the 
	## positions of the elements to transpose in transforming p -> q
	# generate_schedule <- function(p, q, transpositions){
	# 	stopifnot(nrow(transpositions) == 2, is.matrix(transpositions), all(p %in% q))
	# 	stopifnot(all(as.vector(transpositions[,1:2]) %in% p))
	# 	S <- matrix(0L, ncol = 0, nrow = 2)
	# 	for (i in seq(nrow(transpositions))){
	# 		tr_idx <- match(transpositions[i,], p)
	# 		# stopifnot(abs(diff(tr_idx)) == 1)
	# 		S <- cbind(S, tr_idx)
	# 		p <- permute_move(p, i = tr_idx[1], j = tr_idx[2])
	# 	}
	# 	return(unname(S))
	# }
	# impose_total_order <- function(w, eps=.Machine$double.eps){
	# 	cfactor <- seq(1, 2, length.out = length(w))*eps
	# 	f_diffs <- c(diff(cfactor),0)
	# 	noise <- cfactor + runif(n = length(w), min = f_diffs*(1/3), max = f_diffs*(2/3))
	# 	stopifnot(all(order(noise) == seq_along(noise)))
	# 	return(w+noise)
	# }
	# 
	# f_k0 <- impose_total_order(w0, sqrt(.Machine$double.eps))
	# f_k1 <- impose_total_order(w1, sqrt(.Machine$double.eps))
	# 
	# ## F1 in K0 filtration order
	# f_k1_matched <- f_k1[matched_idx]
	
	## Find intersection points

	# indices <- dart:::SearchInStrip(w0, w1, 0, 1)
	# int_points <- dart:::span_intersections(indices, yb, ye, 0, 1)
	
	## Build a straight-line homotopy
	# sf_segms <- lapply(seq_along(K0), function(i){
	# 	sf::st_linestring(rbind(c(0, f_k0[i]), c(1, f_k1_matched[i])))
	# })
	# w1_matched <- w1[matched_idx]
	# sf_segms <- lapply(seq_along(K0), function(i){ 	
	# 	sf::st_linestring(rbind(c(0, w0[i]), c(1, w1_matched[i]))) 
	# })
	
	## Detect the crossings, ordered by x-axis

	# seg_sc <- sf::st_sfc(sf_segms)
	# seg_int <- sf::st_intersection(seg_sc, seg_sc)
	# intersections <- attr(seg_int, "idx")
	# crossings <- intersections[apply(intersections, 1, function(x) x[1] != x[2]),]
	# crossings <- unique(cbind(pmin(crossings[,1], crossings[,2]), pmax(crossings[,1], crossings[,2])))
	# int_pts <- apply(crossings, 1, function(idx){
	# 	sf::st_intersection(sf_segms[[idx[1]]], sf_segms[[idx[2]]])
	# })
	
	# stopifnot(length(unique(int_pts[1,])) == ncol(int_pts))
	
	## Enforce lex order on duplicates 
	# non_duplicates <- int_pts[1, !duplicated(int_pts[1,])]
	# dup_idx <- lapply(non_duplicates, function(x){ which(x == int_pts[1,]) })
	# for (d_idx in dup_idx){
	# 	if (length(d_idx) > 1){
	# 		crossings[d_idx,] <- kdtools::lex_sort(crossings[d_idx,])
	# 	}
	# }
# }
