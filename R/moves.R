## ---- Main interface -----


permute_decomp <- function(x, method=c("moves", "vineyards", "naive")){
	stopifnot(inherits(x, "fbm_decomp"))
	
}

#' Simulate Dynamic Persistent Homology 
#' @description Simulates persistent homology (PH) in dynamic settings
#' @param f function that takes no arguments and yields a filtration on invocation.
#' @param FUN optional function to inspect persistence-related information for each filtration. See details. 
#' @details This function enables the efficient computation of persistence over dynamic or time-varying settings.
#' Given a nullary function which returns a filtration upon each invocation, this function computes the 
#' persistent homology of each filtration in sequence, applying the user-supplied \code{FUN} to the resulting 
#' persistence diagram. The user may specify the type interpolation process to use between adjacent filtrations. 
#' Thus, conceptually this acts as an alternative to collecting all of the filtrations into a list \code{L} and executing 
#' \code{lapply(L, function(filtration){ FUN(ph(filtration)) }) }. \cr
#' \cr     
#' Unlike the simplistic \code{lapply} strategy discussed above, this function requires at most 2 filtrations be 
#' in memory at any given time. Moreover, when adjacent persistence diagrams are relatively "close", it can be more efficient 
#' to update the underlying R=DV decomposition over time instead of simply recomputing the decomposition independently for 
#' each filtration. Optionally, the user may also simulate persistence at intermediate or "interpolated" points as well---the 
#' so-called vineyards approach. 
simulate_dyph <- function(f, FUN = NULL){
	
}


## ---- Helper functions -----
## Helper function to return the lowest non-zero entry in a vector, or 0 otherwise
low_entry <- function(x){ p <- which(x != 0); ifelse(length(p) == 0, 0, tail(p, 1L)) }

## ---- Move Left -----
## Performs reduction on a given set of columns, performing left-to-right additions mod 2
## The reduction is done in order of the largest low entries of the given indices J
restore_left <- function(R, V, J){
	stopifnot(pbgrad::is_matrix(R), pbgrad::is_matrix(V))
	J <- sort(J) # J should be sorted
	additions <- matrix(integer(0L), nrow = 2L)
	if (length(J) > 1){
		low_J <- apply(R[,J], 2, low_entry)
		pair <- sort(head(which(low_J == max(low_J)), 2L)) ## consider changing max to max duplicate
		while(length(pair) == 2 && length(unique(low_J[pair])) == 1L){
			if (all(low_J[pair] == 0)){ break }
			si <- J[pair[1]]
			ti <- J[pair[2]]
			R[,ti] <- (R[,ti] + R[,si]) %% 2
			V[,ti] <- (V[,ti] + V[,si]) %% 2
			low_J[pair[2]] <- low_entry(R[,ti])
			pair <- sort(head(which(low_J == low_J[pair[2]]), 2L))
			additions <- cbind(additions, c(si, ti))
		}
	}
	return(list(R=R,V=V,I=additions))
}

## Helper function to 
move_left_rv <- function(R, V, i, j){
	stopifnot(pbgrad::is_matrix(R), pbgrad::is_matrix(V), i > j, i <= ncol(R), j <= ncol(R))
	
	## If the column of V is already all zeros, just perform the move permutation
	if (all(V[j:(i-1),i] == 0)){
		R <- permute_move(R, i = i, j = j, dims = "cols")
		V <- permute_move(V, i = i, j = j, dims = "both")
		return(list(R=R,V=V,m=0))
	}
	
	## Set the donor column 
	dR <- R[,i,drop=FALSE]
	dV <- V[,i,drop=FALSE]
	
	## Clear the ones in the jth column of V to make sure it stays upper-triangular
	I <- seq(j, i-1L) ## absolute 
	idx <- I[tail(which(V[I, i] != 0), 1L)]
	indices <- c()
	while(length(idx) != 0){
		indices <- c(indices, idx)
		R[,i] <- (R[,idx] + R[,i]) %% 2
		V[,i] <- (V[,idx] + V[,i]) %% 2
		idx <- I[tail(which(V[I, i] != 0), 1L)]
	}
	
	## Perform the move permutation
	V <- permute_move(V, i = i, j = j, dims = "both")
	R <- permute_move(R, i = i, j = j, dims = "cols")
	
	## Collect indices of left-to-right column operations performed on V[,i]
	J <- c(j, rev(indices + 1L)) # J <- J[order(apply(R[,J], 2, low_entry), decreasing = TRUE)]
	RV_res <- restore_left(R, V, J)
	
	## Donate the last column 
	# donor_idx <- tail(J, 1L)
	# RV_res$V[,donor_idx] <- matrix(permute_move(dV, i, j, dims = "rows"))
	# RV_res$R[,donor_idx] <- dR
	
	## Return the resulting matrices, as well as the indices cancelled in the column of V before the move
	# V_additions <- rbind(c(indices, donor_idx), i)
	
	## Count the number of column operations 
	num_additions <- length(indices) + ncol(RV_res$I)
	return(list(R=RV_res$R, V=RV_res$V, m=num_additions))
}

move_left_full <- function(R, V, i, j){
	stopifnot(i > j, i <= ncol(R), j <= ncol(R))
	# stopifnot(as.logical(Matrix::isTriangular(R, upper = TRUE)))
	stopifnot(as.logical(Matrix::isTriangular(V, upper = TRUE)))
	
	## Save the donor column 
	dR <- R[,i,drop=FALSE]
	dV <- V[,i,drop=FALSE]
	
	## Clear the ones in the ith column of V to make sure it stays upper-triangular
	I <- seq(j, i-1L) ## absolute 
	idx <- I[tail(which(V[I, i] != 0), 1L)]
	indices <- c()
	while(length(idx) != 0){
		indices <- c(indices, idx)
		R[,i] <- (R[,idx] + R[,i]) %% 2
		V[,i] <- (V[,idx] + V[,i]) %% 2
		idx <- I[tail(which(V[I, i] != 0), 1L)]
	}
	
	## Perform the move permutation
	V <- permute_move(V, i = i, j = j, dims = "both")
	R <- permute_move(R, i = i, j = j, dims = "both")
	
	## Indices affected by the operations
	J <- sort(c(j, rev(indices + 1L)))
	RV_res <- restore_left(R, V, J)
	
	## Donate the last column 
	if (length(J) > 1){
		donor_idx <- tail(J, 1L)
		RV_res$V[,donor_idx] <- matrix(permute_move(dV, i, j, dims = "rows"))
		RV_res$R[,donor_idx] <- matrix(permute_move(dR, i, j, dims = "rows"))
	}

	# V_additions <- rbind(c(indices, donor_idx), i)
	num_additions <- length(indices) + ncol(RV_res$I)
	return(list(R=RV_res$R, V=RV_res$V, m=num_additions))
}

#' move_left
#' @export
move_left <- function(R, V, i, j, dims = "all"){
	RV_is_matrix <- pbgrad::is_matrix(R) && pbgrad::is_matrix(V)
	if (RV_is_matrix){
		# R_upper <- as.logical(Matrix::isTriangular(R, upper = TRUE))
		# R_square <- length(unique(dim(R))) == 1
		# if (R_upper && R_square){ ## Full matrix case 
			## Either extract the submatrices and recurse or write the custom code
			
			## Full matrix code -- i, j are absolute indices
		if (missing(dims) || dims == "all"){
			res <- move_left_full(R, V, i, j)
			return(list(R=res$R, V=res$V, m=res$m))
				# res1 <- move_left_rv(R = R, V = V, i = i, j = j)
				# 
				# ## Compute the possible columns involved in the reduction involved R2
				# low_J <- apply(R, 2, low_entry)
				# k <- which(low_J == i)
				# if (length(k) == 0){
				# 	R <- permute_move(R, i = i, j = j, dims = "rows")
				# 	m <- c(r1 = ncol(res1$I) + ncol(res1$J), r2 = 0)
				# 	return(list(R=R, V=V, m=m))
				# } else {
				# 	## Apply the row permutation, then restore the columns affected
				# 	k_low <- (j-1) + low_entry(R[seq(j, i-1L),k])
				# 	J <- sort(c(k, which(low_J <= k_low)))
				# 	R <- permute_move(R, i = i, j = j, dims = "rows") 
				# 	res2 <- restore_left(R = R, V = V, J = J)
				# 	m <- c(r1 = ncol(res1$I) + ncol(res1$J), r2 = ncol(res2$I))
				# 	return(list(R=res2$R, V=res2$V, m=m))
				# }
		} else {									## individual case
			result <- move_left_rv(R = R, V = V, i = i, j = j)
			return(list(R=result$R, V=result$V, m=result$m))
		}
	} else if (is.list(R) && is.list(V)) {
		stopifnot(length(R) == 2, length(V) == 2)
		stopifnot(!is.null(dim(R[[1]])), !is.null(dim(R[[2]])), !is.null(dim(V[[1]])), !is.null(dim(V[[2]])))
		stopifnot(ncol(R[[1]]) == ncol(V[[1]]), ncol(R[[2]]) == ncol(V[[2]]))
		R1 <- R[[1]]; R2 <- R[[2]]
		V1 <- V[[1]]; V2 <- V[[2]]
		
		## Perform the move on (R1, V1)
		res1 <- move_left_rv(R = R1, V = V1, i = i, j = j)
		
		## Compute the possible columns involved in the reduction involved R2
		low_J <- apply(R2, 2, low_entry)
		k <- which(low_J == i)
		if (length(k) != 0){
			k_low <- (j-1) + low_entry(R2[seq(j, i-1L),k])
			J <- sort(c(k, which(low_J <= k_low))) ## can we do better?
			print(J)
			R2 <- permute_move(R2, i = i, j = j, dims = "rows") 
			res2 <- restore_left(R = R2, V = V2, J = J)
		} else {
			R2 <- permute_move(R2, i = i, j = j, dims = "rows") 
			res2 <- list(R=R2, V=V2, I=matrix(numeric(0L), ncol = 0))
		}
		
		## Returns the new matrices 
		num_additions <- res1$m + ncol(res2$I)
		return(list(R=list(R1=res1$R, R2=res2$R), V=list(V1=res1$V, V2=res2$V), m=num_additions))
	} else { stop("Unknown input") }
}

## ---- Move Right -----
restore_right <- function(R, V, I){
	{ d_low <- low_entry(R[,I[1]]); d_colR <- R[,I[1],drop=FALSE]; d_colV <- V[,I[1],drop=FALSE] }
	if (length(I) <= 1){ return(list(R=R, V=V, d_colR=d_colR, d_colV=d_colV)) }
	for (k in seq(2, length(I))){
		{ new_d_low <- low_entry(R[,I[k]]); new_d_colR <- R[,I[k],drop=FALSE]; new_d_colV <- V[,I[k],drop=FALSE] }
		R[,I[k]] <- (R[,I[k]] + d_colR) %% 2
		V[,I[k]] <- (V[,I[k]] + d_colV) %% 2 
		if (d_low > new_d_low){ 
			{ d_low <- new_d_low; d_colR <- new_d_colR; d_colV <- new_d_colV }
		}
	}
	return(list(R=R, V=V, d_colR=d_colR, d_colV=d_colV))
}

move_right_local <- function(R, V, i, j){
	stopifnot(is.list(R), is.list(V), length(R) == 2, length(V) == 2)
	stopifnot(i < j)
	R1 <- R[[1]]; R2 <- R[[2]]
	V1 <- V[[1]]; V2 <- V[[2]]
	I <- seq(i,j)[which(V1[i,seq(i,j)] != 0)] ## I should always have length >= 1
	
	low_J <- apply(R2, 2, low_entry) 
	J <- which(low_J >= i & low_J <= j)
	J <- J[R2[i,J] != 0]

	result1 <- restore_right(R1, V1, I)
	R1_new <- permute_move(result1$R, i, j, dims = "cols")
	V1_new <- permute_move(result1$V, i, j, dims = "both")
	
	if (length(J) > 0){
		result2 <- restore_right(R2, V2, J)
		R2_new <- permute_move(result2$R, i, j, dims = "rows")
		V2_new <- result2$V
	} else {
		R2_new <- permute_move(R2, i, j, dims = "rows")
		V2_new <- V2
	}
	
	## Apply the donor column 
	R1_new[,j] <- result1$d_colR
	V1_new[,j] <- permute_move(result1$d_colV, i, j, dims = "rows")
	
	## Count the number of operations (not including donor column)
	num_additions <- (length(I)-1L) + ifelse(length(J) <= 1, 0, length(J) - 1L)
	
	return(list(R=list(R1_new,R2_new), V=list(V1_new,V2_new), m=num_additions))
}

#' move_right
#' @export
move_right <- function(R, V, i, j, dims = "all"){
	if (i == j){ return(list(R=R,V=V,m=0L)) }
	stopifnot(i <= j)
	if (is.list(R)){
		stopifnot(length(R) == 2, length(V) == 2)
		stopifnot(!is.null(dim(R[[1]])), !is.null(dim(R[[2]])), !is.null(dim(V[[1]])), !is.null(dim(V[[2]])))
		stopifnot(ncol(R[[1]]) == ncol(V[[1]]), ncol(R[[2]]) == ncol(V[[2]]))
		RV <- move_right_local(R = R, V = V, i = i, j = j)
		return(RV)
	} else {
		stopifnot(ncol(R) == ncol(V), i < ncol(R), j <= ncol(R))
		I <- seq(i,j)[which(V[i,seq(i,j)] != 0)] ## length(I) >= 1
		stopifnot(length(I) > 0)
		
		## Solution 1 for obtaining J 
		low_J <- apply(R, 2, low_entry) 
		J <- which(low_J >= i & low_J <= j)
		J <- J[R[i,J] != 0]
		#J <- which(R[i,] != 0)
		#J <- intersect(J, with(list(low_j=apply(R[,J,drop=FALSE], 2, low_entry)), { J[low_j >= i & low_j <= j] }))
		
		## Restore R, V
		result1 <- restore_right(R, V, I)
		R_new <- result1$R
		V_new <- result1$V
		
		## Restore R_1, V_1
		if (length(J) > 1){
			result2 <- restore_right(R_new, V_new, J)
			R_new <- result2$R
			V_new <- result2$V
		}
		
		## Generate the permutation 
		Rh <- permute_move(R_new, i, j)
		Vh <- permute_move(V_new, i, j)
		
		## Replace newly moved columns with the donors and return 
		Rh[,j] <- permute_move(result1$d_colR, i, j, dims = "rows")
		Vh[,j] <- permute_move(result1$d_colV, i, j, dims = "rows")
		
		## Count the number of operations (not including donor column)
		num_additions <- (length(I)-1L) + ifelse(length(J) <= 1, 0, length(J) - 1L)
		
		return(list(R=Rh, V=Vh, m=num_additions))
	}
} 

## ---- Related functions -----

#' move_decomp
#' @description Moves a simplex in a boundary matrix decomposition.
#' @export
move_decomp <- function(R, V, i, j, dims = "all"){
	if (i == j){ return(list(R=R,V=V,m=0L)) }
	if (i < j){ return(move_right(R, V, i, j, dims)) } 
	if (i > j){ return(move_left(R, V, i, j, dims)) }
}

#' execute_schedule
#' @description Executes a given schedule of permutations applied to PH-related decomposition. 
#' @param R Reduced matrix. 
#' @param V Upper triangular matrix recording boundaries. 
#' @param S schedule
#' @export
execute_schedule <- function(R, V, S, f = NULL){
	if (!missing(f) && !is.null(f)){ stopifnot(is.function(f)) }
	f_supplied <- (!missing(f) && !is.null(f))
	for (j in 1L:ncol(S)){
		new_rv <- move_decomp(R, V, i = S[1,j], j = S[2,j])
		if (f_supplied){ f(new_rv) }
		R <- new_rv$R
		V <- new_rv$V
	}
}

bubble_sort <- function(A){
	P <- matrix(0, ncol = 0, nrow = 2)
	n <- length(A)
	while(!all(order(A) == seq(n))){
		for (i in seq(n-1)){
			if (A[i] > A[i+1]){
				A[c(i, i+1)] <- A[c(i+1, i)]
				P <- cbind(P, c(i, i+1))
			}
		}
	}
	return(P)
}

#' @param x list of permutations
#' @export 
schedule <- function(x, type = c("bubble", "insertion", "selection", "lcs")){
	matches <- lapply(seq(length(x)-1L), function(i){
		match_idx <- match(x[[i+1L]], x[[i]])
		return(match_idx)
	})
	do.call(cbind, lapply(matches, function(ms){ bubble_sort(ms) }))
}


#' move_pair
#' @description Constructs valid LCS edit operation 
#' @details This function constructs a single, particular type of edit operation (a \emph{move}), and is usually called recursively. 
#' Given a single character \code{symbol} and two strings \code{source} and \code{target} whose longest common subsequence (LCS) is 
#' is \code{lcs}, this function computes a single edit operation using heuristic \code{rule} such that 
#' applying the edit operation \code{source} brings it closer to \code{target}, and increases the size of the LCS
#' by 1. 
#' @return a pair (i,j) of indicies indicating the move.
#' @export 
move_pair <- function(symbol, source, target, lcs, 
											rule = c("earliest", "latest", "closest", "furthest", "all")){
	stopifnot(!(symbol %in% lcs))
	if (missing(rule) || !(rule %in% c("earliest", "latest", "closest", "furthest", "all"))){
		rule <- "earliest"
	}
	n <- length(source)
	
	## Augments the source and target strings to make computing the valid positions easier
	lcs <- c(0, lcs, n+1)
	target_aug <- c(0, target, n+1)
	t_lcs_idx <- match(lcs, target_aug)
	l1 <- target_aug[tail(t_lcs_idx[t_lcs_idx < match(symbol, target_aug)], 1L)]
	l2 <- target_aug[head(t_lcs_idx[t_lcs_idx > match(symbol, target_aug)], 1L)]
	source_aug <- c(0, source, n + 1)
	s_rng <- range(match(c(l1, l2), source_aug))-1L
	
	## Obtain a vector of possible positions to move the symbol to
	si <- match(symbol, source)
	ti <- match(symbol, target)
	if (length(s_rng) != 2 || length(si) != 1 || length(ti) != 1){ stop("Invalid move pair") }
	if (si < s_rng[1]){           ## moving right
		s_rng[2] <- s_rng[2] - 1L
	} else if (si > s_rng[2]){    ## moving left
		s_rng[1] <- s_rng[1] + 1L
	}
	if (s_rng[1] <= 0){ s_rng[1] <- 1 }
	if (s_rng[2] > n){ s_rng[2] <- n }
	possible_idx <- seq(s_rng[1], s_rng[2])

	## Apply the rule
	stopifnot(length(possible_idx) >= 1)
	if (rule == "all"){ return(unname(cbind(si, possible_idx))) }
	index <- switch (rule,
    "earliest" = min(possible_idx),
    "latest" = tail(possible_idx, 1L), #max(possible_idx), 
    "closest" = possible_idx[which.min(abs(possible_idx - si))],
    "furthest" = possible_idx[which.max(abs(possible_idx - si))], 
	)
	return(cbind(si, index))
}


#' Greedy schedule minimization 
#' @description Creates a schedule of move permutation using a greedy heuristic
#' @details Given two permutations, \code{p} and \code{q}
#' @return (2 x m) matrix of move permutations.
#' @export
greedy_min_cross <- function(p, q, lcs, opt=c("minimize", "maximize"), use_lcs=TRUE){
	symbols <- setdiff(p, lcs)
	stopifnot(all(c(symbols, lcs) %in% p))
	
	if (missing(opt) || opt == "minimize"){
		opt <- min
		ignore_val <- Inf
	} else {
		opt <- max
		ignore_val <- -Inf
	}
	
	conflicts <- function(interval, conflicts){ 
		stopifnot(length(interval) == 2)
		interval <- sort(interval)
		sum(conflicts >= interval[1] & conflicts <= interval[2]) 
	}
	## Returns the net cost of performing the move given by int1
	# interval_cost <- function(s, o){
	# 	# overlap <- intersect(seq(s[1],s[2]), seq(o[1],o[2])) 
	# 	# if (length(overlap) == 0){ return(0) }
	# 	min_s <- min(s)
	# 	max_s <- max(s)
	# 	min_o <- min(o)
	# 	max_o <- max(o)
	# 	disjoint <- ifelse(min_s <= min_o, max_s < min_o, max_s > max_o)
	# 	if (disjoint){ return(0) }
	# 	if (min_o > min_s && max_o < max_s){ return(0) }
	# 	if (min_s > min_o && max_s < max_o){ return(0) }
	# 	smr <- s[1] < s[2] ## s interval moving right
	# 	omr <- o[1] < o[2] ## o interval moving right
	# 	if (smr && omr){
	# 		if (s[1] < o[1] && s[2] < o[2]){ return(1) }
	# 		if (s[1] > o[1] && o[2] < s[2]){ return(-1) }
	# 		if (s[2] == o[2]){ return(-1)}
	# 		browser()
	# 		stop("Logic error: this shouldn't occur")
	# 	}
	# 	if (!smr && !omr){
	# 		if (o[2] < s[2] && o[1] < s[1]){ return(1) }
	# 		if (s[2] < o[2] && s[1] < o[1]){ return(-1) }
	# 		if (s[2] == o[2]){ return(-1) }
	# 		browser()
	# 		stop("Logic error: this shouldn't occur")
	# 	}
	# 	if (min_o >= s[1] && s[1] <= max_o){ return(-1) }
	# 	if (min_o >= s[2] && s[2] <= max_o){ return(1) }
	# 	browser()
	# 	stop("Logic error: this shouldn't occur")
	# }
	
	moves <- matrix(0, nrow = 2, ncol = 0)
	ms <- length(symbols); 
	pb <- txtProgressBar(min = 0, max = ms, style = 3)
	while(length(symbols) > 0){
		
		I <- sapply(symbols, function(cs){ 
			move_pair(cs, source = p, target = q, lcs = lcs, rule = "closest") 
		})
		
		## Append LCS 
		if (use_lcs){
			I <- cbind(I, rbind(match(lcs, p), match(lcs, q)))
		}
		
		## Count how many pairwise intervals change
		if (ncol(I) <= 1){
			s_cost <- 0
		} else {
			# o_diff <- sapply(setdiff(seq(ncol(I)), si), function(oi){ interval_cost(I[,si], I[,oi]) })
			# sum(o_diff)
			s_cost <- sapply(seq(ncol(I)), function(si){
				interval_cost_rcpp(I[,si], I[,-si])
			})
		}
		if (use_lcs){ s_cost[tail(seq(ncol(I)), length(lcs))] <- ignore_val }
		opt_idx <- which(s_cost == opt(s_cost))
		if (length(opt_idx) == 1){
			mi <- I[,opt_idx]
		} else {
			I_size <- apply(I[,opt_idx], 2, function(x){ abs(diff(x)) })
			mi <- I[,head(opt_idx[I_size == opt(I_size)], 1L)]
		}
		# mi <- I[,opt(s_cost)]
		moves <- cbind(moves, mi)
		lcs <- c(lcs, p[mi[1]])
		p <- permute_move(p, i = mi[1], j = mi[2])
		lcs <- lcs[order(match(lcs, p))]
		symbols <- setdiff(symbols, lcs)
		setTxtProgressBar(pb, value = ms-length(symbols))
	}
	close(pb)
	return(unname(moves))
}


#' Greedy Scheduling 
#' @export
greedy_min_cross2 <- function(
	x, y, lcs = perm_lcs(x,y), opt=c("minimize", "maximize"), 
	strategy=c("target displacement", "pairwise", "local spearman")
){
	d <- length(x) - length(lcs)
	M <- matrix(0L, nrow = 2, ncol = d)
	cc <- 1L
	Lx <- match(lcs, x)
	Ly <- match(lcs, y)
	## Given a displacement vector and a pair (i,j), returns the updated displacement distance form the move 
	move_cost <- function(i, j, disp_vector){
		unaffected <- sum(abs(disp_vector[-(i:j)]))
		dis_change <- ifelse(i < j, sum(abs(disp_vector[(i+1):j]-1L)), sum(abs(disp_vector[j:(i-1)]+1)))
		change_in_i <- abs(ifelse(i < j, disp_vector[i]+abs(i-j), disp_vector[i]-abs(i-j)))
		return(dis_change + unaffected + change_in_i)
	}
	x_to_move <- setdiff(x, lcs)
	pos_in_y <- match(x_to_move, y)
	while(any(x != y)){
		pos_in_x <- match(x_to_move, x)
		matched_y <- findInterval(pos_in_y, Ly)
		mr <- findInterval(pos_in_x, Lx) < matched_y
		lower_bounds <- matched_y
		targets <- sapply(seq_along(mr), function(ii){ ifelse(mr[ii], Lx[lower_bounds[ii]], Lx[lower_bounds[ii]+1L]) })
		moves <- rbind(pos_in_x, targets)
		
		## Pick strategy heuristic 
		if (!missing(strategy)){
			stopifnot(is.character(strategy))
			strat_num <- agrep(strategy, c("target displacement", "pairwise", "local spearman"))
		} else {
			strat_num <- 1
		}
		stopifnot(length(strat_num) > 0, strat_num %in% 1:3)
		
		## Compute move costs 
		if (strat_num == 1){
			mc <- local({
				displacement <- seq_along(x) - match(x, y)
				apply(moves, 2, function(m){ move_cost(m[1], m[2], displacement) })
			})
		} else if (strat_num == 2){
			mc <- pairwise_cost(moves)
			# mc <- local({
			# 	sapply(1:ncol(moves), function(i){ dart:::interval_cost_rcpp(moves[,i], moves[,-i]) })
			# })
		} else {
			mc <- apply(moves, 2, function(m){
	    	spearman_dist(x, permute_move(x, i = m[1], j = m[2]))
			})
		}
		
		## Choose greedily based on heuristic cost
		greedy_move <- ifelse(missing(opt) || opt == "minimize", which.min(mc), which.max(mc))
		i <- moves[1,greedy_move]
		j <- moves[2,greedy_move]
		
		## Perform the move, update the relevent variables
		x <- permute_move(x, i = i, j = j)
		lcs <- union(lcs, x[j])
		lcs <- lcs[order(match(lcs, x))]
		Lx <- match(lcs, x)
		Ly <- match(lcs, y)
		
		## Remove moved symbol 
		remove_idx <- match(x[j], x_to_move)
		x_to_move <- x_to_move[-remove_idx]
		pos_in_y <- pos_in_y[-remove_idx]
		
		lcs_is_ordered <- all(order(Lx) == order(Ly))
		stopifnot(lcs_is_ordered)
		M[,cc] <- c(i, j)
		cc <- cc + 1L
	}
	return(M)
}
	
	# # stopifnot(length(intersect(J, Ly)) == 0L)
		# for (ii in seq_along(x_to_move)){
		# 	k <- findInterval(pos_in_y[ii], Ly) # match(max(Ly[Ly < J[ii]]), Ly)
		# 	l <- k+1
		# 	ifelse(mr[ii], Lx[k], Lx[l])
		# 	## [k,l) if mr[ii] = TRUE  
		# 	# if (mr[ii]){
		# 	# 	possible_idx <- seq(Lx[k], ifelse(l > length(Lx), Lx[k], Lx[l]-1))
		# 	# } else {
		# 	# 	possible_idx <- seq(ifelse(k == 0, Lx[l], Lx[k]+1), Lx[l])
		# 	# }
		# 	print(possible_idx)
		# }
	# move_cost <- function(i, A){
	# 	j <- i - dis[i]
	# 	unaffected <- sum(abs(dis[-(i:j)]))
	# 	dis_change <- ifelse(dis[i] < 0, sum(abs(dis[(i+1):j]-1L)), sum(abs(dis[j:(i-1)]+1)))
	# 	return(dis_change + unaffected)
	# }
	# # xx <<- x
	# while (any(x != y)){
	# 	## Form displacement vector
	# 	dis <- seq_along(x) - match(x, y)
	# 	move_idx <- setdiff(seq_along(x), match(lcs, x))
	# 	dis_after_move <- sapply(move_idx, function(i){ move_cost(i, dis) })
	# 	
	# 	lcs_idx_y <- match(lcs, y)
	# 	valid_moves <- sapply(move_idx, function(mi){
	# 		# bin0 <- findInterval(match(x[mi], y), lcs_idx_y)
	# 		# bin1 <- findInterval(mi - dis[mi], lcs_idx_y)
	# 		if (dis[mi] == 0){ return(TRUE) }
	# 		## Moving left
	# 		if (dis[mi] > 0){
	# 			## Find the LCS bin where the displaced symbol would go in x
	# 			bin_x <- findInterval(mi - dis[mi], match(lcs, x), left.open = TRUE)
	# 			## Find the LCS bin where the symbol should be in y
	# 			bin_y <- findInterval(match(x[mi], y), match(lcs, y), left.open = TRUE)
	# 		} else {
	# 			bin_x <- findInterval(mi - dis[mi], match(lcs, x), left.open = FALSE)
	# 			## Find the LCS bin where the symbol should be in y
	# 			bin_y <- findInterval(match(x[mi], y), match(lcs, y), left.open = FALSE)
	# 		}
	# 		return(bin_x == bin_y)
	# 	})
	# 	i <- move_idx[valid_moves][which.min(dis_after_move[valid_moves])]
	# 	
	# 	# dis_after_move <- sapply(move_idx, function(i){ spearman_dist(permute_move(x, i = i, j = i - dis[i]), y) })
	# 	# i <- move_idx[which.min(dis_after_move)]
	# 	x <- permute_move(x, i = i, j = i - dis[i])
	# 	lcs <- union(lcs, x[i - dis[i]])
	# 	print(order(match(lcs, y)) == order(match(lcs, x)))
	# 	lcs <- lcs[order(match(lcs, x))]
	# 	M[,cc] <- c(i, i - dis[i])
	# 	cc <- cc + 1L
	# }
	
	# lcs_idx <- match(lcs, y)
	# i - dis[i]
	# 
# }

move_sequence_recursive <- function(symbols, source, target, lcs, vn, moves = matrix(0L, nrow = 2, ncol = 0), ordered=FALSE, rule = "all"){
	if (all(source == target)){ 
		ms <- list(final=source, moves=moves)
		assign(x = vn, value = append(.GlobalEnv[[vn]], list(ms)), envir = .GlobalEnv)
	} else if (length(symbols) == 0){
		print(source)
		print(target)
		stop("Invalid sequence")
	} else {
		if (any(symbols %in% lcs)){ browser() }
		if (!ordered){
			for (cs in symbols){
				mi <- move_pair(cs, source = source, target = target, lcs = lcs, rule = rule)
				if (any(source[mi[,1]] %in% lcs)){ browser() }
				for (ii in seq(nrow(mi))){
					new_moves <- cbind(moves, mi[ii,])
					new_s <- permute_move(source, i = mi[ii,1], j = mi[ii,2])
					new_lcs <- c(cs, lcs)[order(match(c(cs, lcs), new_s))]
					if (match(new_lcs, target) != sort(match(new_lcs, target))){
						browser()
						stop("Invalid operation")
					}
					move_sequence_recursive(setdiff(symbols, cs), new_s, target, new_lcs, vn, new_moves,
																	ordered, rule=rule)
				}
			}
		} else {
			cs <- head(symbols, 1)
			if (length(cs) == 0){ browser() }
			mi <- move_pair(cs, source = source, target = target, lcs = lcs, rule = rule) 
			mi <- matrix(mi, ncol = 2)
			if (any(source[mi[,1]] %in% lcs)){ browser() }
			for (ii in seq(nrow(mi))){
				new_moves <- cbind(moves, mi[ii,])
				new_s <- permute_move(source, i = mi[ii,1], j = mi[ii,2])
				new_lcs <- c(cs, lcs)[order(match(c(cs, lcs), new_s))]
				move_sequence_recursive(setdiff(symbols, cs), new_s, target, new_lcs, vn, new_moves,
																ordered, rule=rule)
			}
		}
	}
}

#' move_sequence
#' @description Recursively produces sequences of move permutations
#' @details Given a set \code{symbols} lying in the set difference \code{s} \ \code{lcs}, this function 
#' recursively constructs sequences of edit operations which transform \code{s} -> \code{t}. If \code{ordered}
#' is \code{TRUE}, the edit operations move \emph{symbols} in the order they were given, otherwise all move 
#' sequences using the given \code{rule} are generated recursively. \cr
#' \cr 
#' NOTE: When \code{ordered} is FALSE, the output of this function is exponential in the size of \code{symbols} 
#' and \code{s}, and thus may be intractable for even modest input sizes.  
#' @export 
move_sequence <- function(symbols, s, t, lcs, ordered=FALSE, rule="all"){
	stopifnot(length(intersect(symbols, lcs)) == 0)
	stopifnot(length(s) == length(t))
	random_string <- function(n = 1) {
	  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
	  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
	}
	varname <- random_string()
	.GlobalEnv[[varname]] <- list()
	move_sequence_recursive(symbols, s, t, lcs, varname, ordered=ordered, rule=rule)
	res <- .GlobalEnv[[varname]] 
	base::remove(list = c(varname), envir = .GlobalEnv)
	res <- Filter(function(x) !is.null(x), res)
	return(res)
}



# greedy_min_cross2 <- function(
# 	x, y, lcs = perm_lcs(x,y), opt=c("minimize", "maximize"), 
# 	strategy=c("target displacement", "pairwise", "local spearman")
# ){
# 	d <- length(x) - length(lcs)
# 	M <- matrix(0L, nrow = 2, ncol = d)
# 	cc <- 1L
# 	Lx <- match(lcs, x)
# 	Ly <- match(lcs, y)
# 	## Given a displacement vector and a pair (i,j), returns the updated displacement distance form the move 
# 	move_cost <- function(i, j, disp_vector){
# 		unaffected <- sum(abs(disp_vector[-(i:j)]))
# 		dis_change <- ifelse(i < j, sum(abs(disp_vector[(i+1):j]-1L)), sum(abs(disp_vector[j:(i-1)]+1)))
# 		change_in_i <- abs(ifelse(i < j, disp_vector[i]+abs(i-j), disp_vector[i]-abs(i-j)))
# 		return(dis_change + unaffected + change_in_i)
# 	}
# 	while(any(x != y)){
# 		x_to_move <- setdiff(x, lcs)
# 		pos_in_x <- match(x_to_move, x)
# 		pos_in_y <- match(x_to_move, y)
# 		mr <- findInterval(pos_in_x, Lx) < findInterval(pos_in_y, Ly)
# 		lower_bounds <- findInterval(pos_in_y, Ly)
# 		targets <- sapply(seq_along(mr), function(ii){ ifelse(mr[ii], Lx[lower_bounds[ii]], Lx[lower_bounds[ii]+1L]) })
# 		moves <- rbind(pos_in_x, targets)
# 		
# 		## Pick strategy heuristic 
# 		if (!missing(strategy)){
# 			stopifnot(is.character(strategy))
# 			strat_num <- agrep(strategy, c("target displacement", "pairwise", "local spearman"))
# 		} else {
# 			strat_num <- 1
# 		}
# 		stopifnot(length(strat_num) > 0, strat_num %in% 1:3)
# 		
# 		## Compute move costs 
# 		if (strat_num == 1){
# 			mc <- local({
# 				displacement <- seq_along(x) - match(x, y)
# 				apply(moves, 2, function(m){ move_cost(m[1], m[2], displacement) })
# 			})
# 		} else if (strat_num == 2){
# 			mc <- local({
# 				sapply(1:ncol(moves), function(i){ dart:::interval_cost_rcpp(moves[,i], moves[,-i]) })
# 			})
# 		} else {
# 			mc <- apply(moves, 2, function(m){
# 	    	spearman_dist(x, permute_move(x, i = m[1], j = m[2]))
# 			})
# 		}
# 		# mc <- apply(moves, 2, function(m){ x[m[1]] })
# 		
# 		## Choose greedily based on heuristic cost
# 		greedy_move <- ifelse(missing(opt) || opt == "minimize", which.min(mc), which.max(mc))
# 		i <- moves[1,greedy_move]
# 		j <- moves[2,greedy_move]
# 		
# 		## Perform the move, update the relevent variables
# 		x <- permute_move(x, i = i, j = j)
# 		lcs <- union(lcs, x[j])
# 		lcs <- lcs[order(match(lcs, x))]
# 		Lx <- match(lcs, x)
# 		Ly <- match(lcs, y)
# 		lcs_is_ordered <- all(order(Lx) == order(Ly))
# 		stopifnot(lcs_is_ordered)
# 		M[,cc] <- c(i, j)
# 		cc <- cc + 1L
# 	}
# 	return(M)
# }
