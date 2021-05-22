library(phtools)
rips <- pbgrad::r_geometric_complex(14, radius = 0.20, filtered = TRUE)

P <- unlist(simplex_to_str(rips$simplices))
Q <- phtools:::sample_filtration(P)

p <- seq(sum(rips$n_simplices))
q <- sample(p)
lcs <- q[phtools:::longest_subseq(match(q,p))]
M <- greedy_min_cross(p, q, lcs)


RV <- pbgrad::reduce(pbgrad::boundary_matrix(rips))

double_mr <- sapply(1:(ncol(M)-1), function(i){ (M[1,i] < M[1,i+1]) & (M[2,i] < M[2,i+1]) & (M[2,i] > M[1,i+1]) })

M <- cbind(M[,1:3], c(18,29), c(28,45), c(28,43), M[,6:7], c(19,32), c(31,42), c(31, 36), M[,10:ncol(M)])

rt <- RV$R
vt <- RV$V
move_costs <- c()
for (i in seq(ncol(M))){
	if (M[1,i] < M[2,i]){
		res <- move_right(R = rt, V = vt, i = M[1,i], j = M[2,i])
	} else {
		res <- move_left(R = rt, V = vt, i = M[1,i], j = M[2,i])
	}
	move_costs <- c(move_costs, res$m)
	rt <- res$R
	vt <- res$V
}
sum(move_costs)