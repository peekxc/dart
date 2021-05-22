// #include <Rcpp.h>
// using namespace Rcpp;

// #include <stdlib>
// #include <span>
// // algorithm welzl is[8]
// //     input: Finite sets P and R of points in the plane |R|≤ 3.
// //     output: Minimal disk enclosing P with R on the boundary.
// // 
// //     if P is empty or |R| = 3 then
// //         return trivial(R)
// //     choose p in P (randomly and uniformly)
// //     D := welzl(P - { p }, R)
// //     if p is in D then
// //         return D
// // 
// //     return welzl(P - { p }, R ∪ { p })
// 
// using ball = std::pair< std::span< double >, double >;
// 
// template < typename Iter >
// auto welzl(Iter ps, Iter pe, Iter rs, Iter re, const size_t n, const size_t d) -> ball {
// 	const size_t pn = std::distance(ps, pe); 
// 	const size_t rn = std::distance(rs, re); 
// 	if (pn == 0 || rn == 3){
// 		
// 		return; 
// 	}
// 	size_t p_idx = rand() % pn;
// 	size_t p = *(ps + p_idx);
// 	pe = std::partition(ps, pe, [](size_t q) -> bool { return(q != p); })
// 	welzl(ps,pe,rs,re,n);
// 	
// }
// 
// // [[Rcpp::export]]
// NumericVector welzl_miniball(NumericMatrix x, const size_t idx) {
// 	if (idx > x.nrow()){ Rcpp::stop("invalid point"); }
// 	
//   return x * 2;
// }



/*** R

# algorithm welzl is[8]
#     input: Finite sets P and R of points in the plane |R|≤ 3.
#     output: Minimal disk enclosing P with R on the boundary.
# 
#     if P is empty or |R| = 3 then
#         return trivial(R)
#     choose p in P (randomly and uniformly)
#     D := welzl(P - { p }, R)
#     if p is in D then
#         return D
# 
#     return welzl(P - { p }, R ∪ { p })

rownames(P) <- seq(nrow(P))

euc_dist <- function(x,y){ dist(rbind(x,y)) }

circle_3pt <- function(x,y,z){
	A <- cbind(do.call(rbind, list(x,y,z)), 1)
	B <- cbind(apply(A[,1:2], 1, function(x) sum(x^2)), A[,2], 1.0)
	C <- cbind(B[,1], A[,1], 1.0)
	D <- cbind(C[,1], A[,1:2])
	ad <- det(A)
	bd <- -det(B)
	cd <- det(C)
	dd <- det(D)
	point <- c(-bd/(2*ad), -cd/(2*ad))
	c(point, radius=sqrt((bd^2 + cd^2 - 4*ad*dd)/(4*ad^2)))
}

smallest_enclosing_ball <- function(P, R){
	# if (class(P) == "numeric" || class(R) == "numeric"){
	# 	print(P)
	# 	print(R)
	# }
	if (nrow(P) == 0 || nrow(R) <= 3){
		if (nrow(R) == 1){ return(list(center=R, radius=0))}
		else if (nrow(R) == 2){ 
			return(list(center=colSums(R)/2, radius=max(dist(R))))
		} else if (nrow(R) == 3){
			max_pw <- max(dist(R))
			print("here")
			return(list(center=circle_3pt(R[1,],R[2,],R[3,]), radius=max(dist(R))))
		} else { stop("invalid") }
	}
	p <- sample(1:nrow(P), size = 1L)
	D <- smallest_enclosing_ball(P[-p,,drop=FALSE], R)
	if (euc_dist(D$center, P[p,,]) <= D$radius){
		return(D)
	}
	return(smallest_enclosing_ball(P[-p,,drop=FALSE], rbind(R, P[p,,drop=FALSE])))
}
seb <- smallest_enclosing_ball(X, X[1:3,])

center <- colSums(X)/nrow(X)
plot(X, xlim = center[1] + c(-1,1)*seb$radius, ylim = center[2] + c(-1,1)*seb$radius)
plotrix::draw.circle(seb$center[1], seb$center[2], radius = seb$radius)
points(seb$center[1], seb$center[2], pch = 20, col = "red")
*/
