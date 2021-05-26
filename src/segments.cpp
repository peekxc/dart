// enum INTERSECTION_TYPE { COLINEAR_DISJOINT = 0, COLINEAR_INTERSECT = 1, PARALLEL_DISJOINT = 2, DISJOINT = 3};

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"

// [[Rcpp::plugins(cpp17)]]
#include <optional>
#include <cinttypes>
using std::optional;
using std::make_optional;
using std::nullopt;
using std::vector; 

using segment = arma::vec4; 
using point2 = arma::vec2; 
static double cross_2(point2 u, point2 v) { return u[0]*v[1] - u[1]*v[0]; };
static constexpr bool in_unit(double x){ return(x <= 1.0 && x >= 0.0); };
static const auto pt0 = arma::span(0,1);
static const auto pt1 = arma::span(2,3);
	
using point = std::pair< double, double >; 
using std::make_pair;
	
// Checks intersection of two segments defined by (x0,y0,x1,y1)
optional< point2 > intersect_segments(point2 p, point2 pr, point2 q, point2 qs)  {
	point2 r = pr - p, s = qs - q;
	double rs = cross_2(r,s); 
	double t = cross_2(q-p, s)/rs;
	if (rs != 0 && in_unit(t) && in_unit(cross_2(q-p, r)/rs)){
		return(make_optional(p+t*r));
	}
	return nullopt;
}		

// Overload to allow checking (s0, s1), where s0=(x0,y0), s1=(x1,y1)
optional< point2 > intersect_segments(segment s0, segment s1)  {
	return intersect_segments(s0(pt0), s0(pt1), s1(pt0), s1(pt1)); 
}	

// [[Rcpp::export]]
arma::mat pairwise_segment_intersections(const arma::mat& S, const arma::mat& L, bool one_based=false) {
	if (S.n_rows != 4 || L.n_rows != 4){ throw std::invalid_argument("Invalid segment matrix given."); }
	const size_t n = S.n_cols;
	const size_t m = L.n_cols; 
	auto P = arma::mat(4, 0);
	
	size_t ii = 0;
	for (size_t i = 0; i < n; ++i){
		const segment s = S.col(i);
		for (size_t j = 0; j < m; ++j){
			auto point = intersect_segments(s, L.col(j));
			if (point){
				auto p = *point;
				P.insert_cols(ii++, arma::vec({ double(i + one_based), double(j + one_based), p[0], p[1] }));
			}
		}
	}
	return(P);
}

// Input: Expects a matrix where each column represents a segment, specified as (x0,y0,x1,y1)
// Output: Returns a matrix where each column represents an intersection between two segments, 
// specified as (s0, s1, x, y), where (s0,s1) are 0-based indices of the segments
// [[Rcpp::export]]
arma::mat all_segment_intersections(const arma::mat& S, bool one_based=false) {
	if (S.n_rows != 4){ throw std::invalid_argument("Invalid segment matrix given."); }
	const size_t n = S.n_cols;
	auto P = arma::mat(4, 0);
	size_t m = 0;
	size_t cc = 0; 
	for_each_combination_idx(n, size_t(2), [&S, &P, &m, &cc, one_based](vector< size_t > idx){
		if (cc % 10000 == 0){ Rcpp::checkUserInterrupt(); Rcpp::Rcout << cc << std::endl; }
		const size_t i = idx[0], j = idx[1];
		auto point = intersect_segments(S.col(i), S.col(j));
		if (point){
			auto p = *point;
			P.insert_cols(m++, arma::vec({ double(i + one_based), double(j + one_based), p[0], p[1] }));
		}
		++cc;
	});
	return(P);
}

// Given an set of segments spanning the strip [b,e] ordered w.r.t b, partitions
// greedily partitions L into staircase segments Q and non-staircase segments L' = L / Q
// auto split_staircase(const arma::mat& L) -> std::pair< arma::mat, arma::mat > {
// 	arma::mat Q = arma::mat(L.n_rows, 0);
// 	arma::mat Lp = arma::mat(L.n_rows, 0);
// 	if (L.n_cols == 0){ return std::make_pair(Q, Lp); }
// 	Q.insert_cols(0, L.col(0));
// 	for (size_t i = 1; i < L.n_cols; ++i){
// 		auto pt = intersect_segments(Q.col(Q.n_cols-1), L.col(i));
// 		if (!pt){
// 			Q.insert_cols(Q.n_cols, L.col(i));
// 		} else {
// 			Lp.insert_cols(Lp.n_cols, L.col(i));
// 		}
// 	}
// 	return(std::make_pair(Q, Lp));
// }
// 
// // Partitions L into two submatrices [Q, Lp] where Q represents a 
// // set of maximal disjoint segments and Lp = L / Q
// auto split_staircase_inplace(arma::mat& L) -> arma::mat::col_iterator {
// 	if (L.n_cols == 0){ return L.begin(); }
// 	size_t c_tail = 0; 
// 	for (size_t i = 1; i < L.n_cols; ++i){
// 		auto pt = intersect_segments(L.col(c_tail), L.col(i));
// 		if (!pt){
// 			L.swap_cols(i, ++c_tail);
// 		} 
// 	}
// 	return L.begin_col(c_tail+1);
// }

// bins a sequence { x1, x2, ..., xn } of non-decreasing values
// into m+1 non-decreasing bins given separated by values { b1, b2, ..., bm }
template < typename It, typename It2,  typename OutputIt, typename Compare = std::less<> >
void bincode(It x_a, const It x_e, It2 b_a, const It2 b_e, OutputIt out, Compare f = Compare()){
	for(size_t i = 0; x_a != x_e; ++x_a, ++out){
		while (f(*b_a, *x_a) && b_a != b_e){ ++b_a; ++i; }
		*out = i;
	}
}

// [[Rcpp::export]]
Rcpp::IntegerVector bin_break(Rcpp::NumericVector x, Rcpp::NumericVector bins){
	Rcpp::IntegerVector binned(x.size()); 
	bincode(x.begin(), x.end(), bins.begin(), bins.end(), binned.begin());
	return binned;
}

using std::pair; 

struct segment_t {
	double yb; 
	double ye;
	size_t id; 
};

constexpr auto less_b = [](const segment_t& s1, const segment_t& s2){ return s1.yb < s2.yb; };
constexpr auto less_be = [](const segment_t& s1, const segment_t& s2){ return s1.yb < s2.ye; }; // Note first on lhs, second on rhs
constexpr auto less_e = [](const segment_t& s1, const segment_t& s2){ return s1.ye < s2.ye; };

// Given two segments (s1,s2) spanning [b, e], check if they intersect
constexpr bool intersect_span(const segment_t& s1, const segment_t& s2) noexcept {
	// if ((s1.yb < s2.yb) != (s1.ye < s2.ye)){ Rprintf("(%d: %.3g, %.3g) intersects (%d: %.3g, %.3g)\n", s1.id+1, s1.yb, s1.ye, s2.id+1, s2.yb, s2.ye); }
	return (s1.yb < s2.yb) != (s1.ye < s2.ye);
}	

// Accepts an iterator of pairs
template< typename It, typename OutputIt >
void SearchInStrip_r(It Lb, It Le, const double b, const double e, OutputIt out){
	//using segment_t = typename It::value_type; // pair< double, double >
	const size_t n_segments = std::distance(Lb, Le);
	
	// Split into stairs [Q, Lp] where the range [Q, Lp) contains the stair segments, 
	// while the range [Lp, Le) contains the non-stair segments
	segment_t tail = *Lb;  
	auto Qb = Lb; 
	size_t c = 1; 
	for (size_t i = 1; i < n_segments; ++i){
		auto s_it = Qb + i;
		if (!intersect_span(*s_it, tail)){
			tail = *s_it; 
			std::rotate(Lb + c, Qb + i, Qb + i + 1);
			++c;
		}
	}
	auto Lp = Qb + c; 
	
	// Convenient renaming
	const size_t n_stairs = std::distance(Qb, Lp);
	const size_t n_non_stairs = std::distance(Lp, Le);
	const auto Q_end = Lp;

	// Base case: Q == R == Lb -> Le
	if (Lp == Le){ return; }

	// O(|Q| + |S|) binning procedure
	auto b_loc = arma::uvec(n_non_stairs);
	bincode(Lp, Le, Qb, Q_end, b_loc.begin(), less_b);
	
	// Reports an intersection to out
	const auto save_out = [&out, &Qb, &Lp](const size_t i, const size_t j){
		*out = std::make_pair((*(Qb + i)).id, (*(Lp + j)).id);
		++out;
	}; 
		
	// Test Int(Q, Lp): Compare the non-stair segments w/ the stair segments pairwise		
	for (size_t si = 0; si < n_non_stairs; ++si){

		// Variables
		bool continue_checking = false; 
		const segment_t non_stair = *(Lp + si);
		size_t i = std::min(size_t(b_loc[si]), size_t(n_stairs - 1)), j; 

		// Check at the index
		if (intersect_span(*(Qb + i), non_stair)){ save_out(i, si); }
		
		// Check below 
		j = i;
		continue_checking = (j > 0) && intersect_span(*(Qb + (j-1)), non_stair);
		while (continue_checking){
			save_out(--j, si);
			continue_checking = (j > 0) && intersect_span(*(Qb + (j-1)), non_stair);
		}
		
		// Check above
		j = i; 
		continue_checking = (j < (n_stairs - 1)) && intersect_span(*(Qb + (j+1)), non_stair);
		while (continue_checking){
			save_out(++j, si);
			continue_checking = (j < (n_stairs - 1)) && intersect_span(*(Qb + (j+1)), non_stair);
		}
	}
	
	// Recurse on the non-stair segments
	SearchInStrip_r(Lp, Le, b, e, out);
	
	// (Lp_b -> Lb) now sorted as Rp
	std::inplace_merge(Qb, Q_end, Qb + n_segments, less_e);
	return; 
}

// Given y-values 'y_b' and 'y_e' of points at x = b and x = e, respectively, compute all segment intersections
// Precondition: 'y_b' values are sorted, segment i is defined by (b, y_b[i]) -> (e, y_e[i])
// [[Rcpp::export]]
arma::umat SearchInStrip(std::vector< double > y_b, std::vector< double > y_e, double b, double e, bool one_based = false){
	if (y_b.size() != y_e.size()){ throw std::invalid_argument("Invalid segment endpoints given."); }
	bool y_b_sorted = std::is_sorted(y_b.begin(), y_b.end());
	if (!y_b_sorted){ throw std::invalid_argument("Initial endpoints must be sorted."); }
	
	// Create the segments
	const size_t m = y_b.size(); 
	auto L = vector< segment_t >();
	for (size_t i = 0; i < m; ++i){
		L.push_back(segment_t({ y_b[i], y_e[i], i }));
	}
	
	// Apply Balaban's recursive algorithm
	auto int_out = vector< pair< size_t, size_t > >();
	SearchInStrip_r(L.begin(), L.end(), b, e, std::back_inserter(int_out));
	
	// Fill the results
	arma::umat A_out(2, int_out.size()); 
	size_t i = 0; 
	for (auto& p_int: int_out){
		A_out.col(i++) = arma::uvec({ arma::uword(p_int.first + one_based), arma::uword(p_int.second + one_based) });
	}
	return(A_out);
}

// [[Rcpp::export]]
arma::mat span_intersections(const arma::umat indices, const Rcpp::NumericVector& y_b, const Rcpp::NumericVector& y_e, double b, double e, bool one_based = false){
	auto coords = arma::mat(2, indices.n_cols);
	size_t i = 0; 
	const double inf = std::numeric_limits< double >::infinity();
	indices.each_col([&](const arma::uvec& idx){
		auto s0 = segment{ b, y_b[idx[0] + one_based], e,  y_e[idx[0] + one_based] };
		auto s1 = segment{ b, y_b[idx[1] + one_based], e,  y_e[idx[1] + one_based] };
		auto pt = intersect_segments(s0, s1);
		coords.col(i++) = pt ? *pt : point2{ inf, inf };
	});
	return(coords);
}


/*** R
# X <- matrix(rnorm(10*5), ncol = 10, nrow = 5)
# test_arma(X)
# x <- c(0, 0.0099, 0.012, 0.30, 0.6, 0.80,0.999, 1, 1.1)
# b <- c(0.00, 0.01, 0.50, 0.99, 1.00)
# bin_break(x, b)
# test_staircase(S)
#217
# for (jj in 1:100){
	# i <- 18
	set.seed(18)
	ns <- c(3, 3, 6, 4)
	# c(sample(runif(4)), sample(runif(3, min = 1, max = 2)), sample(runif(3, min = 2, max = 3)))
	fa <- unlist(sapply(seq(length(ns)), function(i){ runif(ns[i], min=i, max=i+1) }))
	fb <- unlist(sapply(seq(length(ns)), function(i){ runif(ns[i], min=i, max=i+1) }))
	N <- length(fa)
	
	order_idx <- order(fa)
	fa <- fa[order_idx]
	fb <- fb[order_idx]
	
	plot.default(NULL, xlim = c(-0.1, 1.1), ylim = range(c(fa,fb)))
	abline(v = 0, lty = 20)
	abline(v = 1, lty = 20)
	points(cbind(0, fa), pch = 21)
	points(cbind(1, fb), pch = 21)
	
	for (i in seq(length(fa))){
	  segments(x0 = 0, x1 = 1, y0 = fa[i], y1 = fb[i])
	}
	text(x = 0, y = fa, labels = seq(N), pos = 2)
	text(x = 1, y = fb, labels = seq(N), pos = 4)
	
	L <- rbind(fa, fb)
	colnames(L) <- seq(N)
	stairs <- split_staircase(L, 0, 1)
	for (i in seq(ncol(stairs$Q))){
	  segments(x0 = 0, x1 = 1, y0 = stairs$Q[1,i], y1 = stairs$Q[2,i], lwd = 2.1)
	}

	# segments(x0 = 0, x1 = 1, y0 = fa[order_idx][order_idx==2], y1 = fb[order_idx][order_idx==2], lwd = 2)
	# segments(x0 = 0, x1 = 1, y0 = fa[order_idx][order_idx==5], y1 = fb[order_idx][order_idx==5], lwd = 2)
	# segments(x0 = 0, x1 = 1, y0 = fa[order_idx][order_idx==10], y1 = fb[order_idx][order_idx==10], lwd = 2)

	
	int_res <- pbgrad:::SearchInStrip(L[1,], L[2,], b = 0.0, e = 1.0, one_based = TRUE)
	print(L)
	L_test <- matrix(as.integer(colnames(L))[int_res], nrow = 2)
	L_test <- rbind(pmin(L_test[1,], L_test[2,]), pmax(L_test[1,], L_test[2,]))
	L_test <- t(kdtools::lex_sort(t(L_test)))
	L_test <- t(unique(t(L_test)))
	
	S_truth <- all_segment_intersections(rbind(rbind(L[1,], 0.0), rbind(L[2,], 1.0)), one_based = TRUE)
	S_truth <- matrix(as.integer(colnames(L))[as.integer(S_truth[1:2,])], nrow = 2)
	S_truth <- rbind(pmin(S_truth[1,], S_truth[2,]), pmax(S_truth[1,], S_truth[2,]))
	S_truth <- t(kdtools::lex_sort(t(S_truth)))
	
	matches <- (ncol(L_test) == ncol(S_truth)) && all(L_test == S_truth)
# 	print(matches)
# 	print(jj)
# 	if (!matches){ break }
# }
*/


enum POINT_TYPE { LEFT_END = 0, RIGHT_END = 1, INT_PT = 2 };

// Events for the X-structure / scheduler
struct event {
	uint_fast64_t index; 
	double x; 
	double y; 
	POINT_TYPE type; 
	bool operator< (const event& e) const { return x < e.x ; }
};

// Active markers for the Y-structure / sweep-line
struct active {
	uint_fast64_t id; 
	double y;
	bool operator< (const active& a) const { return y < a.y; }
};

template< typename SortedContainer, class UnaryPredicate, typename T = typename SortedContainer::value_type >
auto predecessor(const SortedContainer& c, UnaryPredicate f) -> optional< T > {
	auto elem_it = std::find_if(c.begin(), c.end(), f); 
	optional< T > pred;
	if (elem_it != c.end() && elem_it != c.begin()){
		pred.emplace(*std::prev(elem_it, 1));
	}
	return(pred);
}

template< typename SortedContainer, class UnaryPredicate, typename T = typename SortedContainer::value_type >
auto successor(const SortedContainer& c, UnaryPredicate f) -> optional< T > {
	auto elem_it = std::find_if(c.begin(), c.end(), f); 
	optional< T > succ;
	if (elem_it != c.end()){
		auto next_it = std::next(elem_it,1);
		if (next_it != c.end()){ succ.emplace(*next_it); }
	}
	return(succ);
}

// 'ancestry' returns the predecessor and the successor of a given value
template< typename SortedContainer, class UnaryPredicate, typename T = typename SortedContainer::value_type >// typename It = typename SortedContainer::iterator
auto ancestry(const SortedContainer& c, UnaryPredicate f) -> pair< optional< T >, optional< T > > {
	auto elem_it = std::find_if(c.begin(), c.end(), f); 
	optional< T > pred; 
	optional< T > succ;
	if (elem_it != c.end()){
		if (elem_it != c.begin()){ pred.emplace(*std::prev(elem_it, 1)); }
		auto next_it = std::next(elem_it,1);
		if (next_it != c.end()){ succ.emplace(*next_it); }
	}
	return(std::make_pair(pred, succ));
}

// One-liner to erase from a set using a predicate
// https://stackoverflow.com/questions/24263259/c-stdseterase-with-stdremove-if
template <class T, class Comp, class Alloc, class Predicate>
void discard_if(std::set<T, Comp, Alloc>& c, Predicate pred) {
  for (auto it{c.begin()}, end{c.end()}; it != end; ) {
    if (pred(*it)) { it = c.erase(it); }
    else { ++it; }
  }
}

using seg_id_t = uint_fast64_t; 
using std::set; 

// Intersects the segments given by S(i), S(j). If they intersect past the sweeping-line, 
// then their intersection point is removed from the event scheduler. 
void intersect_and_discard(seg_id_t i, seg_id_t j, const arma::mat& S, set< event >& X, const double SL){
	auto int_pt = intersect_segments(S.col(i), S.col(j));
	if (int_pt && (*int_pt)[1] > SL){ // If they intersect right of sweep line
		const auto id = szudzik_pair(i, j);
		discard_if(X, [id](const event& x){ return(x.index == id && x.type == INT_PT); });
	}
}

// Intersects the segments given by S(i), S(j). If they intersect past the sweeping-line, 
// then their intersection point is removed from the event scheduler. 
void intersect_and_insert(seg_id_t i, seg_id_t j, const arma::mat& S, set< event >& X, const double SL = -std::numeric_limits<double>::infinity()){
	auto int_pt = intersect_segments(S.col(i), S.col(j)); 
	// if (int_pt){ Rcpp::Rcout << "found intersection between: " << i+1 << ", " << j+1 << ", past SL: " << ((*int_pt)[1] > SL) << std::endl; }
	if (int_pt && (*int_pt)[1] > SL){ // if intersection point exists past the sweep-line
		X.insert({ szudzik_pair(i, j), (*int_pt)[0], (*int_pt)[1], INT_PT });
	}
}


// [[Rcpp::export]]
arma::mat bentley_ottmann(const arma::mat& S){
	if (S.n_rows != 4){ throw std::invalid_argument("Invalid segment matrix given."); }
	
	auto X = std::set< event >();		// event scheduler 
	auto Y = std::set< active >();  // active segments
		
	// Initialize with x-coordinates
	uint_fast64_t i = 0; 
	S.each_col([&X, &i](const arma::vec& v){ 
		X.insert({ i, v[0], v[1], LEFT_END }); 
		X.insert({ i, v[2], v[3], RIGHT_END }); 
		++i;
	});
	const double NEG_INF = -std::numeric_limits<double>::infinity(); 
	double SL = NEG_INF; 
	
	// While the schedule is not empty, process the next point
	size_t ii = 0, m = 0; 
	auto R = arma::mat(4, 0);
	while (!X.empty()){
		Rcpp::checkUserInterrupt();
		// Rcpp::Rcout << "i: " << ii++ << std::endl;
		
		Rcpp::Rcout << "X size: " << X.size() << ", ";
		Rcpp::Rcout << "Y size: " << Y.size() << std::endl;
		
		auto sv = *X.begin(); // Get pt w/ minimum x-value 
		X.erase(X.begin()); 	// Remove from scheduler
		SL = sv.x;						// Move sweep line forward 
		
		Rcpp::Rcout << "Y= ";
		for (auto& as: Y){ Rcpp::Rcout << "(" << as.id << ", " << as.y << ") "; }
		Rcpp::Rcout << std::endl;
		
		switch(sv.type){
			case LEFT_END: {
				Rcpp::Rcout << "LEFT END: " << sv.index+1 << std::endl;
				// A.1
				Y.insert({ sv.index, sv.y }); 
				auto ps = ancestry(Y, [&sv](const active& seg){ return seg.id == sv.index; });
				
				// A.2
				if (ps.first){ intersect_and_insert((*ps.first).id, sv.index, S, X); }
				if (ps.second){ intersect_and_insert(sv.index, (*ps.second).id, S, X); }
				
				// A.3
				if (ps.first && ps.second){ intersect_and_discard((*ps.first).id, (*ps.second).id, S, X, SL); }
				break;
			}
			case RIGHT_END: {
				Rcpp::Rcout << "RIGHT END: " << sv.index+1 << std::endl;
				// B.1 
				auto ps = ancestry(Y, [&sv](const active& seg){ return seg.id == sv.index; });
				
				// B.2 
				discard_if(Y, [&sv](const active& s){ return(s.id == sv.index); });
				
				// B.3 
				if (ps.first && ps.second){ intersect_and_insert((*ps.first).id, (*ps.second).id, S, X, SL); }
				break; 
			} // RIGHT_END
			case INT_PT: {
				Rcpp::Rcout << "INT_PT: " << sv.index+1 << ", ";
				// C.1 
				auto L = szudzik_unpair(sv.index);
				auto L1 = L.second, L2 = L.first;
				R.insert_cols(m++, arma::vec({ double(L2), double(L1), sv.x, sv.y }));
				
				// C.2 
				auto L1_succ = successor(Y, [L1](const active& s){ return s.id == L1; });
				auto L2_pred = predecessor(Y, [L2](const active& s){ return s.id == L2; });
				
				// C.3.1
				auto L1_it = std::find_if(Y.begin(), Y.end(), [L1](const active& s){ return s.id == L1; });
				auto L2_it = std::find_if(Y.begin(), Y.end(), [L2](const active& s){ return s.id == L2; });
				if (L1_it != Y.end() && L2_it != Y.end()){
					active s1 = *L1_it, s2 = *L2_it;
					Y.erase(L1_it), Y.erase(L2_it); // L2 should still valid after erase
					
					// L2 < L1 by pairing assumption, swap so L1 < L2
					s1.y = sv.y - std::numeric_limits< double >::epsilon();
					s2.y = sv.y + std::numeric_limits< double >::epsilon();
					
					// Reinsert
					Y.insert(s1);
					Y.insert(s2);
				}
				
				// C.3.2: Handle cases where predecessor/successor intersect to the right of the sweep-line
				if (L1_succ){ intersect_and_insert(L2, (*L1_succ).id, S, X, SL); }
				if (L2_pred){ intersect_and_insert((*L2_pred).id, L1, S, X, SL); }
				
				// C.4
				if (L1_succ){ intersect_and_discard(L1, (*L1_succ).id, S, X, SL); }
				if (L2_pred){ intersect_and_discard((*L2_pred).id, L2, S, X, SL); }
				break; 
			} // INT_PT
		} // switch 
	} // while 
	
	return(R);
}


	// arma::vec cc1 = S.col(idx[0]);
		// arma::vec cc2 = S.col(idx[1]);
	// auto P = arma::mat(4, intersections.size());
	// size_t i = 0;
	// for (auto& ip: intersections){ P.col(i++) = ip; }

					// if (int_pt && (*int_pt)[1] > SL){ // If they intersect right of sweep line
					// 	const size_t id = szudzik_pair(pr, pq);
					// 	discard_if(X, [id](const event& x){ return(x.index == id && x.type == INT_PT); });
					// }
// double test_vec1(arma::vec x, arma::vec y){
// 	return test_vec2(arma::conv_to< arma::vec2 >::from(x), arma::conv_to< arma::vec2 >::from(y));
// }

// // (x0,y0,x0,y1)
// using segment_t = std::array< double, 4 >; 
// using point_t = std::array< double, 2 >; 
// 
// 
// NumericVector segment_intersect(segment_t&& s0, segment_t&& s1) {
// 	// dot <- function(x, y){ sum(x*y) }
// 	// vector_cross <- function(v, w){ v[1]*w[2] - v[2]*w[1] }
// 	const auto dot2 = [](point_t x, point_t y){ return std::inner_product(x.begin(), x.end(), y.begin(), 0.0); }
// 	const auto cross2 = [](point_t v, point_t w) -> double { v[0]*w[1] - v[1]*w[0]; };
// 	
// 	auto p = point_t { s0[0], s0[1] };
// 	auto pr = point_t { s0[2], s0[3] };
// 	auto r = point_t { pr[0] - p[0], pr[1] - p[1] };
// 	
// 	auto q = point_t { s1[0], s1[1] };
// 	auto qs = point_t { s1[2], s1[3] };
// 	auto s = point_t { qs[0] - q[0], qs[1] - q[1] };
// 	
// 	double rs = cross2(r,s)
// 	t <- cross2(q-p, s)/vector_cross(r,s)
// 	u <- vector_cross(q-p, r)/vector_cross(r,s)
// 	
// 	if (rs == 0 && vector_cross(q-p, r) == 0){
// 		t0 <- dot((q-p), r)/(dot(r,r))
// 		t1 <- t0 + dot(s,r)/(dot(r,r))
// 		rng <- sort(c(t0,t1))
// 		intersects_unit <- (rng[1] <= 1.0 && rng[1] >= 0.0) || (rng[2] <= 1.0 && rng[2] >= 0.0)
// 		## Do somethign for colinear and overlapping and colinear and disjoint
// 		return(NULL)
// 	} else if (rs == 0 && vector_cross(q-p, r) != 0){
// 		## do parallel and non-intersecting
// 		return(NULL)
// 	} else if (rs != 0 && (t >= 0.0 && t <= 1.0) && (u >= 0.0 && u <= 1.0)){
// 		# They intersect!
// 		intersection_pt <- p+t*r
// 		return(intersection_pt)
// 	} else {
// 		## Line segments not parallel and do not intersect
// 		return(NULL)
// 	}
// }
// 
// // If not found, idx_t points to the lower bound 
// if (idx_it == Y.end() || (idx_it->index != sv.index)){
// 
// // If it's not at the beginning, insert the predecessor
// if (idx_it != Y.begin()){
// 	size_t pr = std::prev(idx_it,1)->index;
// 	auto int_pt = intersect_segments(S(pt0, pr), S(pt1, pr), S(pt0, sv.index), S(pt1, sv.index));
// 	if (int_pt.first){
// 		X.insert({ to_natural_2(sv.index, pr, m), int_pt.second[0], INT_PT });
// 	}
// }
// // If it's not at the end, insert the successor
// if (idx_it != Y.end()){
// 	size_t pq = std::next(idx_it,1)->index;
// 	auto int_pt = intersect_segments(S(pt0, pq), S(pt1, pq), S(pt0, sv.index), S(pt1, sv.index));
// 	if (int_pt.first){
// 		X.insert({ to_natural_2(sv.index, pq, m), int_pt.second[0], INT_PT });
// 	}
// }
	// Rcpp::Rcout << "Stairs: " << std::endl; 
	// std::for_each(Qb, Lp, [](const segment_t& s){ Rprintf("(%d: %.3g, %.3g)", s.id+1, s.yb, s.ye); });
	// Rcpp::Rcout << std::endl; 
	// 
	// Rcpp::Rcout << "Non-stairs: " << std::endl;
	// std::for_each(Lp, Le, [](const segment_t& s){ Rprintf("(%d: %.3g, %.3g)", s.id+1, s.yb, s.ye); });
	// Rcpp::Rcout << std::endl; 
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// auto lower_it = Y.lower_bound(sv.index);
				// auto upper_it = Y.upper_bound(sv.index);
				// 
				// bool lower_exists = lower_it != Y.begin(); 
				// bool upper_exists = upper_it != Y.end();
				// 
				// if (lower_exists){
				// 	size_t pr = std::prev(lower_it, 1)->index;
				// 	auto int_pt = intersect_segments(S(pt0, pr), S(pt1, pr), S(pt0, sv.index), S(pt1, sv.index));
				// 	if (int_pt.first){
				// 		X.insert({ to_natural_2(sv.index, pr, m), int_pt.second[0], INT_PT });
				// 	}
				// }
				// if (upper_exists){
				// 	size_t pq = upper_it->index;
				// 	auto int_pt = intersect_segments(S(pt0, pq), S(pt1, pq), S(pt0, sv.index), S(pt1, sv.index));
				// 	if (int_pt.first){
				// 		X.insert({ to_natural_2(sv.index, pq, m), int_pt.second[0], INT_PT });
				// 	}
				// }
				
				// // A.3
				// if (lower_exists && upper_exists){
				// 	size_t pr = std::prev(lower_it, 1)->index;
				// 	size_t pq = upper_it->index;
				// 	auto int_pt = intersect_segments(S(pt0, pr), S(pt1, pr), S(pt0, pq), S(pt1, pq));
				// 	if (int_pt.first){
				// 		size_t id = to_natural_2(pr, pq, m);
				// 		auto it = std::find_if(X.begin(), X.end(), [id](seg_value& x){ return(x.index == id && x.type == INT_PT); })
				// 		if (it != X.end()){
				// 			X.erase(it);
				// 		}
				// 	}
				// }				// if (ps.first){ // predecessor exists
				// 	auto pr = (*ps.first).id; 
				// 	auto ip = intersect_segments(S.col(pr), S.col(sv.index)); 
				// 	if (ip){ // if intersection point exists 
				// 		pt int_pt = *ip;
				// 		X.insert({ szudzik_pair(pr, sv.index), int_pt[0], int_pt[1], INT_PT });
				// 	}
				// }
				// if (ps.second){ // successor exists
				// 	auto pq = (*ps.second).id; 
				// 	auto ip = intersect_segments(S.col(pq), S.col(sv.index)); 
				// 	if (ip){ // if intersection point exists 
				// 		pt int_pt = *ip;
				// 		X.insert({ szudzik_pair(sv.index, pq), int_pt[0], int_pt[1], INT_PT });
				// 	}
				// }

/*** R
# set.seed(219)
# ns <- c(3, 3, 4)
# # c(sample(runif(4)), sample(runif(3, min = 1, max = 2)), sample(runif(3, min = 2, max = 3)))
# fa <- unlist(sapply(seq(length(ns)), function(i){ runif(ns[i], min=i, max=i+1) }))
# fb <- unlist(sapply(seq(length(ns)), function(i){ runif(ns[i], min=i, max=i+1) }))
# N <- length(fa)
# 
# plot.default(NULL, xlim = c(-0.1, 1.1), ylim = range(c(fa,fb)))
# abline(v = 0, lty = 20)
# abline(v = 1, lty = 20)
# points(cbind(0, fa), pch = 21)
# points(cbind(1, fb), pch = 21)
# 
# for (i in seq(length(fa))){
#   segments(x0 = 0, x1 = 1, y0 = fa[i], y1 = fb[i])
# }
# text(x = 0, y = fa, labels = seq(N), pos = 2)
# text(x = 1, y = fb, labels = seq(N), pos = 4)
# 
# remove_duplicates <- function(x){
# 	idx <- order(x)
# 	for (i in idx){
# 		idx_to_rep <- which(x == x[i])
# 		if (length(idx_to_rep) > 1){
# 			dupped <- x[idx_to_rep[-1]]
# 			dupped <- dupped + .Machine$double.eps^(1/1.2)*length(dupped)
# 			x[idx_to_rep[-1]] <- dupped
# 		}
# 	}
# 	return(x)
# }
# 
# S <- rbind(rbind(0, fa), rbind(1, fb))
# S <- t(apply(S, 1, remove_duplicates))

## bentley_ottmann(S)
# P <- bentley_ottmann(S)
# points(t(P[3:4,]), col = "red")



*/
