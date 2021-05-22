#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"
#include "utility/combinations.h"
#include "utility/discrete.h"
#include "combinadic.h"
#include <tuple> 
#include <cstdint>
#include <iterator>
#include <chrono>

template< class Iter, class Index >
auto cycle_leader_shortcut(Iter t, Index start){
  auto cur = t[start];
  while(cur != start){
    if (cur < start) return false;
    cur = t[cur];
  }
  return true;
}
template< class Iter, class Index >
void reverse_cycle(Iter t, Index start){
  auto cur = t[start];
  auto prev = start;
  while( cur != start ){
    auto next = t[cur];
    t[cur] = prev;
    prev = cur;
    cur = next;
  }
  t[start] = prev;
}

// O(n log n) in-place permutation inverse algorithm
// Based on: https://stackoverflow.com/questions/56603153/how-to-invert-a-permutation-represented-by-an-array-in-place
template< class Iter >
void inverse_permutation_n(Iter first, size_t n){
  for (size_t i = 0; i != n; ++i){
		if (cycle_leader_shortcut(first, i)){ reverse_cycle(first, i); };
  }
}


// Schwartzian transform
// From: https://devblogs.microsoft.com/oldnewthing/20170106-00/?p=95135
template<typename Iter, typename UnaryOperation, typename Compare>
void sort_by_with_caching(Iter first, Iter last, UnaryOperation op, Compare comp) {
	using Diff = typename std::iterator_traits<Iter>::difference_type;
	using T = typename std::iterator_traits<Iter>::value_type;
	using Key = decltype(op(std::declval<T>()));
	using Pair = std::pair<T, Key>;
	Diff length = std::distance(first, last);
	
	// Construct the "cache" by applying op to every T 
	std::vector< Pair > pairs;
	pairs.reserve(length);
	std::transform(first, last, std::back_inserter(pairs), [&](T& t) { 
		return std::make_pair(std::move(t), op(t)); 
	});
	
	// Now do the cached comparison-based sort  
	std::sort(pairs.begin(), pairs.end(), [&](const Pair& a, const Pair& b) { 
		return comp(a.second, b.second); 
	});
	// ... and move the resulting sorted keys back into the range
	std::transform(pairs.begin(), pairs.end(), first, [](Pair& p) { return std::move(p.first); });
}

// Applies a permutation (indices) in O(n) time in-place - Essentially performs cycle sort
// Based on: https://devblogs.microsoft.com/oldnewthing/20170104-00/?p=95115
template< typename Iter1, typename Iter2 >
void apply_permutation(Iter1 b, Iter1 e, Iter2 p) {
	using T = typename std::iterator_traits< Iter1 >::value_type;
	using Diff = typename std::iterator_traits< Iter2 >::value_type;
	const Diff length = std::distance(b, e);
	for (Diff i = 0; i < length; ++i) {
		Diff current = i;
		if (i != p[current]) {				// Don't process trivial cycles
			T element{std::move(b[i])};
			while (i != p[current]) {   // Sort the current non-trivial cycle
				Diff next = p[current];
				b[current] = std::move(b[next]);
				p[current] = current;
				current = next;
			}
			b[current] = std::move(element);
			p[current] = current;
		}
	}
}

template< typename Iter1, typename Iter2 >
void apply_reverse_permutation(Iter1 first, Iter1 last, Iter2 indices) {
	using T = typename std::iterator_traits<Iter1>::value_type;
	using Diff = typename std::iterator_traits<Iter2>::value_type;
	const Diff length = std::distance(first, last);
	for (Diff i = 0; i < length; ++i) {
		while (i != indices[i]) {
			Diff next = indices[i];
			if (next < 0 || next >= length) { throw std::range_error("Invalid index in permutation"); }
	    if (next == indices[next]) { throw std::range_error("Not a permutation"); }
			std::swap(first[i], first[next]);
			std::swap(indices[i], indices[next]);
		}
	}
}

// For comparing very close values
template< typename T > [[nodiscard]]
constexpr bool within_machine_prec(T d1, T d2) noexcept {
	static_assert(std::is_floating_point< T >::value, "Must be floating point value.");
	return(std::abs(d1 - d2) < std::numeric_limits< T >::epsilon());
}


// Implicit Filtration class
// Stores a low-dimensional simplexwise filtration of simplices using a light encoding (combinatorial number system). 
// Filtration grades are stored minimally. 
// After construction, allows O(d log m) access to individual simplex < ranks / grades / weights >, 
// or alternatively O(m) access to stream the full filtration
struct ImplicitFiltration {
	const size_t m; 									 // number of total simplices 
	const size_t n; 									 // number of vertices
	const size_t d; 									 // dimension of complex + 1
	const vector< double > grades;		 // filtration grades; size depends on complex
	const vector< idx_t > labels;			 // vertex original labels
	vector< uint_fast64_t > ranks;     // simplex ranks 
	vector< size_t > indexes;          // simplex positions
	vector< uint_least8_t > dims;			 // simplex dimensions
	vector< size_t > n_simplexes; 		 // number of simplices of each dimension
	vector< size_t > cum_ns; 
		
	enum class filt_t { LOWER_STAR, FLAG, GENERIC };
	filt_t type; 
	
	// Return types for () and []
	enum class elem_t { RANK, GRADE, SIMPLEX, ALL };
	
	// Rcpp constructor
	ImplicitFiltration(SEXP stree, NumericVector weights) 
		: ImplicitFiltration(*Rcpp::XPtr< SimplexTree >(stree), weights.begin(), weights.end()){}; 	
	
	template< typename Iter >
	ImplicitFiltration(SimplexTree& st, Iter w_begin, const Iter w_end) 
		: m(std::accumulate(st.n_simplexes.begin(), st.n_simplexes.end(), 0)),
    	n(st.n_simplexes.at(0)),
    	d(st.dimension()+1),
    	grades(w_begin, w_end), 
    	labels(st.get_vertices()), 
      n_simplexes(st.n_simplexes.begin(), st.n_simplexes.end()) {
		
		std::partial_sum(n_simplexes.begin(), n_simplexes.begin()+d, std::back_inserter(cum_ns));
		// TODO: check choose(m, d) < 2^64 - 1 for encoding 
		
		// Check the type of filtration to construct
		if (grades.size() == n){
			type = filt_t::LOWER_STAR; 
		} else if (grades.size() == st.n_simplexes.at(1)){
			type = filt_t::FLAG; 
		} else if (grades.size() == m){
			type = filt_t::GENERIC; 
		} else {
			Rprintf("Weights size: %d, n: %d, m: %d\n", grades.size(), n, m);
			throw std::invalid_argument("Filtration weights must match size of vertices, edges, or all simplices.");
		}
		auto start1 = std::chrono::steady_clock::now();
    
		// Populate initial ranks + dims + indexes in shortlex order 
		auto tr = st::level_order< true >(&st);
		ranks.reserve(m); 
		dims.reserve(m);
		st::traverse(tr, [this](node_ptr p, size_t d, simplex_t simplex){
			ranks.push_back(dart::lex_rank(reindex(std::span{simplex}), n));
			dims.push_back(static_cast< uint_least8_t >(d-1));
			return true; 
		});
		indexes = vector< size_t >(m);
		std::iota(indexes.begin(), indexes.end(), 0);
	   
	  auto end1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> ms1 = (end1-start1);
    Rcout << "initial traversal took: " << ms1.count()*1000 << " ms" << std::endl;
    
		// Construct the filtration 
		// Postcondition: indexes array is sorted in filtration order
		auto start2 = std::chrono::steady_clock::now();
		// if (type != filt_t::FLAG){ throw std::invalid_argument("Filtration weights must match size of vertices, edges, or all simplices."); }
		switch(type){
			case filt_t::LOWER_STAR:{ build_lower_star(); break; }
			case filt_t::FLAG: { build_flag(); break; }
			default: { build_generic(); break; }
		}
		auto end2 = std::chrono::steady_clock::now();
		std::chrono::duration<double> ms2 = (end2-start2);
    Rcout << "building complex took: " << ms2.count()*1000 << " ms" << std::endl;

		// Sort ranks and dims vectors according to filtration order
		auto start3 = std::chrono::steady_clock::now();
		{
			auto p = indexes;
			apply_permutation(ranks.begin(), ranks.end(), p.begin());
			p = indexes;
			apply_permutation(dims.begin(), dims.end(), p.begin());
		}
		auto end3 = std::chrono::steady_clock::now();
		std::chrono::duration<double> ms3 = (end3-start3);
    Rcout << "Applying permutations took: " << ms3.count()*1000 << " ms" << std::endl;

		// Invert indexes to get O(log n) lookup using shortlex order
		// ranks[indexes[i]] := i'th simplex in shortlex order
		auto start4 = std::chrono::steady_clock::now();
		// inverse_permutation_n(indexes.begin(), m);
		{
			auto ip = vector< size_t >(m);
			std::iota(ip.begin(), ip.end(), 0);
			apply_reverse_permutation(ip.begin(), ip.end(), indexes.begin());
			indexes = std::move(ip);
		}
		auto end4 = std::chrono::steady_clock::now();
		std::chrono::duration<double> ms4 = (end4-start4);
    Rcout << "Inverse permutations took: " << ms4.count()*1000 << " ms" << std::endl;
	}
	
	void build_lower_star(){
		
		// Given a simplex index 'i', compute its grade via its vertices maximum grades
		const auto grade = [this](size_t i){
			if (dims[i] == 0){ return(grades[i]); }
			auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
			return std::reduce(s.begin(), s.end(), 0.0, [this](auto i, auto j){ return std::max(grades.at(i), grades.at(j)); });
		};
	
		// Created the lexicographically-refined ordering
		std::sort(indexes.begin(), indexes.end(), [this, &grade](auto i, auto j){
			auto wi = grade(i), wj = grade(j);
			if (within_machine_prec(wi, wj)){ return(dims[i] == dims[j] ? ranks[i] < ranks[j] : dims[i] < dims[j]); }
			return(wi < wj);
		});
	}
	
	// Created the lexicographically-refined ordering for a generic filtration
	void build_generic(){
		std::sort(indexes.begin(), indexes.end(), [this](auto i, auto j){
			auto wi = grades[i], wj = grades[j];
			if (within_machine_prec(wi, wj)){ return(dims[i] == dims[j] ? ranks[i] < ranks[j] : dims[i] < dims[j]); }
			return(wi < wj);
		});
	}
	
	// Builds a filtration from a flag complex (from edge weights) 
	// Preconditions: 
	// 1. grades.size() == st.n_simplexes[1] 
	// 2. indexes, ranks, and dims vectors are populated in shortlex order
	void build_flag(){
		
		// Extract edge ranks since they define the flag complex 
		// Ranks will be in shortlex order, which is naturally sorted
		auto edge_ranks = vector< uint_fast64_t >(); 
		edge_ranks.reserve(n_simplexes.at(1));
		for (size_t i = 0; i < m; ++i){
			if (dims[i] == 1){ edge_ranks.push_back(ranks[i]); }
		}
		
		// Given an index 'i', obtains the diameter of the simplex at position 'i' in the shortlex order
		const auto diameter = [this, &edge_ranks](size_t i) -> double {
			if (dims[i] == 0){ return 0.0; }
			else if (dims[i] == 1){ 
				auto e_it = std::lower_bound(edge_ranks.begin(), edge_ranks.end(), ranks[i]);
				return(grades[std::distance(edge_ranks.begin(), e_it)]);
			} else {
				auto s = dart::lex_unrank(ranks[i], n, dims[i]+1); 
				double diam = 0.0;
				for_each_combination(s.begin(), s.begin() + 2, s.end(), [this, &diam, &edge_ranks](auto b, auto e){
					auto it = std::lower_bound(edge_ranks.begin(), edge_ranks.end(), dart::lex_rank_2(*b, *(b+1), n));
					diam = std::max(diam, grades[std::distance(edge_ranks.begin(), it)]);
					return false;
				});
				return(diam);
			}
		};
		
		// Cache diameters before sort
		auto diameters = vector< double >(); 
		diameters.reserve(m);
		for (size_t i = 0; i < m; ++i){ diameters.push_back(diameter(i)); }
			
		// Created the lexicographically-refined ordering
		std::sort(indexes.begin(), indexes.end(), [this, &diameters](auto i, auto j){
			auto wi = diameters[i], wj = diameters[j];
			if (within_machine_prec(wi, wj)){ return(dims[i] == dims[j] ? ranks[i] < ranks[j] : dims[i] < dims[j]); }
			return(wi < wj);
		});
	}
	
	// indexes := i -> j := ( short lex i ) -> filtration index j such that ranks[indexes[i]] < ranks[indexes[i+1]]
	double explicit_grade(const size_t i) const {
		if (i >= m){ throw std::invalid_argument("Index larger than filtration"); }
		if (type == filt_t::FLAG){
			if (dims[i] == 0){ return 0.0; }
			else if (dims[i] == 1){
				// Rcout << "Searching for: " << i << " with rank " << ranks[i] << std::endl;
				const auto lb = cum_ns.at(0), ub = cum_ns.at(1);
				const auto it = std::lower_bound(indexes.begin()+lb, indexes.begin()+ub, ranks.at(i), [this](auto j, auto r){
					return(ranks.at(j) < r);
				});
				if (it != indexes.begin()+ub){
					auto gr_idx = std::distance(indexes.begin()+lb, it);
					return(grades.at(gr_idx));
				} else {
					throw std::invalid_argument("Invalid index/dimension given. Edge rank not found.");
				}
			} else {
				auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
				const auto lb = cum_ns.at(0), ub = cum_ns.at(1);
				double diam = 0.0;
				const auto ib = indexes.begin()+lb;
				const auto ie = indexes.begin()+ub;
				for_each_combination(s.begin(), s.begin() + 2, s.end(), [this, &ie, &ib, &diam](auto b, auto e){
					auto it = std::lower_bound(ib, ie, dart::lex_rank_2(*b, *(b+1), n), [this](auto j, auto r){
						return(ranks[j] < r);
					});
					diam = std::max(diam, grades.at(std::distance(ib, it)));
					return false;
				});
				return(diam);
			}
		} else if (type == filt_t::LOWER_STAR){
			// Given a simplex index 'i', compute its grade via its vertices maximum grades
			if (dims[i] == 0){ return(grades[i]); }
			auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
			return std::reduce(s.begin(), s.end(), 0.0, [this](auto i, auto j){ return std::max(grades.at(i), grades.at(j)); });
	
		}
		return(-1.0);
	}
	
	template< typename T >
	void grading(std::span< T > G) const {
		if (G.size() != m){ throw std::invalid_argument("Output range must match filtration size."); }
		for (size_t i = 0; i < m; ++i){
			G[i] = explicit_grade(i);
		}
	}
	
	// Returns the k-simplices in filtration order (unranked)
	// template< typename OutputIt >
	// void k_simplices(const size_t k, OutputIt out) const {
	// 	if (k >= d){ return; }
	// 	dart::lex_unrank(ranks[k].begin(), ranks[k].end(), n, k+1, out);
	// } // todo: change to match below
	
	// T must accept a span< size_t > as a constructor
	template< typename T >
	void simplices(std::span< T > S) const {
		if (S.size() != m){ throw std::invalid_argument("Output range must match filtration size."); }
		for (size_t i = 0; i < m; ++i){
			auto sx = dart::lex_unrank(ranks[i], n, dims[i]+1);
			auto ss = relabel(std::span{sx});
			S[i] = T(ss.begin(), ss.end());
		}
	}
	
	
	// Looks for a simplex in the filtration. If found, returns the filtration index
	// auto find(std::span< size_t > simplex) const noexcept -> std::optional< size_t > {
	// 	if (simplex.size() > d || simplex.size() == 0){ return(std::nullopt); }
	// 	const size_t ki = simplex.size()-1; 
	// 	auto it = std::lower_bound(ranks[ki].begin(), ranks[ki].end(), dart::lex_rank(simplex, n));
	// 	return(it != ranks[ki].end() ? std::make_optional(std::distance(ranks[ki].begin(), it)) : std::nullopt);
	// }
		
	// Returns the 0-based (contiguous) index of a vertex label ( O(log n) time )
	auto vertex_index(idx_t label) const noexcept -> size_t {
		return std::distance(labels.begin(), std::lower_bound(labels.begin(), labels.end(), label));
	};
	
	template< typename Size_t > [[nodiscard]]
	auto relabel(std::span< Size_t > s) const -> std::span< Size_t > {
		std::transform(s.begin(), s.end(), s.begin(), [this](auto i){ return(this->labels[i]); });
		return(s);
	}
	
	template< typename Size_t > [[nodiscard]]
	auto reindex(std::span< Size_t > s) const -> std::span< Size_t > {
		std::transform(s.begin(), s.end(), s.begin(), [this](auto l){ return(this->vertex_index(l)); });
		return(s);
	}
};


SEXP as_XPtr(ImplicitFiltration* fi){
  Rcpp::XPtr< ImplicitFiltration > p(fi, false);
  return(p);
}

List fi_simplices(ImplicitFiltration* fi) {
	ImplicitFiltration& f = *fi; 
	auto S = vector< simplex_t >(f.m); 
	f.simplices(std::span{S});
	return(Rcpp::wrap(S));
}
NumericVector fi_grades(ImplicitFiltration* fi) {
	ImplicitFiltration& f = *fi; 
	auto G = vector< double >(f.m); 
	f.grading(std::span{G});
	return(Rcpp::wrap(G));
}

RCPP_MODULE(implicit_filtration_module) {
  using namespace Rcpp;
  Rcpp::class_< ImplicitFiltration >("ImplicitFiltration")
    .constructor< SEXP, NumericVector >()
    .method( "as_XPtr", &as_XPtr)
    .field_readonly("size", &ImplicitFiltration::m)
    .field_readonly("num_vertices", &ImplicitFiltration::n)
    .field_readonly("cum_sizes", &ImplicitFiltration::cum_ns)
    .field_readonly("minimal_grades", &ImplicitFiltration::grades)
    .field_readonly("vertex_labels", &ImplicitFiltration::labels)
    .field_readonly("simplex_ranks", &ImplicitFiltration::ranks)
    .field_readonly("simplex_dims", &ImplicitFiltration::dims)
    .field_readonly("shortlex_perm", &ImplicitFiltration::indexes)
    .property("simplices", &fi_simplices, "Returns the simplices in filtration order as a List.")
	  .property("grading", &fi_grades, "Returns the grading in filtration order as a vector.")
    // .method("build_flag", &build_flag)
    ;
}

// [[Rcpp::export]]
IntegerVector inverse_permutation(IntegerVector p){
	auto a = IntegerVector(p.size());
	std::iota(a.begin(), a.end(), 0);
	apply_reverse_permutation(a.begin(), a.end(), p.begin());
	return(a);
}

// [[Rcpp::export]]
IntegerVector inverse_permutation2(IntegerVector p){
	inverse_permutation_n(p.begin(), p.size());
	return(p);
}
// [[Rcpp::export]]
IntegerVector test_unrank(const size_t r, const size_t n, const size_t k){
	auto s = dart::lex_unrank(r, n, k);
	IntegerVector v(s.begin(), s.end());
	return(v);
};
 

/*** R
unrank_cc <- function(r, k, n){
	N <- choose(n, k)
	ri <- (N-1L)-r
  s <- vector("integer", length = k)
  C <- lapply(seq(k), function(ki){ choose(0:n, k = ki) })
  for (i in seq(length(s), 1L)){
  	s[i] <- findInterval(x = ri, vec = C[[k]])-1L
  	ri <- ri - C[[k]][s[i]+1L]
  	k <- k - 1L
  }
  s <- (n-1L) - s
  return(rev(s))
}

A <- combn(5,3)-1L
unrank_cc(r = 5, k = 3, n = 5)
dart:::test_unrank(r = 5, n = 5, k = 3)
r = 5; n = 5; k = 3

A <- combn(10,3)-1L
tests <- lapply(0:(ncol(A)-1L), function(r){
	dart:::test_unrank(r = r, n = 10, k = 3)
})
sum(abs(do.call(cbind, tests) - A)) == 0

*/

// int test_binom_search(const size_t r, const size_t n, const size_t k){
// 	int out = binary_search(r, n, [k](const auto& key, int index) -> int {
// 		auto c = dart::BinomialCoefficient(index, k);
// 		return(key == c ? 0 : (key < c ? -1 : 1));
// 	});
// 	return(out);
// } 


// IntegerVector test_unrank(const size_t r, const size_t n, const size_t k){
// 	const size_t N = dart::BinomialCoefficient(n, k);
// 	size_t ri = (N-1) - r; 
// 	auto S = vector< size_t >(k);
// 	for (size_t ki = k; ki > 0; --ki){
// 		S.at(ki-1) = binary_search(ri, n, [ki](const auto& key, int index) -> int {
// 			auto c = dart::BinomialCoefficient(index, ki);
// 			return(key == c ? 0 : (key < c ? -1 : 1));
// 		});
// 		Rprintf("ri: %d, ki: %d, s[i]: %d, bc: %d\n", ri, ki, S.at(ki-1), dart::BinomialCoefficient(S.at(ki-1), ki));
// 		ri -= dart::BinomialCoefficient(S.at(ki-1), ki); 
// 	}
// 	std::transform(S.begin(), S.end(), S.begin(), [n](auto el){ return((n-1) - el); });
// 	std::reverse(S.begin(), S.end());
// 	return(wrap(S));
// }


// int test_binary_search(const IntegerVector& a, int key){
// 	int out = binary_search(key, a.size(), [&a](const int k, int idx) -> int {
// 		return(k < a.at(idx) ? -1 : (k > a.at(idx) ? 1 : 0));
// 	});
// 	return(out);
// }
// for (auto it = indexes.begin()+lb; it != indexes.begin()+ub; ++it){
// 	Rcout << *it << ", ";
// }
// Rcout << "Ranks: "; 
// // for (auto it = indexes.begin()+lb; it != indexes.begin()+ub; ++it){
// // 	Rcout << ranks.at(*it) << ", ";
// // }
// 
// vector< uint_least64_t > ER = vector< uint_least64_t >();
// for (auto it = indexes.begin()+lb; it != indexes.begin()+ub; ++it){
// 	ER.push_back(ranks.at(*it));
// }
// for (auto er: ER){ Rcout << er << ", "; }
// Rcout << std::endl;
// Rcout << "is sorted? " << std::is_sorted(ER.begin(), ER.end()) << std::endl;
// // const auto it = std::lower_bound(ER.begin(), ER.end(), ranks.at(i), [this](auto j1, auto j2){
// // 	return(ranks.at(j1) < ranks.at(j2));
// // });
// for (size_t ki = 0; ki < d; ++ki){
// 	auto& ix = indexes[ki];
// 	auto& rk = ranks[ki];
// 	for (size_t j = 0; j < ix.size(); ++j){
// 		auto sx = relabel(dart::lex_unrank(rk[j], n, ki+1));
// 		S[ix[j]] = T(sx);
// 	}
// }
// 
	// template< elem_t ET >
	// constexpr auto filt_at(const size_t i) const {
	// 	auto k = std::optional< size_t >();
	// 	auto relative_idx = std::optional< size_t >(); 
	// 	for (size_t ki = 0; ki < d; ++ki){
	// 		auto& idx = indexes[ki];
	// 		auto it = std::lower_bound(idx.begin(), idx.end(), i);
	// 		if (it != idx.end()){
	// 			k = std::make_optional(ki); 
	// 			relative_idx = std::make_optional(std::distance(it, idx.end()));
	// 			break; 
	// 		}
	// 	}
	// 	if (!k || !relative_idx){ return std::nullopt;	}
	// 	if constexpr (ET == RANK){
	// 		auto r = ranks[k.value()][relative_idx.value()]; 
	// 		return(std::make_pair(r, k.value()));
	// 	} else if constexpr (ET == GRADE){
	// 		
	// 	} else if constexpr (ET == SIMPLEX){
	// 		auto r = ranks[k.value()][relative_idx.value()]; 
	// 		return(relabel(dart::lex_unrank(r, n, k.value()+1)));
	// 	} else {
	// 			
	// 		}
	// 	}
	// 	if constexpr (ET == RANK){
	// 		
	// 	} else if constexpr (ET == GRADE){
	// 		
	// 			if (it != idx.end()){
	// 				return relabel(dart::lex_unrank(ranks[ki][std::distance(it, idx.end())], n, ki+1));
	// 			}
	// 		}
	// 	} else if constexpr (ET == SIMPLEX){
	// 		for (size_t ki = 0; ki < d; ++ki){
	// 			auto& idx = indexes[ki];
	// 			auto it = std::lower_bound(idx.begin(), idx.end(), i);
	// 			if (it != idx.end()){
	// 				return relabel(dart::lex_unrank(ranks[ki][std::distance(it, idx.end())], n, ki+1));
	// 			}
	// 		}
	// 		return(std::span< size_t >());
	// 	} else {
	// 		
	// 	}
	// }