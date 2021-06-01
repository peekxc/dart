#include <Rcpp.h>
using namespace Rcpp;

#include "implicit_filtration.h"

SEXP as_XPtr(ImplicitFiltration* fi){
  Rcpp::XPtr< ImplicitFiltration > p(fi, false);
  return(p);
}

// Rcpp constructor
ImplicitFiltration* IF_constructor(SEXP stree, NumericVector weights){
	auto fi = new ImplicitFiltration(*Rcpp::XPtr< SimplexTree >(stree),  weights.begin(), weights.end());
	return(fi);
}

List fi_simplices(ImplicitFiltration* fi) {
	ImplicitFiltration& f = *fi; 
	auto S = vector< simplex_t >(f.m); 
	f.simplices(std::span{S});
	return(Rcpp::wrap(S));
}
List fi_k_simplices(ImplicitFiltration* fi, const size_t k) {
	ImplicitFiltration& f = *fi; 
	auto S = vector< simplex_t >(); 
	if (f.n_simplexes.size() > k){
		S.resize(f.n_simplexes.at(k));
		f.k_simplices(std::span{S}, k);
	}
	return(Rcpp::wrap(S));
}
NumericVector fi_grades(ImplicitFiltration* fi) {
	ImplicitFiltration& f = *fi; 
	auto G = vector< double >(f.m); 
	f.grading(std::span{G});
	return(Rcpp::wrap(G));
}
NumericVector fi_k_grades(ImplicitFiltration* fi, const size_t k) {
	ImplicitFiltration& f = *fi; 
	auto G= vector< double >(); 
	if (f.n_simplexes.size() > k){
		G.resize(f.n_simplexes.at(k));
		f.k_grades(std::span{G}, k);
	}
	return(Rcpp::wrap(G));
}
IntegerVector simplex_at(ImplicitFiltration* fi, const size_t i) {
	ImplicitFiltration& f = *fi; 
	return(Rcpp::wrap(f.explicit_simplex(i)));
}
double grade_at(ImplicitFiltration* fi, const size_t i) {
	ImplicitFiltration& f = *fi; 
	return(f.explicit_grade(i));
}

size_t dimension(ImplicitFiltration* fi){
	ImplicitFiltration& f = *fi; 
	return(f.d - 1L);
}

vector< size_t > boundary(ImplicitFiltration* fi, const size_t i){
	ImplicitFiltration& f = *fi; 
	vector< size_t > indices; 
	f.boundary(i, true, std::back_inserter(indices));
	return(indices);
}

RCPP_MODULE(implicit_filtration_module) {
  using namespace Rcpp;
  Rcpp::class_< ImplicitFiltration >("ImplicitFiltration")
  	.factory(IF_constructor)
    .method( "as_XPtr", &as_XPtr)
    .field_readonly("size", &ImplicitFiltration::m)
    .field_readonly("num_simplices", &ImplicitFiltration::m)
    .field_readonly("num_vertices", &ImplicitFiltration::n)
    .property("dimension", &dimension, "Returns the dimension.")
    .field_readonly("cum_sizes", &ImplicitFiltration::cum_ns)
    .field_readonly("minimal_grades", &ImplicitFiltration::grades)
    .field_readonly("vertex_labels", &ImplicitFiltration::labels)
    .field_readonly("simplex_ranks", &ImplicitFiltration::ranks)
    .field_readonly("simplex_dims", &ImplicitFiltration::dims)
    .field_readonly("shortlex_perm", &ImplicitFiltration::indexes)
    .property("simplices", &fi_simplices, "Returns the simplices in filtration order as a List.")
	  .property("grading", &fi_grades, "Returns the grading in filtration order as a vector.")
  	.method( "k_simplices", &fi_k_simplices)
  	.method( "k_grades", &fi_k_grades)
  	.method("simplex_at", &simplex_at)
  	.method("grade_at", &grade_at)
  	.method("boundary", &boundary)
    ;
}

// [[Rcpp::export]]
IntegerVector inverse_permutation(const IntegerVector& p){
	return(wrap(inverse_permutation(p.begin(), p.end())));
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