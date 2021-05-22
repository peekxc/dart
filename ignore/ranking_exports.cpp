#include <Rcpp.h>
using namespace Rcpp;

#include "combinadic.h"
#include <cmath>

// [[Rcpp::export]]
size_t binomial_coefficient(const size_t n, const size_t k){
	return BinomialCoefficient(n,k);
}

// 0-based conversion of natural number to (n choose k) combinadic
// [[Rcpp::export]]
IntegerMatrix unrank_R(const IntegerVector& ranks, const size_t n, const size_t k){
	auto result = std::vector< size_t >(); 
	auto out = std::back_inserter(result);
	lex_unrank(ranks.begin(), ranks.end(), n, k, out);
	auto results = IntegerMatrix(k, ranks.size(), result.begin());
	return(results);
}

// 0-based conversion of (n choose k) combinadic subscripts to natural numbers. Expects column-matrix
// [[Rcpp::export]]
IntegerVector rank_R(const IntegerMatrix& T, const size_t n){
  const size_t k = T.nrow();
  auto result = std::vector< size_t >(); 
	auto out = std::back_inserter(result);
	lex_rank(T.begin(), T.end(), n, k, out);
  return(wrap(result));
}


// IntegerVector unrank(IntegerVector ranks, const size_t n, const size_t k){
// 	auto result = std::vector< size_t >(); 
// 	auto out = std::back_inserter(result);
// 	for (auto r: ranks){
// 		lex_unrank_k(r, k, n, out);
// 	}
// 	return(wrap(result));
// }

// [[Rcpp::export]]
size_t benchmark_bc(const size_t n){
	size_t result = 0; 
	for (size_t i = 1; i < n; ++i){
		for (size_t j = i+1; j < n; ++j){
			result += binomial_coefficient(j, i);
		}
	}
	return(result);
}

template < bool column_major = true, typename InputIt, typename OutputIt >
void unrank_grid_(InputIt s, const InputIt e, const size_t nr, OutputIt out){
	using value_t = typename std::iterator_traits< InputIt >::value_type; 
	value_t r, c;
	double x; 
	for(; s != e; ++s){ 
		x = double(*s); 
		r = (value_t(x) % nr);
		c = std::floor(x / double(nr));
		if constexpr (!column_major){ std::swap(r,c); }
		*out++ = r; 
		*out++ = c;
		// Rprintf("(x,r,c):(%d,%d,%d)\n",x,r,c);
	}
}

template < typename InputIt, typename OutputIt >
void unrank_grid(InputIt s, const InputIt e, const size_t nr, const size_t nc, bool column_major, OutputIt out){
	if (column_major){
		unrank_grid_< true >(s, e, nr, out);
	} else {
		unrank_grid_< false >(s, e, nc, out);
	}
}
	
template < bool column_major = true, typename InputIt, typename OutputIt >
void rank_grid_(InputIt s, const InputIt e, const size_t nc, OutputIt out){
	if constexpr (column_major){
		for(; s != e; s += 2){ *out++ = (*(s+1))*nc + (*s); }
	} else {
		for(; s != e; s += 2){ *out++ = (*s)*nc + *(s+1); }
	}
}

template <  typename InputIt, typename OutputIt >
void rank_grid(InputIt s, const InputIt e, const size_t nr, const size_t nc, bool column_major, OutputIt out){
	if (column_major){
		rank_grid_< true >(s, e, nr, out);
	} else {
		rank_grid_< false >(s, e, nc, out);
	}
}


// [[Rcpp::export]]
IntegerVector rank_gridR(const IntegerMatrix& T, const size_t nr, const size_t nc, bool column_major=true){
	std::vector< I > output; 
	output.reserve(T.ncol());
	auto ins = std::back_inserter(output);
	rank_grid(T.begin(), T.end(), nr, nc, column_major, ins);
	return(wrap(output));
}

// [[Rcpp::export]]
IntegerMatrix unrank_gridR(const IntegerVector& T, const size_t nr, const size_t nc, bool column_major=true){
	std::vector< I > output; 
	output.reserve(T.size() * 2);
	auto ins = std::back_inserter(output);
	unrank_grid(T.begin(), T.end(), nr, nc, column_major, ins);
	IntegerMatrix result(2, T.size(), output.begin());
	return(wrap(result));
}

/*** R
time_cached <- sapply(5:25, function(i){ 
	b <- microbenchmark::microbenchmark({ phtools:::benchmark_bc(i) }, times = 1L)
	b$time
})

time_uncached <- sapply(5:25, function(i){ 
	b <- microbenchmark::microbenchmark({ phtools:::benchmark_bc(i) }, times = 1L)
	b$time
})

*/
