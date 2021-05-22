#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector unique_numeric(std::vector< double > x, const double eps){
	auto it = std::unique(x.begin(), x.end(), [eps](double x1, double x2){
		return(std::abs(x1 - x2) < eps);
	}); 
	x.resize(std::distance(x.begin(), it));
	return(wrap(x));
}

typedef int index_type;

// From: https://www.r-bloggers.com/2014/09/compute-longest-increasingdecreasing-subsequence-using-rcpp/
// [[Rcpp::export]]
Rcpp::IntegerVector longest_inc_subseq(SEXP X){
	Rcpp::IntegerVector x(X);
	std::vector< int > P(x.size(), 0);
	std::vector< int > M(x.size()+1, 0);
	index_type L(0), newL;
	for(index_type i=0; i < x.size(); ++i) {
		index_type lo(1), hi(L), mid;
		while( lo <= hi) { mid = (lo + hi) / 2;
			if (x[M[mid]] < x[i] ) { lo = mid + 1; } else { hi = mid - 1; } 
		}
	  newL = lo;
	  P[i] = M[newL - 1];
	  if (newL > L) {
	  	M[newL] = i;
	  	L = newL;
	  } else if (x[i] < x[M[newL]]) {
	  	M[newL] = i;
	  }
	}    
	std::vector< index_type > re(L);
	index_type k(M[L]);
	for(index_type i=L-1; i>=0; --i){
		re[i] = k + 1;
		k = P[k];
	}    
	return(Rcpp::wrap(re));
}
