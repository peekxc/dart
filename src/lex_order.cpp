#include <Rcpp.h>
using namespace Rcpp;

// Returns a vector of indices 'idx' such that simplices[idx] is sorted in increasing order, with faces before 
// cofaces, breaking ties in a lexicographically-increasing manner. 
// [[Rcpp::export]]
IntegerVector order_simplices(List simplices, NumericVector weights) {
	if (weights.size() != simplices.size()){ Rcpp::stop("weights length must match simplices length."); }
	const size_t m = simplices.size(); 
  std::vector< size_t > ids(m);
  std::iota(ids.begin(), ids.end(), 0);
  
  const auto machine_eq = [](double d1, double d2){
  	return std::abs(d1 - d2) < std::numeric_limits< double >::epsilon();
  };
  
	std::sort(ids.begin(), ids.end(), [&machine_eq, &weights, &simplices](size_t i, size_t j){
  	if (machine_eq(weights[i], weights[j])){
  		IntegerVector si = simplices[i];
  		IntegerVector sj = simplices[j];
  		return(si.size() != sj.size() ? si.size() < sj.size() : std::lexicographical_compare(si.begin(), si.end(), sj.begin(), sj.end()));
  	}
  	else if (weights[i] < weights[j]){ return true; }
		else { return false; }
  });
	return(wrap(ids));
}

/*** R
si <- list(3, 2, 1, c(1,2), c(2,3), c(1,2,3))
wi <- c(0, 0, 0, 0.3, 0, 0.50)
si[order_simplices(si, wi)+1]

si <- list(1,2,3,c(1,2),c(2,3),c(1,3),c(1,2,3))
si[order_simplices(si,weights = rep(0.0, 7))+1]
*/
