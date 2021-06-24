#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"
#include "combinadic.h"

// simplex_birth_time <- function(simplex){
// 	if (length(simplex) == 1){ return(0.0) }
// 	r <- lambda_inorder[simplex]*((1+epsilon)^2)/epsilon
// 	ifelse(all(alpha <= r), min(r), Inf)
// }

// [[Rcpp::export]]
NumericVector sparse_complex(SEXP st_ptr, SEXP st_out, const NumericVector& lambda,	const double alpha, const double epsilon) {
	SimplexTree& st = *XPtr< SimplexTree >(st_ptr);
	SimplexTree& ss = *XPtr< SimplexTree >(st_out);
	const double eps = std::pow(1+epsilon, 2)/epsilon;
	auto edges = st::k_simplices< true >(&st, st.root.get(), 1);
	NumericVector weights; 
  traverse(edges, [&lambda, &weights, &ss, eps, alpha](node_ptr cn, idx_t depth, simplex_t sigma){
    double eps1 = lambda.at(sigma[0]-1)*eps;
  	double eps2 = lambda.at(sigma[1]-1)*eps;
  	if (alpha <= eps1 && alpha <= eps2){
  		ss.insert(sigma);
  		weights.push_back(std::min(eps1, eps2));
  	}
    return true; 
  });
  return(weights);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
