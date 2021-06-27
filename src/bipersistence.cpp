#include <Rcpp.h>
using namespace Rcpp;


// Given a set of points (column-oriented) in the plane in a line parameterized by y = m*x + b, 
// applies the 'push_map' projecting each point onto its least upper bound (LUB) on the line
// [[Rcpp::export]]
NumericMatrix push_map(const NumericMatrix& x, const double m, const double b) {
	const size_t n_pts = x.ncol(); 
	auto xb = x.begin(); 
	auto out = NumericMatrix(n_pts, 3);
	size_t i = 0; 
	for (double xc, yc; xb != x.end(); xb += 2){
		xc = *xb; yc = *(xb+1);
		NumericVector tmp = (xc*m + b < yc) ? NumericVector::create((yc-b)/m, yc) : NumericVector::create(xc, m*xc + b);
		double dist_to_end = std::sqrt(std::pow(tmp[0], 2) + std::pow(b - tmp[1], 2));
		tmp.push_back(dist_to_end);
		out(i++, _) = tmp;
	}
	return(out);
}


// NumericMatrix dist_along_line(const NumericMatrix& x, const double m, const double b) {
// 	const size_t n_pts = x.ncol(); 
// 	auto xb = x.begin(); 
// 	auto out = NumericMatrix(n_pts, 2);
// 	size_t i = 0; 
// 	for (double xc, yc; xb != x.end(); xb += 2){
// 		xc = *xb; yc = *(xb+1);
// 		out(i++, _) = (xc*m + b < yc) ? NumericVector::create((yc-b)/m, yc) : NumericVector::create(xc, m*xc + b);
// 	}
// 	return(out);
// }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
