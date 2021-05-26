#include <Rcpp.h>
using namespace Rcpp;

// void greedy_spearman(IntegerVector displacement){
// 	size_t sd = sum(abs(displacement));
// 
// }

template< typename Iter >
int interval_cost(Iter s, Iter ob, Iter oe) {
  const int s1 = s[0];
	const int s2 = s[1];
	int sum = 0;
	if (s1 == s2){ return(0); }
	if (s1 < s2){
		for (auto o_it = ob; o_it != oe; o_it += 2){
			const int o1 = *o_it;
			const int o2 = *(o_it+1);
			if (o1 == o2){ continue; }
			else if (o1 < o2){
				if (s1 < o1){
					sum += ((s2 < o1 || s2 >= o2) ? 0 : 1);
				} else if (s1 > o1 && s1 < o2) { 
					sum += s2 < o2 ? 0 : -1;
				} 
			} else {
				if (s2 == o2){ sum += 1; } 
				else if (s2 > o2){ sum += s2 > o1 ? 0 : (s2 == o1 ? -1 : 1); } 
			}
		}
	} else {// s1 > s2 
		for (auto o_it = ob; o_it != oe; o_it += 2){
			const int o1 = *o_it;
			const int o2 = *(o_it+1);
			if (o1 == o2){ continue; }
			else if (o1 < o2){
				if (s1 > o1){ sum += s1 <= o2 ? -1 : 0; } 
			} else { 
				if (s1 > o2){ sum += s1 < o1 ? -1 : 0; }
			}
		}
	}
	return(sum);
}

// Given interval [s1, s2] and a set of other intervals O = { [o1,o2], [o3,o4], ..., [ok,ol] }, computes 
// net change in length moving s1 to position s2 would have on all intervals in O
// [[Rcpp::export]]
int interval_cost_rcpp(IntegerVector s, IntegerVector O) {
	return(interval_cost(s.begin(), O.begin(), O.end()));
//   const int s1 = s[0];
// 	const int s2 = s[1];
// 	int sum = 0;
// 	if (s1 == s2){ return(0); }
// 	if (s1 < s2){
// 		for (auto o_it = O.begin(); o_it != O.end(); o_it += 2){
// 			const int o1 = *o_it;
// 			const int o2 = *(o_it+1);
// 			if (o1 == o2){ continue; }
// 			else if (o1 < o2){
// 				if (s1 < o1){
// 					sum += ((s2 < o1 || s2 >= o2) ? 0 : 1);
// 				} else if (s1 > o1 && s1 < o2) { 
// 					sum += s2 < o2 ? 0 : -1;
// 				} 
// 			} else {
// 				if (s2 == o2){ sum += 1; } 
// 				else if (s2 > o2){ sum += s2 > o1 ? 0 : (s2 == o1 ? -1 : 1); } 
// 			}
// 		}
// 	} else {// s1 > s2 
// 		for (auto o_it = O.begin(); o_it != O.end(); o_it += 2){
// 			const int o1 = *o_it;
// 			const int o2 = *(o_it+1);
// 			if (o1 == o2){ continue; }
// 			else if (o1 < o2){
// 				if (s1 > o1){ sum += s1 <= o2 ? -1 : 0; } 
// 			} else { 
// 				if (s1 > o2){ sum += s1 < o1 ? -1 : 0; }
// 			}
// 		}
// 	}
// 	return(sum);
}

// sapply(1:ncol(moves), function(i){ dart:::interval_cost_rcpp(moves[,i], moves[,-i]) })
// [[Rcpp::export]]
IntegerVector pairwise_cost(const IntegerMatrix& M){
	IntegerVector indices(M.begin(), M.end());
	const size_t n = indices.size()/2; 
	IntegerVector res(n, 0);
	for (size_t k = 0; k < n; ++k){
		res[k] = interval_cost(indices.begin(), indices.begin() + 2, indices.end());
		if (k < (n - 1)){
			std::swap(indices[0], indices[2*(k+1)]);
			std::swap(indices[1], indices[2*(k+1)+1]);
		}
	}
	return(res);
}


/*** R

*/
