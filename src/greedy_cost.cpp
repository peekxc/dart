#include <Rcpp.h>
using namespace Rcpp;


void greedy_spearman(IntegerVector displacement){
	size_t sd = sum(abs(displacement));
	
}

// [[Rcpp::export]]
int interval_cost_rcpp(IntegerVector s, IntegerVector O) {
  const int s1 = s[0];
	const int s2 = s[1];
	int sum = 0;
	if (s1 == s2){ return(0); }
	if (s1 < s2){
		for (auto o_it = O.begin(); o_it != O.end(); o_it += 2){
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
		for (auto o_it = O.begin(); o_it != O.end(); o_it += 2){
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


/*** R

*/
