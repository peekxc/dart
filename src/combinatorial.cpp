#include <Rcpp.h>
using namespace Rcpp;

#include <set>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
using std::set; 
using std::string;
using std::vector;

#include "combinadic.h"
#include "combinations.h"

using int_v = vector< int >;

// Returns set containing all LCS for X[0..m-1], Y[0..n-1]
set< int_v > generate_lcs(int_v X, int_v Y, int m, int n, vector< vector< int > >& L) { 
    
  // If we reaches end of either string, return  a empty set 
  set< int_v > s; 
  if (m == 0 || n == 0){ return set< int_v >{ int_v() }; } 

  // If the last characters of X and Y are same 
  // Append current character to all possible LCS of substring X[0..m-2] and Y[0..n-2]. 
  // recurse for X[0..m-2] and Y[0..n-2] in the matrix 
  if (X[m - 1] == Y[n - 1]) { 
    set< int_v > tmp = generate_lcs(X, Y, m - 1, n - 1, L); 
    for (int_v str : tmp){ 
    	str.push_back(X[m - 1]);
    	s.insert(str); 
    }
    // Rprintf("1: s size: %d, temp size: %d\n", s.size(), tmp.size());
  } else { 
    // If LCS can be constructed from top side of the matrix, recurse for X[0..m-2] and Y[0..n-1] 
    if (L[m - 1][n] >= L[m][n - 1]){ s = generate_lcs(X, Y, m - 1, n, L); }
      
    // If LCS can be constructed from left side of the matrix, recurse for X[0..m-1] and Y[0..n-2] 
    // merge two sets if L[m-1][n] == L[m][n-1]
    if (L[m][n - 1] >= L[m - 1][n]) { 
      set< int_v > tmp = generate_lcs(X, Y, m, n - 1, L); 
    	for (int_v v: tmp){ s.insert(v); }
    } 
  } 
  return s; 
} 
  
/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
int lcs_costs(int_v X, int_v Y, int m, int n, vector< vector< int > >& L) { 
  // Build L[m+1][n+1] in bottom up fashion 
  for (int i = 0; i <= m; i++) { 
    for (int j = 0; j <= n; j++) { 
    	L[i][j] = (i == 0 || j == 0) ? 0 : ((X[i - 1] == Y[j - 1]) ? L[i - 1][j - 1] + 1 : std::max(L[i - 1][j], L[i][j - 1]));
    } 
  } 
  return L[m][n]; 
} 

// [[Rcpp::export]]
ListOf< IntegerVector > all_lcs(IntegerVector a, IntegerVector b) {
	vector< vector< int > > L(a.size()+1, vector<int>(b.size()+1, 0));
	auto av = int_v(a.begin(), a.end());
	auto bv = int_v(b.begin(), b.end());
	auto lcs_sz = lcs_costs(av, bv, av.size(), bv.size(), L); // populates L
	set< int_v > all_lcs = generate_lcs(av, bv, av.size(), bv.size(), L);
  return wrap(all_lcs);
}

// [[Rcpp::export]]
List LIS(const IntegerVector& a) {
	const int inf = std::numeric_limits< int >::max(); // note int has no natural infinity()
  int n = a.size();
  vector< int > L(n+1, inf); // L[i] := index of the elements in d[i] 
  vector< int > d(n+1, inf); // d[i] := element in a at which subsequence of length i terminates (T) 
  vector< int > p(n, -1);    // p[i] := index of second last element in LIS ending in i.
  d[0] = -inf;               

  // Build length and ancestor arrays
  for (int i = 0; i < n; i++) {
    int j = std::distance(d.begin(), std::upper_bound(d.begin(), d.end(), a[i]));
    if (d[j-1] < a[i] && a[i] < d[j]){
    	d[j] = a[i];  // record element 
    	p[i] = j; 		// record ancestry 
    }
  }
	return(List::create(_["p"]=p, _["d"]=d));
 
  // Recover the sequence using ancestry array 
  auto last_it = std::max_element(d.begin(), d.end(), [inf](auto di, auto dj){
  	if (abs(di) == inf){ return(true); }
  	if (abs(dj) == inf){ return(false); }
  	return(di < dj);
  });
  // auto last_it = std::find_if(d.begin(), d.end(), [inf](auto di){ return(di < inf); });
  auto pos = std::distance(d.begin(), last_it);
  vector< int > subsequence;
  while (pos != -1) {
  	Rcpp::checkUserInterrupt();
    subsequence.push_back(a[pos]);
    pos = p[pos];
  }
  reverse(subsequence.begin(), subsequence.end());
  return(wrap(subsequence));
  // return(IntegerVector(10));
}

// Returns the number of LIS's in the vector x. Runs in O(n^2).
// [[Rcpp::export]]
int countNumLIS(std::vector<int> x) {
  if (x.size() == 0){ return 0; } // Base Case
  int n = x.size();
  vector< int > dp_l(n, 1);
  vector< int > dp_c(n, 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      if (x[i] <= x[j]){ continue; }
      if (dp_l[j] + 1 > dp_l[i]){
        dp_l[i] = dp_l[j] + 1;
        dp_c[i] = dp_c[j];
      }
      else if (dp_l[j] + 1 == dp_l[i]){
      	dp_c[i] += dp_c[j];
      }
    }
  }
  int count = 0, max_length = *std::max_element(dp_l.begin(), dp_l.end());
  for(int i = 0; i < n; i++){
    if (dp_l[i] == max_length){
      count += dp_c[i];
    }
  }
  return count;
}


// [[Rcpp::export]]
IntegerMatrix inversions(const IntegerVector& a){
	vector< int > A(a.begin(), a.end()); 
	auto inv = vector< int >();
	const size_t n = A.size(); 
	size_t i = 1, j = 0; 
	while(i < n){
		j = i; 
		while(j > 0 && A[j-1] > A[j]){
			inv.push_back(A[j]);
			inv.push_back(A[j-1]);
			std::swap(A[j], A[j-1]);
			j -= 1; 
		}
		i += 1;
	}
	IntegerMatrix out(2, inv.size()/2, inv.begin());
	return(out);
}

[[nodiscard]]
size_t count_inversions(std::span< int > a) noexcept {
  std::multiset< int > set { a[0] };
  size_t count = 0;
  for (size_t i = 1; i < a.size(); ++i) {
    set.insert(a[i]);
    count += std::distance(set.upper_bound(a[i]), set.end());
  }
  return(count);
}

// [[Rcpp::export]]
int inversion_count(IntegerVector iv){
  std::vector< int > v(iv.begin(), iv.end());
  return count_inversions(v);
}

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

// [[Rcpp::export]]
NumericVector perm_dist_mat(IntegerMatrix P, const bool kendall = true, const bool normalize=false){
	
	// Prepare variables
	const size_t n = P.nrow();
	const size_t m = P.ncol(); 
	auto dist_p = NumericVector(dart::BinomialCoefficient(m, 2), 0.0); 
	
	// Indices to track pairwise index comparisons
	auto indices = std::vector< size_t >(m);
	std::iota(indices.begin(), indices.end(), 0);
	
	// Compute all the pairwise distances
	if (kendall){
		// Kendall distance computation
		const size_t max_dist = dart::BinomialCoefficient(n, 2);
		for_each_combination(indices.begin(), indices.begin() + 2, indices.end(), [&](auto b, auto e){
			const size_t i = *b, j = *(b+1);
			IntegerVector matched = Rcpp::match(P.column(j), P.column(i));
			dist_p[dart::lex_rank_2(i, j, m)] = normalize ? double(inversion_count(matched))/double(max_dist) : inversion_count(matched);
			return false; 
		});
	} else {
		// Spearman distance computation
		const IntegerVector offset = Rcpp::seq(1, n); 
		const double max_dist = 0.50*std::pow(n, 2);
		for_each_combination(indices.begin(), indices.begin() + 2, indices.end(), [&](auto b, auto e){
			const size_t i = *b, j = *(b+1);
			size_t sd = Rcpp::sum(Rcpp::abs(offset - Rcpp::match(P.column(i), P.column(j)))); // sugar ftw
			dist_p[dart::lex_rank_2(i, j, m)] = normalize ? sd/max_dist : sd;
			return false; 
		});
	}
	
	// Make it a valid 'dist' object
	dist_p.attr("Size") = m;
	dist_p.attr("Upper") = false; 
	dist_p.attr("class") = "dist";
	dist_p.attr("method") = kendall ? "Kendall" : "Spearman";
	return(dist_p);
}

// [[Rcpp::export]]
size_t fast_choose(const size_t n, const size_t k){
	return(dart::BinomialCoefficient(n,k));
}

// Computes the Spearman distance between two permutations. 
// Precondition: x and y are permutations of the same set. 
// [[Rcpp::export]]
size_t spearman_perm(const IntegerVector& x, const IntegerVector& y){
	const size_t x_len = x.size();
	size_t sd = 0; 
	vector< size_t > matched(x_len);
	for (size_t i = 0; i < x_len; ++i){ matched[y[i]] = i; }
	for (size_t i = 0; i < x_len; ++i){ sd += std::abs(int(i - matched[x[i]])); }
	return(sd);
}

// [[Rcpp::export]]
IntegerVector rank_combn_rcpp(const IntegerMatrix& C, const size_t n){
	const size_t k = C.nrow(); 
	auto out = IntegerVector(C.ncol());
	dart::lex_rank(C.begin(), C.end(), n, k, out.begin());
	return(out);
}
// [[Rcpp::export]]
int rank_combn_single_rcpp(const IntegerVector& C, const size_t n){
	auto rank = dart::lex_rank(std::span(C.begin(), C.end()), n);
	return(rank);
}
// [[Rcpp::export]]
IntegerMatrix unrank_combn_rcpp(const IntegerVector& R, const size_t n, const size_t k){
	auto out_v = vector< int >();
	out_v.reserve(R.size()*k);
	dart::lex_unrank(R.begin(), R.end(), n, k, std::back_inserter(out_v));
	auto out = IntegerMatrix(k, R.size(), out_v.begin());
	return(out);
}
// [[Rcpp::export]]
IntegerVector unrank_combn_single_rcpp(const int64_t rank, const size_t n, const size_t k){
	auto combn = dart::lex_unrank(rank, n, k);
	return(IntegerVector(combn.begin(), combn.end()));
}

// https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
template<typename ... Args>
std::string string_format( const std::string& format, Args ... args ) {
  int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  auto size = static_cast<size_t>( size_s );
  auto buf = std::make_unique<char[]>( size );
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

template< typename RandomIter > [[nodiscard]]
inline std::string pretty_simplex(RandomIter b, const RandomIter e){
	switch(std::distance(b,e)){
		case 0: return("");
		case 1: return(string_format("(%d)", *b ));
		case 2: return(string_format("(%d,%d)", *b, *(b+1) ));
		case 3: return(string_format("(%d,%d,%d)", *b, *(b+1), *(b+2) ));
		case 4: return(string_format("(%d,%d,%d,%d)", *b, *(b+1), *(b+2), *(b+3) ));
		case 5: return(string_format("(%d,%d,%d,%d,%d)", *b, *(b+1), *(b+2), *(b+3), *(b+4) ));
		default: {
			std::stringstream ss;
			ss << "(";
			for (; b != e; ++b){ ss << *b << ","; }
			ss.seekp(-1,ss.cur);
			ss << ")";
			return(ss.str());
		}
	}
}

// [[Rcpp::export]]
StringVector simplex_to_str_rcpp(const ListOf< const IntegerVector >& x){
	std::vector< std::string > out; 
	out.reserve(x.size());
	for (const IntegerVector& s: x){
		out.push_back(pretty_simplex(s.begin(), s.end()));
	}
	return(wrap(out));
}


/*** R
a <- paste0(letters[1:7], collapse = "")
b <- paste0(sample(strsplit(a, "")[[1]]), collapse = "")
phtools:::all_lcs(a, b)


a <- 1:7
b <- sample(a)
phtools:::all_lcs(a, b)

a_sep <- strsplit(a, "")[[1]]
b_sep <- strsplit(b, "")[[1]]
ab_lis <- b_sep[phtools:::longest_subseq(match(b_sep, a_sep))]



*/
