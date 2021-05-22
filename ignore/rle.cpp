#include <Rcpp.h>
using namespace Rcpp;

#include <type_traits>

template<typename T, typename = void>
struct is_equality_comparable : std::false_type
{ };

template<typename T>
struct is_equality_comparable<T,
    typename std::enable_if<
        true, 
        decltype(std::declval<T&>() == std::declval<T&>(), (void)0)
        >::type
    > : std::true_type
{
};

// template< typename T >
// struct RLE {
// 	static_assert(is_equality_comparable< T >::value, "not equality comparable");
// 	
// 	const size_t n;
// 	std::vector< size_t > lens;  
// 	std::vector< T > vals; 
// 	
// 	template< typename Iter >
// 	RLE(Iter b, const Iter e) : n(std::distance(b, e)){
//     for (size_t i = 0; i < n; ++i) {
//       size_t count = 1;
//       while (i < n - 1 && *b == *(b+1)) {
//       	++count;
//         ++i;
//       }
//       vals.push_back(*b);
//       lens.push_back(count);
//     }
// 	}
// 	
// 	// Maps uncompressed range (0, n-1) -> (val) in the compressed range in O(log c) time where c < n
// 	void operator()(size_t i) const {
// 		if (i >= n){ throw std::invalid_argument("i out of range"); }
// 		
// 	}
// 	
// 	
// 	
// 	std::vector< T > inverse() const {
// 		
// 	}
// };


// void printRLE(string str)
// {
//     
// }

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
