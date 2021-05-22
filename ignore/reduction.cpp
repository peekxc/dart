#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <array>
#include "ref.hpp"
#include "iterator_facade.h"
#include "short_alloc.h"
#include "combinadic.h"
#include "combinations.h"


using std::vector;
using std::array;
using std::unique_ptr;
using std::size_t;

// template < size_t d, typename Field = bool >
// struct Reduction {
// 	using F = Field; 
// 	using boundary_ptr = unique_ptr< vector< size_t > >;
// 	template < typename T > using k_array = array< T, d >;
// 	
// 	size_t n; 
// 	k_array< vector< I > > simplex_ranks;					// each dimensions simplices
// 	k_array< vector< Field > > coeff;							// coefficients associated with each boundary chain / column 
// 	k_array< boundary_ptr > V;     								// upper-triangular portion of V 
// 	k_array< vector< size_t > > row_p;  					// permutation vector to store the order of the rows
// 	
// };


template < typename InputIt, typename Lambda >
void facets(InputIt s, const InputIt e, Lambda&& f){
	const size_t d = std::distance(s, e);
	for_each_combination(s, s+d-1, e, [d, &f](auto fs, auto fe){
		f(fs, fe);
		return false;
	}); 
}

// template < typename InputIt, typename OutputIt >
// void faces(InputIt s, const InputIt e, OutputIt out){
// 	const size_t d = std::distance(s, e);
// 	for (size_t i = 0; i < d; ++i){
// 		for (size_t j = 0; j < d; ++j){
// 			if (i != j){ *out++ = *(s+j); }
// 		}
// 		// std::copy_if(s, e, out, [value](auto label){ return label != value; });
// 	}
// }

// Converts a set of (d-1)-simplex ranks to their faces ranks  
template< size_t d, typename InputIt, typename OutputIt >
void boundary_ranks_(InputIt s, const InputIt e, const size_t n, OutputIt out){
	auto simplex = std::array< size_t, d >();
	// auto face = std::array< size_t, d*(d-1) >();
	for (; s != e; s++){
		lex_unrank(s, s+1, n, d, simplex.begin());	         // labels of current (d-1)-simplex 
		// faces(simplex.begin(), simplex.end(), n, out);
		// faces(simplex.begin(), simplex.end(), face.begin()); // labels of (d-2) faces 
		// lex_rank(face.begin(), face.end(), n, d-1, out);     // ranks of faces
	}
}

template< typename InputIt, typename OutputIt >
void boundary_ranks(InputIt s, const InputIt e, const size_t d, const size_t n, OutputIt out){
	switch(d){
		case 1: 
			return; 
		case 2: 
			boundary_ranks_< 2 >(s, e, n, out);
			break; 
		case 3: 
			boundary_ranks_< 3 >(s, e, n, out);
			break; 
		case 4: 
			boundary_ranks_< 4 >(s, e, n, out);
			break; 
		default: 
			std::invalid_argument("Dimensions supported are [1, 4].");
			break; 
	}
}

// Unranks an ordered set of d-length simplices, applying the boundary operator to each, 
// template < size_t d >
// void d_boundaries(vector< size_t > simplices, ){
// 	
// }

// [[Rcpp::export]]
IntegerMatrix boundaries(IntegerVector simplices, const size_t k, const size_t n){
	auto output_chains = vector< size_t >();
	output_chains.reserve(simplices.size()*(k+1));
	boundary_ranks(simplices.begin(), simplices.end(), k+1, n, std::back_inserter(output_chains));
	return IntegerMatrix(k+1, simplices.size(), output_chains.begin());
}


// Returns either (1) an R = DV (or D = RU) decomposition, or the persistnce pairing defined by R
NumericVector reduce(vector< size_t > simplices, vector< double > fv) {
  
	return(NumericVector::create());
}

// 
template< size_t d, bool compressed, typename InputIt, typename OutputIt >
void facets(InputIt b, const InputIt e, const size_t n, OutputIt out){
	if constexpr(compressed){
		auto c_simplex = array< size_t, d >();
		for(; b != e; ++b){
			lex_unrank(b, b+1, n, d, c_simplex.begin()); // Unrank the current index
			facets(c_simplex.begin(), c_simplex.end(), [n, &c_simplex, &out](auto face_begin, auto face_end){
				lex_rank(face_begin, face_end, n, d-1, out);
			});
		}
	} else {
		for(; b != e; b += d){
			facets(b, b+d, [&out](auto face_begin, auto face_end){
				std::copy(face_begin, face_end, out);
			});
		}
	}
}

template< bool compressed = false, typename InputIt, typename OutputIt >
void facets(InputIt s, const InputIt e, const size_t d, const size_t n, OutputIt out){
	switch(d){
		case 1: 
			return; 
		case 2: 
			facets< 2, compressed >(s, e, n, out);
			break; 
		case 3: 
			facets< 3, compressed >(s, e, n, out);
			break; 
		case 4: 
			facets< 4, compressed >(s, e, n, out);
			break; 
		default: 
			std::invalid_argument("Dimensions supported are [1, 4].");
			break; 
	}
}


// [[Rcpp::export]]
IntegerMatrix facets(IntegerMatrix simplices){
	const size_t d = simplices.nrow();
	
	// Reserve the output
	auto face_labels = vector< size_t >();
	face_labels.reserve((d-1)*simplices.ncol()*d);
	auto output = std::back_inserter(face_labels);
	
	// Generate the facets
	facets(simplices.begin(), simplices.end(), d, simplices.ncol(), output);
	
	// Return as a integer matrix
	IntegerMatrix result = IntegerMatrix(d-1, simplices.ncol()*d, face_labels.begin());
	return(result);
}

// [[Rcpp::export]]
IntegerMatrix facets_ranked(IntegerVector simplices, const size_t d, const size_t n){
	const size_t m = simplices.size();
	
	// Reserve the output
	auto face_labels = vector< size_t >();
	face_labels.reserve(m*d);
	auto output = std::back_inserter(face_labels);
	
	// Generate the facets
	facets< true >(simplices.begin(), simplices.end(), d, n, output);
	
	IntegerMatrix result = IntegerMatrix(d, m, face_labels.begin());
	return(result);
}


/*** R
simplices <- phtools::rank_lex(combn(5,3), n = 5)
phtools:::simplex_faces(combn(5,3))

phtools:::boundaries(simplices = simplices-1L, k = 2L, n = 5L)


phtools:::simplex_faces_compressed(0:9, d = 3, n = 10)


st <- pbgrad::r_geometric_complex(10, radius = 0.25, dim = 2)

e_ids <- phtools::rank_lex(t(st$edges), n = st$n_simplices[1])-1L
t_ids <- phtools:::rank_lex(t(st$triangles), st$n_simplices[1])-1L
facets <- phtools:::facets_ranked(t_ids, d = 3, n = st$n_simplices[1])




*/

	// vector< size_t > simplex(d, 0);
	// for (size_t i = 0; i < m; ++i){
	// 	lex_unrank(simplices.begin()+i, simplices.begin()+i+1, n, d, simplex.begin());
	// 	faces(simplex.begin(), simplex.end(), [n, d, &output](auto face_begin, auto face_end){
	// 		lex_rank(face_begin, face_end, n, d-1, output);
	// 	});
	// }
	// for (size_t i = 0; i < simplices.ncol(); ++i){
	// 	auto simplex = simplices.column(i);
	// 	faces(simplex.begin(), simplex.end(), [&output](auto face_begin, auto face_end){
	// 		std::copy(face_begin, face_end, output);
	// 	});
	// }