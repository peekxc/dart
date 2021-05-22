#include <string>
#include <cstddef>
#include <concepts>
#include <iostream>
#include <optional>
#include <tuple>

using std::pair; 
using std::make_pair;
using std::size_t; 
using std::optional;
using std::tuple;
using std::make_tuple;

#include "reduction_concepts.h"

template< bool clearing = true, ReducibleMatrix Matrix, typename Iter >
void pHcol_local(Matrix& R1, Matrix& V1, Matrix& R2, Matrix& V2, Iter b1, const Iter e1, Iter b2, const Iter e2){
	using field_t = typename Matrix::value_type;
	using entry_t = typename std::pair< size_t, field_t >;

	// Reduce (R2,V2)
	pHcol(R2, V2, b2, e2);

	// Apply clearing optimization
	if constexpr(clearing){
		optional< entry_t > low_j;
		for (size_t j; b2 != e2; ++b2){
			j = *b2; 
			if (!R2.column_empty(j) && (low_j = R2.low(j))){
				R1.clear_column(low_j->first);
				auto R_j = std::vector< entry_t >();
				R2.column(j, [&](auto row_idx, auto v){ R_j.push_back(std::make_pair(row_idx, v)); });
				V1.assign_column(low_j->first, R_j.begin(), R_j.end());
			}
		}
	}
	
	// Now reduce (R1, V1)
	pHcol(R1, V1, b1, e1);
	return;
}


template< ReducibleMatrix Matrix, typename Iter >
void pHcol(Matrix& R, Matrix& V, Iter b, const Iter e){
	using field_t = typename Matrix::value_type;
	using entry_t = typename std::pair< size_t, field_t >;
	
	// Reduction algorithm
	const size_t m = R.dim().first;
	for (size_t j; b != e; ++b){
		j = *b; 
		std::cout << "Reducing " << j << std::endl; 
		optional< entry_t > low_i, low_j;
		while((low_j = R.low(j)) && (low_i = R.find_low(j, low_j->first))){
			// Rprintf("quotienting(%d/%d): %.2g / %.2g\n", low_j->first, low_i->first, low_j->second, low_i->second);
			size_t i = low_i->first;
			auto lambda = low_j->second / low_i->second;
			std::cout << "i = " << i << ", j = " << j << ", lambda = " << lambda << std::endl; 
			std::cout << "pivot_i = " << low_i->second << ", pivot_j = " << low_j->second << std::endl; 
			R.add_scaled_col(i, j, -lambda);
			V.add_scaled_col(i, j, -lambda);
		}
	}
	return; 
}


// Square transposition framework
template< PermutableMatrix Matrix, typename Iter, typename Lambda >
void transpose_schedule_full(Matrix& R, Matrix& V, Iter sb, const Iter se, Lambda f){
	using entry_t = typename Matrix::entry_t;
	const size_t nc = R.dim().first;
	const size_t nr = R.dim().second;
	if (nc == 0 || nr == 0){ return; }
	if (nc != nr){ throw std::invalid_argument("R must be square."); }
	// auto max_el = std::max_element(sb, se);
	// if (*max_el >= (nc-1)){ throw std::invalid_argument("Given indices exceed matrix dimensions"); }

	// Perform the transpositions
	size_t status = 0;
	size_t line = 1;
	optional< entry_t > low_i, low_j;
	for (size_t i = *sb; sb != se; ++sb, i = *sb){
		auto j = i + 1;
		if (R.column_empty(i) && R.column_empty(j)){
			if (V(i,j) != 0){ V.cancel_lowest(i,j); } // cancel lowest of j using lowest entry in i
			if ((low_i = R.low(i)) && (low_j = R.low(j)) && (R(i, low_j->second) != 0)){
				size_t k = low_i->first, l = low_j->second;
				status = k < l ? 1 : 2; 
				if (status == 1){ // Cases 1.1.1
					R.cancel_lowest(k,l); R.swap(i,j);
					V.cancel_lowest(k,l); V.swap(i,j);
				} else {          // Cases 1.1.2
					R.cancel_lowest(l,k); R.swap(i,j);
					V.cancel_lowest(l,k); V.swap(i,j);
				}
			} else {
				status = 3; // Case 1.2
			}
		} else if (!R.column_empty(i) && !R.column_empty(j)){
			if (V(i,j) != 0){ // Case 2.1
				if (R.low_index(i) < R.low_index(j)){
					R.cancel_lowest(i,j); R.swap(i,j);
					V.cancel_lowest(i,j); V.swap(i,j);
					status = 4; // Case 2.1.1
				} else {
					R.cancel_lowest(i,j); R.swap(i,j); R.cancel_lowest(i,j);
					V.cancel_lowest(i,j); V.swap(i,j); V.cancel_lowest(i,j);
					status = 5; // Case 2.1.2
				}
			} else {
				status = 6; // Case 2.2
			}
		} else if (!R.column_empty(i) && R.column_empty(j)){
			if (V(i,j) != 0){
				R.cancel_lowest(i,j); R.swap(i,j); R.cancel_lowest(i,j);
				V.cancel_lowest(i,j); R.swap(i,j); V.cancel_lowest(i,j);
				status = 7; // Case 3.1
			} else {
				status = 8; // Case 3.2
			}
		} else {
			status = 9; // Case 4
			if (V(i,j) != 0){ V.cancel_lowest(i,j); }
		}
		f(status); // Apply user function
		line++;
	}
}