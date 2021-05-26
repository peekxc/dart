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
#include "combinadic.h"

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


#include "Rcpp.h"
using namespace Rcpp;

using std::vector;

// Restore right 
template< PermutableMatrix Matrix, typename Iter >
void restore_right(Matrix& R, Matrix& V, Iter bi, const Iter ei, vector< typename Matrix::entry_t >& dr, vector< typename Matrix::entry_t >& dv){
	using entry_t = typename Matrix::entry_t;

	// Start with donor columns
	size_t donor_idx = *bi;
	auto d_low_index = R.low_index(donor_idx);
	dr.clear(); dv.clear();
	R.write_column(donor_idx, std::back_inserter(dr));
	V.write_column(donor_idx, std::back_inserter(dv));
	if (std::distance(bi, ei) <= 1){ return; }

	// Apply the donor concept
	std:advance(bi, 1);
	auto new_dr = vector< entry_t >(), new_dv = vector< entry_t >();
	for (size_t k = *bi; bi != ei; ++bi){
		auto new_low_index = R.low_index(k);

		// Copy column into new donor column
		new_dr.clear(); new_dv.clear();
		R.write_column(donor_idx, std::back_inserter(new_dr));
		V.write_column(donor_idx, std::back_inserter(new_dv));

		// Add the columns
		R.add_cols(k, donor_idx); V.add_cols(k, donor_idx);

		// Replace the donor columns
		if ((d_low_index && new_low_index) && *d_low_index > *new_low_index){
			d_low_index = std::move(new_low_index);
			dr.swap(new_dr);
			dv.swap(new_dv);
		}
	}
}

// Performs a sequence of moves (i,j) on the pair (R1, V1), correcting (R2, V2) as necessary. 
template< PermutableMatrix Matrix, typename Iter, typename Lambda >
void move_schedule_local(Matrix& R1, Matrix& V1, Matrix& R2, Matrix& V2, Iter sb, const Iter se, Lambda f){
	using entry_t = typename Matrix::entry_t;
	using value_t = typename Matrix::value_type;
	const size_t nr1 = R1.dim().first;
	const size_t nc1 = R1.dim().second;
	const size_t nr2 = R2.dim().first;
	const size_t nc2 = R2.dim().second;
	if (nc1 != nr2){ throw std::invalid_argument("Number of rows in R2 must match number of columns in R1."); }
	if (nc1 == 0 || nr1 == 0){ return; }
	if (std::distance(sb,se) < 2 || std::distance(sb,se) % 2 != 0){ throw std::invalid_argument("Pairs of indices must be passed."); }

	for (size_t i, j; sb != se; sb += 2){
		i = *sb, j = *(sb+1);
		if (i == j){ continue; }
		if (i > j || i >= (nc1-1) || j >= nc1){  throw std::invalid_argument("Invalid pairs (i,j) passed."); }
		
		// Collect indices I
		auto I = vector< size_t >();
		V1.row(i, [&](auto col_idx, auto v){
			if ((col_idx >= i) && (col_idx <= j)){ I.push_back(col_idx); }
		});

		// Collect indices J
		auto J = vector< size_t >();
		for (size_t c = 0; c < nc2; ++c){
			auto low_idx = R2.low_index(c);
			if (low_idx && (*low_idx) >= i && (low_idx) <= j && R2(i,c) != 0){
				J.push_back(c);
			}
		}
		
		// Restore invariants to columns affected by I
		auto dr1 = vector< entry_t >();
		auto dv1 = vector< entry_t >();
		restore_right(R1, V1, I.begin(), I.end(), dr1, dv1); // don't wrap in if condition

		// Restore invariants to columns/rows affected by J
		if (!J.empty()){
			auto dr2 = vector< entry_t >();
			auto dv2 = vector< entry_t >();
			restore_right(R2, V2, J.begin(), J.end(), dr2, dv2);
		}

		// Apply permutations
		vector< size_t > p(nc1);
		std::iota(p.begin(), p.end(), 0);
		std::rotate(p.begin()+i, p.begin()+i+1, p.begin()+j);

		R1.permute_cols(p.begin(), p.end());
		V1.permute(p.begin(), p.end());
		R2.permute_rows(p.begin(), p.end());

		// Apply permutations and assign donors
		apply_permutation(dr1.begin(), dr1.end(), p.begin());
		apply_permutation(dv1.begin(), dv1.end(), p.begin());
		R1.assign_column(j, dr1.begin(), dr1.end());
		V1.assign_column(j, dv1.begin(), dv1.end());
	}
	
}

// Given a contiguous schedule of integer indices (i,j) via [sb, se), performs a sequence of move operations
// on the iterator range
template< PermutableMatrix Matrix, typename Iter, typename Lambda >
void move_schedule_full(Matrix& R, Matrix& V, Iter sb, const Iter se, Lambda f){
	using entry_t = typename Matrix::entry_t;
	using value_t = typename Matrix::value_type;
	const size_t nc = R.dim().first;
	const size_t nr = R.dim().second;
	if (nc == 0 || nr == 0){ return; }
	if (nc != nr){ throw std::invalid_argument("R must be square."); }
	if (std::distance(sb,se) < 2 || std::distance(sb,se) % 2 != 0){  throw std::invalid_argument("Pairs of indices must be passed."); }

	for (size_t i, j; sb != se; sb += 2){
		i = *sb, j = *(sb+1);
		if (i == j){ continue; }
		if (i > j || i >= (nc-1) || j >= nc){  throw std::invalid_argument("Invalid pairs (i,j) passed."); }

		// Collect indices I
		auto I = vector< size_t >();
		V.row(i, [&](auto col_idx, auto v){
			if ((col_idx >= i) & (col_idx <= j)){ I.push_back(col_idx); }
		});

		// Collect indices J
		auto J = vector< size_t >();
		for (size_t c = 0; c < nc; ++c){
			auto low_idx = R.low_index(c);
			if (low_idx && (*low_idx) >= i && (low_idx) <= j && R(i,c) != 0){
				J.push_back(c);
			}
		}

		Rcout << "I: ";
		for (auto ii: I){ Rcout << ii << ", "; }
		Rcout << std::endl;

		Rcout << "J: ";
		for (auto ii: J){ Rcout << ii << ", "; }
		Rcout << std::endl;


		// Restore invariants
		auto dr1 = vector< entry_t >();
		auto dv1 = vector< entry_t >();
		restore_right(R, V, I.begin(), I.end(), dr1, dv1); // don't wrap in if condition

		// Restore invariants
		if (!J.empty()){
			auto dr2 = vector< entry_t >();
			auto dv2 = vector< entry_t >();
			restore_right(R, V, J.begin(), J.end(), dr2, dv2);
		}

		// Apply permutations
		vector< size_t > p(nc);
		std::iota(p.begin(), p.end(), 0);
		std::rotate(p.begin()+i, p.begin()+i+1, p.begin()+j);

		for (auto pe: p){ Rcout << pe << " "; }
		Rcout << std::endl;

		R.permute(p.begin(), p.end());
		V.permute(p.begin(), p.end());

		Rcout << "Donor dr: ";
		for (auto pe: dr1){ Rcout << pe.second << " "; }
		Rcout << std::endl;

		Rcout << "Donor dv: ";
		for (auto pe: dv1){ Rcout << pe.second << " "; }
		Rcout << std::endl;

		// Perform donor replacement
		apply_permutation(dr1.begin(), dr1.end(), p.begin());
		apply_permutation(dv1.begin(), dv1.end(), p.begin());

		Rcout << "Donor dr: ";
		for (auto pe: dr1){ Rcout << pe.second << " "; }
		Rcout << std::endl;

		Rcout << "Donor dv: ";
		for (auto pe: dv1){ Rcout << pe.second << " "; }
		Rcout << std::endl;

		R.assign_column(j, dr1.begin(), dr1.end());
		V.assign_column(j, dv1.begin(), dv1.end());
	}
}
