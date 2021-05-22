#include <Rcpp.h>
using namespace Rcpp;

#include "reduction_concepts.h"
// #include "PspMatrix.h"




// // Square transposition framework
// template< typename Lambda >
// void transpose_schedule_full(PspBoolMatrix& R, PspBoolMatrix& V, const vector< size_t >& S, Lambda f){
// 	const size_t nc = R.n_cols();
// 	const size_t nr = R.n_rows();
// 	if (nc == 0 || nr == 0){ return; }
// 	if (nc != nr){ throw std::invalid_argument("R must be square."); }
// 	auto max_el = std::max_element(S.begin(), S.end());
// 	if (*max_el >= (nc-1)){ throw std::invalid_argument("Given indices exceed matrix dimensions"); }
// 	
// 	// // Returns the current row index of the lowest entry in column j
// 	// const auto low_index = [&R](size_t j) -> optional< size_t >{
// 	// 	if (R.column_empty(j)){ return(std::nullopt); }
// 	// 	auto& c = *R.columns[j];
// 	// 	auto me = std::max_element(c.begin(), c.end(), [&R](auto& e1, auto& e2){
// 	// 		return R.otc[e1.first] < R.otc[e2.first];
// 	// 	});
// 	// 	return(R.otc[(*me).first]);
// 	// };
// 	
// 	// Construct three O(n)-sized maps to make various checks in the loop O(1)
// 	auto positive = vector< bool >(nc); 					// whether a column is a creator or destroyer 
// 	auto Low = vector< optional< size_t > >(nc);	// map from column index --> low index 
// 	auto RLow = vector< optional< size_t > >(nc);	// map from low index --> column index
// 	for (size_t j = 0; j < nc; ++j){
// 		positive[j] = R.column_empty(j);
// 		Low[j] = R.lowest_nonzero(j);
// 		if (Low[j]){
// 			RLow[Low[j].value()] = make_optional(j);
// 		}
// 	}
// 	
// 	// Perform the transpositions
// 	size_t status = 0; 
// 	size_t line = 1; 
// 	for (auto i: S){
// 		auto j = i + 1;
// 		Rcout << line << ": swapping(0) <=> " << i << "," << j << std::endl; 
// 		if (positive[i] && positive[j]){
// 			if (V(i,j) != 0){ V.add_cols(j, i); } 
// 			if (RLow[i] && RLow[j] && (R(i, *RLow[j]) != 0)){ 
// 				size_t k = RLow[i].value(), l = RLow[j].value();
// 				// Rprintf("k: %d, l: %d \n", k, l);
// 				status = k < l ? 1 : 2; // Cases 1.1.1 and 1.1.2
// 				if (l < k){ std::swap(k,l); }
// 				// R.swap(i,j); R.add_cols(l,k);
// 				// V.swap(i,j); V.add_cols(l,k);
// 				R.add_cols(l,k); R.swap(i,j);
// 				V.add_cols(l,k); V.swap(i,j); 
// 			} else {
// 				status = 3; // Case 1.2
// 			}
// 		} else if (!positive[i] && !positive[j]){
// 			if (V(i,j) != 0){ // Case 2.1
// 				if (!bool(Low[i]) || !bool(Low[j])){ Rcpp::stop("Positive not maintained."); }
// 				if (Low[i].value() < Low[j].value()){ 
// 					R.add_cols(j,i); R.swap(i,j);
// 					V.add_cols(j,i); V.swap(i,j); 
// 					status = 4; // Case 2.1.1
// 				} else {
// 					R.add_cols(j,i); R.swap(i,j); R.add_cols(j,i);
// 					V.add_cols(j,i); V.swap(i,j); V.add_cols(j,i);
// 					status = 5; // Case 2.1.2
// 				}
// 			} else { 
// 				status = 6; // Case 2.2
// 			}
// 		} else if (!positive[i] && positive[j]){
// 			if (V(i,j) != 0){
// 				R.add_cols(j,i); R.swap(i,j); R.add_cols(j,i);
// 				V.add_cols(j,i); R.swap(i,j); V.add_cols(j,i);
// 				status = 7; // Case 3.1
// 			} else {
// 				status = 8; // Case 3.2
// 			}
// 		} else {
// 			status = 9; // Case 4
// 			if (V(i,j) != 0){ V.add_cols(j,i); }
// 		}
// 		// Cases 1.2, 2.2, 3.2, and 4 just need a single permutation
// 		if (status == 0 || status == 9 || status == 8 || status == 6 || status == 3){
// 			// Rcout << "swapping: " << i << ", " << j << std::endl;
// 			R.swap(i,j); V.swap(i,j);
// 		}
// 		
// 		// Update all the maps (all in O(1))
// 		// positive[i] = R.column_empty(i);
// 		// positive[j] = R.column_empty(j);
// 		
// 		Rcout << " --- Updating maps --- " << std::endl;
// 		// Columns to update low entries of 
// 		// auto ck = RLow[i];
// 		// auto cl = RLow[j];
// 		// 
// 		// Rcout << " 1" << std::endl;
// 		// Low[i] = R.lowest_nonzero(i);
// 		// Low[j] = R.lowest_nonzero(j);
// 		// if (ck){ Low[*ck] = R.lowest_nonzero(*ck); }
// 		// if (cl){ Low[*cl] = R.lowest_nonzero(*cl); }
// 		// 
// 		// Rcout << " 2" << std::endl;
// 		// positive[i] = !bool(Low[i]);
// 		// positive[j] = !bool(Low[j]);
// 		// if (ck){ positive[*ck] = !bool(Low[*ck]); }
// 		// if (cl){ positive[*cl] = !bool(Low[*cl]); }
// 		
// 		// Rcout << " 3" << std::endl;
// 		// if (Low[i]){ RLow[Low[i].value()] = make_optional(i); }
// 		// if (Low[j]){ RLow[Low[j].value()] = make_optional(j); }
// 		// if (ck){ RLow[Low[*ck].value()] = make_optional(*ck); }
// 		// if (cl){ RLow[Low[*cl].value()] = make_optional(*cl); }
// 		
// 		
// 		for (size_t c = 0; c < nc; ++c){
// 			positive[c] = R.column_empty(c);
// 		}
// 		for (size_t c = 0; c < nc; ++c){
// 			Low[c] = R.lowest_nonzero(c);
// 		}
// 		RLow = vector< optional< size_t > >(nc);
// 		for (size_t c = 0; c < nc; ++c){
// 			if (Low[c]){
// 				RLow[Low[c].value()] = make_optional(c);
// 			}
// 		}
// 		R.clean(false);
// 		V.clean(false);
// 		// Low[i] = R.lowest_nonzero(i); 
// 		// Low[j] = R.lowest_nonzero(j); 
// 		// for (size_t c = 0; c < nc; ++c){
// 		// 	auto row_entry = std::find_if(Low.begin(), Low.end(), [c](optional< size_t >& le){
// 		// 		return(bool(le) && (*le == c));
// 		// 	});
// 		// 	RLow[c] = (row_entry != Low.end()) ? make_optional(std::distance(Low.begin(), row_entry)) : std::nullopt; 
// 		// }
// 		// if (low_i){ 
// 		// 	auto row_entry = std::find_if(Low.begin(), Low.end(), [&low_i](optional< size_t >& le){
// 		// 		return(bool(le) && (*le == *low_i));
// 		// 	});
// 		// 	RLow[Low[i].value()] = make_optional(std::distance(Low.begin(), row_entry)); 
// 		// }
// 		// if (low_j){ 
// 		// 	auto row_entry = std::find_if(Low.begin(), Low.end(), [&low_j](optional< size_t >& le){
// 		// 		return(bool(le) && (*le == *low_j));
// 		// 	});
// 		// 	RLow[Low[j].value()] = make_optional(std::distance(Low.begin(), row_entry)); 
// 		// }
// 		
// 		// if (!positive[i]){
// 		// 	const auto& col_i = *R.columns[i];
// 		// 	auto me1 = std::max_element(col_i.begin(), col_i.end(), [&R](auto& e1, auto& e2){
// 		// 		return(R.otc[e1.first] < R.otc[e2.first]);
// 		// 	});
// 		// 	Low[i] = R.otc[(*me1).first];
// 		// } else {
// 		// 	Low[i] = std::nullopt;
// 		// }
// 		// if (!positive[j]){
// 		// 	const auto& col_j = *R.columns[j];
// 		// 	auto me2 = std::max_element(col_j.begin(), col_j.end(), [&R](auto& e1, auto& e2){
// 		// 		return(R.otc[e1.first] < R.otc[e2.first]);
// 		// 	});
// 		// 	Low[j] = R.otc[(*me2).first];
// 		// } else {
// 		// 	Low[j] = std::nullopt;
// 		// }
// 		// Low[i] = positive[i] ? std::nullopt : make_optional(R.cto[R.columns[i]->back().first]);
// 		// Low[j] = positive[j] ? std::nullopt : make_optional(R.cto[R.columns[j]->back().first]);
// 		
// 		// if (low_i){
// 		// 	if (Low[i] && Low[i].value() == *low_i){
// 		// 		RLow[*low_i] = i;
// 		// 	} else if (Low[j] && Low[j].value() == *low_i){
// 		// 		RLow[*low_i] = j;
// 		// 	} else {
// 		// 		RLow[*low_i] = std::nullopt; 
// 		// 	}
// 		// }
// 		// if (low_j){
// 		// 	if (Low[i] && Low[i].value() == *low_j){
// 		// 		RLow[*low_j] = i;
// 		// 	} else if (Low[j] && Low[j].value() == *low_j){
// 		// 		RLow[*low_j] = j;
// 		// 	} else {
// 		// 		RLow[*low_j] = std::nullopt; 
// 		// 	}
// 		// }
// 		// if (low_j){ RLow[Low[j].value()] = make_optional(j); }
// 		
// 		Rcout << line << ": Positive: ";
// 		for (size_t c = 0; c < nc; ++c){
// 			Rcout << (positive[c] ? 1 : 0) << ",";
// 		}
// 		Rcout << std::endl;
// 		
// 		Rcout << line << ": Low: ";
// 		for (size_t c = 0; c < nc; ++c){
// 			Rcout << (bool(Low[c]) ? std::to_string(*Low[c]) : std::string("NA")) << std::string(",");
// 		}
// 		Rcout << std::endl;
// 		
// 		Rcout << line << ": Low (reported): ";
// 		for (size_t c = 0; c < nc; ++c){
// 			auto le = R.lowest_nonzero(c);
// 			Rcout << (bool(le) ? std::to_string(*le) : std::string("NA")) << std::string(",");
// 		}
// 		Rcout << std::endl;
// 		
// 		Rcout << line << ": RLow: ";
// 		for (size_t c = 0; c < nc; ++c){
// 			Rcout << (bool(RLow[c]) ? std::to_string(*RLow[c]) : std::string("NA")) << std::string(",");
// 		}
// 		Rcout << std::endl;
// 		Rcout << std::endl;
// 		
// 		// Apply user function
// 		f(R, V, status);
// 		line++;
// 	}
// }