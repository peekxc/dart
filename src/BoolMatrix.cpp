
#include "PspMatrix.h" // Permutable Sparse Matrix template
#include "reduction.h" // reduction templates

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Define the field operation and specialize
struct Mod2 { 
	constexpr bool operator()(bool b1, bool b2) { return(b1 ^ b2); }
};
typedef PspMatrix< bool, Mod2 > PspBoolMatrix;

// TODO: switch to CRTP-facade like technique like https://vector-of-bool.github.io/2020/06/13/cpp20-iter-facade.html
// or use a type-erasure / duck-typing strategy 
struct BoolMatrix {
	using F = PspBoolMatrix::value_type;
	using value_type = F; 
	using entry_t = PspBoolMatrix::entry_t;
	// using optional< pair< size_t, entry_t > >
	
	PspBoolMatrix& m;
	
	constexpr size_t n_rows() const { return m.n_rows(); };
	constexpr size_t n_cols() const { return m.n_cols(); };
	
	
	// Given row index i, j = col_of_low[i] yields the column j 
	// whose lowest non-zero row entry is i, if it exists
	// vector< std::optional< size_t > > row_to_col; 
	// vector< size_t > c_otc;  
	
	// Stores sorted ( original )
	// vector< size_t > col_otc;
	// vector< size_t > col_cto;
	// vector< optional< size_t > > Lr2c; // maps original low row indices --> original col indices
	// vector< optional< size_t > > Lc2r; // maps original col indices --> original low row indices
		
	BoolMatrix(PspBoolMatrix& pm) : m(pm) { //col_otc(pm.n_cols()), col_cto(pm.n_cols()), row_to_col(pm.n_cols()) 
		// std::iota(col_otc.begin(), col_otc.end(), 0);
		// std::iota(col_cto.begin(), col_cto.end(), 0);
		// 
		// // Initialize low entry cache
		// for (size_t j = 0; j < pm.n_cols(); ++j){
		// 	auto low_j = m.lowest_nonzero(j);
		// 	if (low_j){
		// 		const size_t ri = m.cto[low_j->first];     // original row index
		// 		Lr2c[ri] = std::make_optional(col_cto[j]); // store original column index (same in this case)
		// 		Lc2r[col_cto[j]] = std::make_optional(ri);
		// 	}
		// }
	};
	
	// Returns the lowest entry of column j 
	// { a.low(size_t(0)) } -> std::same_as< optional< pair< size_t, F > > >;
	auto low(size_t j) { 
		if (j >= n_cols()){ throw std::invalid_argument("column index out of range"); }
		// if (Lc2r[col_cto[j]]){
		// 	const size_t ri = Lc2r[col_cto[j]].value(); // original row index 
		// 	return(std::make_optional(pm.otc[ri], true));
		// }
		return m.lowest_nonzero(j);
	}

	// { a.low_index(size_t(0)) } -> std::same_as< optional< size_t > >;
	auto low_index(size_t j){
		auto le = low(j);
		return(le ? std::make_optional(le->first) : std::nullopt);
	}
	
	// { a.low_value(size_t(0)) } -> std::same_as< optional< F > >;
	auto low_value(size_t j){
		auto le = low(j);
		return(le ? std::make_optional(le->second) : std::nullopt);
	}
	// a.add_scaled_col(size_t(0), size_t(0), F(0));
	// (lambda * col(i)) + col(j) -> col(j)
	// auto add_scaled_col(size_t i, size_t j, size_t k, F lambda){
	// 	m.add_cols(i, j, k);
	// }
	
	// Use column s to cancel lowest entry of t, if it exists
	// { a.cancel_lowest(size_t(0), size_t(0)) } -> std::same_as< void >;
	void cancel_lowest(size_t t, size_t s){
		auto low_t = m.lowest_nonzero(t);
		auto low_s = m.lowest_nonzero(s);
		if (low_t && low_s && low_s->first == low_t->first){
			m.add_cols(s, t, t);
			// Rprintf("added columns %d to %d (s low: %d, t low: %d)\n", s, t, low_s->first, low_t->first);
		} else {
			Rprintf("added columns %d to %d (s low: %d, t low: %d)\n", s, t, low_s->first, low_t->first);
			throw std::invalid_argument("Unable to cancel the lowest entry!");
		}
	}
	
	// Zero's out column j
	auto clear_column(size_t j){
		if (m.columns.at(j)){
			m.columns[j]->clear();
		}
	}
	
	// TODO: add cancel_lowest(size_t i, size_t j, low_entry& s_low, low_entry& t_low){
	
	// { a.dim() } -> std::same_as< pair< size_t, size_t > >;
	auto dim() -> pair< size_t, size_t > {
		return std::make_pair(m.size[0], m.size[1]);
	} 
	
	// Given column j which has low row index 'j_low_index', find the column 'i' which has the same low row index
	// { a.find_low(size_t(0), size_t(0)) } -> std::same_as< std::optional< pair< size_t, F > > >; 
	auto find_low(size_t j, size_t j_low_index) -> std::optional< pair< size_t, F > >  {
		for (size_t i = 0; i < j; ++i){
			auto low_i = low(i);
			if (low_i && low_i->first == j_low_index){
				return(make_optional(make_pair(i, low_i->second)));	// Note this is (column index, field value), not (row index, field value)
			}
		}
		return(std::nullopt);
	}
	
	template <typename ... Args>
	void add_col(Args&& ... args){ m.add_col(std::forward<Args>(args)...); }
	
	// { a.swap_rows(size_t(0), size_t(0)) } -> std::same_as< void >;
	template <typename ... Args>
	void swap_rows(Args&& ... args){ m.swap_rows(std::forward<Args>(args)...); }
	
	// { a.swap_cols(size_t(0), size_t(0)) } -> std::same_as< void >;
	template <typename ... Args>
	void swap_cols(Args&& ... args){ m.swap_cols(std::forward<Args>(args)...); }
	
	// { a.swap(size_t(0), size_t(0)) } -> std::same_as< void >;
	template <typename ... Args>
	void swap(Args&& ... args){ m.swap(std::forward<Args>(args)...); }
	
	// { a.permute_rows(std::span< size_t >()) } -> std::same_as< void >;
	template <typename ... Args>
	void permute_rows(Args&& ... args){ m.permute_rows(std::forward<Args>(args)...); }
	
	// { a.permute_cols(std::span< size_t >()) } -> std::same_as< void >;
	template < typename ... Args >
	void permute_cols(Args&& ... args){ 
		m.permute_cols(std::forward<Args>(args)...);
	}
	
	template <typename ... Args>
	void permute(Args&& ... args){ m.permute(std::forward<Args>(args)...); }
	
	template <typename ... Args>
	void row(Args&& ... args){ m.row(std::forward<Args>(args)...); }
	
	// { a.column_empty(size_t(0)) } -> std::same_as< bool >;
	template <typename ... Args>
	bool column_empty(Args&& ... args){ return m.column_empty(std::forward<Args>(args)...); }
	
	template <typename ... Args>
	void column(Args&& ... args){ m.column(std::forward<Args>(args)...); }
	
	template <typename ... Args>
	void write_column(Args&& ... args){ m.write_column(std::forward<Args>(args)...); }
	
	template <typename ... Args>
	auto find_in_col(Args&& ... args) -> std::optional< pair< size_t, F > >{ 
		return m.find_in_col(std::forward<Args>(args)...); 
	}
		
	template <typename ... Args>
	void add_cols(Args&& ... args){ 
		// sm.print(Rcout);
		m.add_cols(std::forward<Args>(args)...); 
		// m.print(Rcout);
		// m.clean(0);
	}
	
	// unneeded
	template <typename ... Args>
	void scale_col(Args&& ... args){ }
	
	void add_scaled_col(size_t i, size_t j, size_t k, F lambda = true){
		if (i != k && j != k){ throw std::invalid_argument("i or j must equal k."); }
		add_cols(i,j,k);
		// m.print(Rcout);
	}
	
	// { a(size_t(0), size_t(0)) } -> std::same_as< F >; 
	template <typename ... Args>
	F operator()(Args&& ... args){ return m(std::forward<Args>(args)...); }
	
	template <typename ... Args>
	void assign_column(Args&& ... args){
		m.assign_column(std::forward<Args>(args)...);
	}
	
};


// ----- BEGIN MOVES ------

#include "Rcpp.h"
using namespace Rcpp;

using std::vector;

// Restore right: given column indices i \in [bi, ei) to restore, apply the donor concept to given indices 
// Postcondition: dr and dv are populated as donor columns
template< PermutableMatrix Matrix, typename Iter >
int restore_right(Matrix& R, Matrix& V, Iter bi, const Iter ei, vector< typename Matrix::entry_t >& dr, vector< typename Matrix::entry_t >& dv){
	using entry_t = typename Matrix::entry_t;
	int nr = 0; 
	const size_t ne = std::distance(bi, ei); 
	if (ne == 0){ return(0); }
	vector< entry_t > dr_new, dv_new; 
	auto d_low_index = R.low_index(*bi);
	auto d_low_index_new = std::optional< size_t >{ std::nullopt };
	dr.clear(); dv.clear();
	R.write_column(*bi, std::back_inserter(dr));
	V.write_column(*bi, std::back_inserter(dv));
	
	
	// Restore columns 
	for (auto b = std::next(bi); b != ei; ++b){
		size_t k = *b;
		d_low_index_new = R.low_index(k);
		
		// Condition when to replace donor columns
		const bool overwrite_donors = d_low_index_new.has_value() 
			&& (d_low_index && d_low_index_new.value() < d_low_index.value())
			|| (!d_low_index_new.has_value());
		
		// Determine if columns will need to be saved
		if (overwrite_donors){
			// Rcout << "overwriting donor w/ column " << k << std::endl; 
			R.write_column(k, std::back_inserter(dr_new));
			V.write_column(k, std::back_inserter(dv_new));
		}
		
		// Do the reductions 
		Rcout << "Adding donor column to column " << k << std::endl; 
		R.add_col(dr.begin(), dr.end(), k);
		V.add_col(dv.begin(), dv.end(), k);
		++nr;
		
		// Rcout << std::endl; 
		// Rcout << "donor column: ";
		// for (auto e: dr){ Rprintf("(%d, %d) ", e.first, e.second); }
		// Rcout << std::endl; 
		
		// Rcout << "column " << k << ": ";
		// R.column(k, [](auto ri, auto val){ Rprintf("(%d, %d) ", ri, val); });
		// Rcout << std::endl; 
		
		// Save the new donor columns if needed
		if (overwrite_donors){
			d_low_index = d_low_index_new;
			dr = dr_new;
			dv = dv_new;
			dr_new.clear();
			dv_new.clear();
		}
	}
	return(nr);
}

// Restore right
// template< PermutableMatrix Matrix, typename Iter >
// int restore_right(Matrix& R, Matrix& V, Iter bi, const Iter ei, vector< typename Matrix::entry_t >& dr, vector< typename Matrix::entry_t >& dv){
// 	using entry_t = typename Matrix::entry_t;
// 
// 	// Start with donor columns
// 	size_t donor_idx = *bi;
// 	auto d_low_index = R.low_index(donor_idx);
// 	dr.clear(); dv.clear();
// 	R.write_column(donor_idx, std::back_inserter(dr));
// 	V.write_column(donor_idx, std::back_inserter(dv));
// 	if (std::distance(bi, ei) <= 1){ return 0; }
// 
// 	// Apply the donor concept
// 	int nr = 0; 
// 	std::advance(bi, 1);
// 	auto new_dr = vector< entry_t >(), new_dv = vector< entry_t >();
// 	for (size_t k = *bi; bi != ei; ++bi){
// 		auto new_low_index = R.low_index(k);
// 
// 		// Copy column into new donor column
// 		new_dr.clear(); new_dv.clear();
// 		R.write_column(donor_idx, std::back_inserter(new_dr));
// 		V.write_column(donor_idx, std::back_inserter(new_dv));
// 
// 		// Add the columns
// 		R.add_cols(donor_idx, k, k); 
// 		V.add_cols(donor_idx, k, k);
// 		++nr;
// 
// 		// Replace the donor columns
// 		if (d_low_index && new_low_index.has_value() && *d_low_index > *new_low_index){
// 			d_low_index = std::move(new_low_index);
// 			dr.swap(new_dr);
// 			dv.swap(new_dv);
// 		} else if (d_low_index && !new_low_index.has_value()){
// 			donor_idx = k;
// 			d_low_index = std::nullopt;
// 			dr.swap(new_dr);
// 			dv.swap(new_dv);
// 		}
// 	}
// 	return(nr);
// }

// Restore left
// Restores columns [b, e) of matrices R and V to valid states 
template< PermutableMatrix Matrix, typename Iter >
int restore_left(Matrix& R, Matrix& V, Iter b, const Iter e){
	using F = typename Matrix::value_type;
	using entry_t = typename Matrix::entry_t;
	using low_entry = typename std::optional< pair< size_t, F > >; 
	const size_t n = std::distance(b, e);
	if (n <= 1){ return 0; }
	
	struct LowPair {
		size_t column_index; 
		optional< pair< size_t, F > > low_entry;
	};
	
	// First sort the (column index, low entry) elements in increasing order of row indices
	// Settle ties by sorting by increasing column index 
	auto K = vector< LowPair >();
	K.reserve(std::distance(b, e));
	std::for_each(b, e, [&](auto j){ K.push_back({ j, R.low(j) }); });
	const auto low_cmp = [](const LowPair& le1, const LowPair& le2){
		if (!le1.low_entry && !le2.low_entry){ return(false); }
		if (!le1.low_entry){ return(true); }
		if (!le2.low_entry){ return(false); }
		if (le1.low_entry->first == le2.low_entry->first){ return(le1.column_index < le2.column_index); }
		return(le1.low_entry->first < le2.low_entry->first); 
	};
	std::sort(K.begin(), K.end(), low_cmp);
	// const bool is_ordered = std::is_sorted(K.begin(), K.end(), low_cmp);
	// Rprintf(" is ordered? %s \n", is_ordered ? "TRUE" : "FALSE");
	// for (auto le: K){
	// 	Rprintf("(c:%d,r:%d) ", le.column_index, le.low_entry ? le.low_entry->first : -1);
	// }
	// Rcout << std::endl; 
	
	// Perform the reductions on (l, r)
	int nr = 0; 
	LowPair lp = K[n-2], rp = K[n-1];
	while(lp.low_entry && rp.low_entry){
		// Rcout << "r: ";
		// for (auto le: K){
		// 	// Rprintf("(c:%d,r:%d) ", le.column_index, le.low_entry ? le.low_entry->first : -1);
		// 	Rprintf("%d ", le.low_entry ? le.low_entry->first : -1);
		// }
		// Rcout << std::endl;
		// Rcout << "c: ";
		// for (auto le: K){
		// 	// Rprintf("(c:%d,r:%d) ", le.column_index, le.low_entry ? le.low_entry->first : -1);
		// 	Rprintf("%d ", le.column_index);
		// }
		// Rcout << std::endl;
		if (lp.low_entry->first == rp.low_entry->first){
			// Rprintf("(%d -> %d) \n", lp.column_index, rp.column_index);
			R.cancel_lowest(rp.column_index, lp.column_index);
			V.add_scaled_col(lp.column_index, rp.column_index, rp.column_index, 1);
			K.pop_back();
			++nr;
			
			// column r is reduced: check to see if need to reinsert updated entry for r
			LowPair updated_lp = { rp.column_index, R.low(rp.column_index) };
			K.insert(std::upper_bound(K.begin(), K.end(), updated_lp, low_cmp), updated_lp);
			// if (updated_lp.low_entry){
			// 	auto it = std::find_if(K.begin(), K.end(), [&updated_lp](const LowPair& p){
			// 		return(p.low_entry && p.low_entry->first == updated_lp.low_entry->first);
			// 	});
			// 	if (it != K.end()){
			// 		K.insert(std::upper_bound(K.begin(), K.end(), updated_lp, low_cmp), updated_lp);
			// 	}
			// }
		} else {
			// Low entries don't match; remove highest row index and continue
			K.pop_back();
		}
		// Update the left and right column entries
		if (K.size() <= 1){ break; }
		rp = K.back();
		lp = *std::prev(std::end(K), 2);
	}
	return(nr);
}

[[nodiscard]]
inline auto move_right_permutation(size_t i, size_t j, const size_t n) -> std::vector< size_t > {
  if (i > j){ throw std::invalid_argument("invalid");}
  std::vector< size_t > v(n);
  std::iota(v.begin(), v.end(), 0);
  std::rotate(v.begin()+i, v.begin()+i+1, v.begin()+j+1);
  return(v);
}

[[nodiscard]]
inline auto move_left_permutation(size_t i, size_t j, const size_t n) -> std::vector< size_t >{
  if (i < j){ throw std::invalid_argument("invalid");}
  std::vector< size_t > v(n);
  std::iota(v.begin(), v.end(), 0);
  std::rotate(v.rbegin()+(n-(i+1)), v.rbegin()+(n-i), v.rbegin()+(n-j));
  return(v);
}

template< PermutableMatrix Matrix >
int move_left_local(Matrix& R1, Matrix& V1, Matrix& R2, Matrix& V2, const size_t i, const size_t j){
	using entry_t = typename Matrix::entry_t; 
	if (i == j){ return 0; }
	if (i < j){ throw std::invalid_argument("Invalid pair (i,j) given (i < j)."); }
	const size_t nc1 = V1.dim().second;
	const size_t nc2 = V2.dim().second;
	
	// Remove non-zero row entries in column i of V1 
	auto K = vector< size_t >();
	for (size_t r = i; r > j; --r){
		const size_t ri = r - 1; // note: ri < i
		std::optional< entry_t > next_low = V1.find_in_col(i, ri);
		if (next_low && next_low.value().second != 0){
			V1.add_scaled_col(ri, i, i); 
			R1.add_scaled_col(ri, i, i);
			K.push_back(ri+1); // Because all row indices will be shifted up one
		}
	}

	// Apply permutation PR1P^T, PV1P^T
	vector< size_t > p = move_left_permutation(i, j, nc1);
	R1.permute_cols(p.begin(), p.end());
	V1.permute(p.begin(), p.end());
	
	// Restore R1 to a valid state
	// Rcout << "Restoring R1: ";
	K.push_back(j);
	// for (auto ki: K){ Rcout << ki << ", "; }
	int nr = restore_left(R1, V1, K.begin(), K.end());
	// Rcout << std::endl;
	
	// Fix columns in R2, if necessary 
	// low_J <- apply(R2, 2, low_entry)
	// k <- which(low_J == i)
	// if (length(k) != 0){
	// 	k_low <- (j-1) + low_entry(R2[seq(j, i-1L),k])
	// 	J <- sort(c(k, which(low_J <= k_low))) ## can we do better?
	// 	R2 <- permute_move(R2, i = i, j = j, dims = "rows") 
	// 	res2 <- restore_left(R = R2, V = V2, J = J)
	
	// To remove once low caches are made
	// Need: columns whose low row indices are <= second lowest index of column k
	// where column k is the column whose lowest row entry is i
	// Rcout << "Restoring R2: ";
	auto J = vector< size_t >(nc2);
	std::iota(J.begin(), J.end(), 0);
	R2.permute_rows(p.begin(), p.end());
	nr += restore_left(R2, V2, J.begin(), J.end());
	
	// auto J = vector< size_t >();
	// std::optional< size_t > k_col;
	// size_t k_low = 0;
	// for (size_t ki = 0; ki < nc2; ++ki){
	// 	auto k_low = R2.low(ki);
	// 	if (k_low && k_low->first == i){
	// 		J.push_back(ki);
	// 		break;
	// 	}
	// }
	
	// if (J.size() >= 1){
	// 	for (size_t ki = 0; ki < nc2; ++ki){
	// 		auto k_low = R2.low(ki);
	// 		if (k_low && k_low->first < i){
	// 			J.push_back(ki);
	// 		}
	// 	}
	// 	std::sort(J.begin(), J.end());
	// 	R2.permute_rows(p.begin(), p.end());
	// 	Rcout << "Restoring R2: " << std::endl;
	// 	nr += restore_left(R2, V2, J.begin(), J.end());
	// } else {
	// 	R2.permute_rows(p.begin(), p.end());
	// }
	// Rcout << std::endl;
	
	
	return(nr + K.size());
}

template< PermutableMatrix Matrix >
int move_right_local(Matrix& R1, Matrix& V1, Matrix& R2, Matrix& V2, const size_t i, const size_t j){
	using entry_t = typename Matrix::entry_t; 
	if (i == j){ return 0; }
	if (i > j){ throw std::invalid_argument("Invalid pair (i,j) given (i > j)."); }
	const size_t nc1 = R1.dim().second;
	const size_t nc2 = R2.dim().second;
	
	// Collect indices I
	auto I = vector< size_t >();
	V1.row(i, [&](auto col_idx, auto v){
		if ((col_idx >= i) && (col_idx <= j)){ I.push_back(col_idx); }
	});
	Rcout << "I: ";
	for (auto jj: I){ Rcout << jj << ", "; }
	Rcout << std::endl;

	// Collect indices J
	auto J = vector< size_t >();
	for (size_t c = 0; c < nc2; ++c){
		auto low_idx = R2.low_index(c);
		if (low_idx && (*low_idx) >= i && (low_idx) <= j && R2(i,c) != 0){
			J.push_back(c);
		}
	}
	Rcout << "J: ";
	for (auto jj: J){ Rcout << jj << ", "; }
	Rcout << std::endl;

	// Restore invariants to columns affected by I
	auto dr1 = vector< entry_t >();
	auto dv1 = vector< entry_t >();
	Rcout << "Restoring R1: ";
	int nr = restore_right(R1, V1, I.begin(), I.end(), dr1, dv1); // don't wrap in if condition
	Rcout << std::endl;
	
	// Restore invariants to columns/rows affected by J
	if (!J.empty()){
		auto dr2 = vector< entry_t >();
		auto dv2 = vector< entry_t >();
		Rcout << "Restoring R2: ";
		nr += restore_right(R2, V2, J.begin(), J.end(), dr2, dv2);
		Rcout << std::endl;
	}

	// Apply permutations
	const vector< size_t > p = move_right_permutation(i, j, nc1);
	Rcout << "move permutation";
	for (auto ii: p){ Rcout << ii << ", "; }
	Rcout << std::endl;
	const vector< size_t > q = inverse_permutation(p.begin(), p.end());
	// Rcout << "inverse move permutation";
	// for (auto ii: q){ Rcout << ii << ", "; }
	// Rcout << std::endl;
	
	R1.permute_cols(p.begin(), p.end());
	V1.permute(p.begin(), p.end());
	R2.permute_rows(p.begin(), p.end());
	
	// Rcout << "dv: ";
	// for (auto e: dv1){ Rcout << e.first << " : " << e.second << ", "; }
	// Rcout << std::endl;
	// Rcout << "dr: ";
	// for (auto e: dr1){ Rcout << e.first << " : " << e.second << ", "; }
	// Rcout << std::endl;
	
	// Apply permutations and assign donors
	// Note: this is not counted towards 'nr' because this can be made O(1) w/ move semantics
	// and because no field operations are counted
	const auto apply_q = [&q](entry_t& e){ return(std::make_pair(q[e.first], e.second)); };
	std::transform(dv1.begin(), dv1.end(), dv1.begin(), apply_q);
	// note dr1 rows need not be modified, only dr2 
	
	// apply_permutation(dr1.begin(), dr1.end(), p.begin());
	// apply_permutation(dv1.begin(), dv1.end(), p.begin());
	
	// Rcout << "dv: ";
	// for (auto e: dv1){ Rcout << e.first << " : " << e.second << ", "; }
	// Rcout << std::endl;
	// Rcout << "dr: ";
	// for (auto e: dr1){ Rcout << e.first << " : " << e.second << ", "; }
	// Rcout << std::endl;
	
	
	R1.assign_column(j, dr1.begin(), dr1.end());
	V1.assign_column(j, dv1.begin(), dv1.end());
	return(nr);
}


// Performs a sequence of moves (i,j) on the pair (R1, V1), correcting (R2, V2) as necessary. 
template< PermutableMatrix Matrix, typename Iter, typename Lambda >
int move_schedule_local(Matrix& R1, Matrix& V1, Matrix& R2, Matrix& V2, Iter sb, const Iter se, Lambda f){
	using entry_t = typename Matrix::entry_t;
	const size_t nr1 = R1.dim().first;
	const size_t nc1 = R1.dim().second;
	const size_t nr2 = R2.dim().first;
	if (nc1 != nr2){ throw std::invalid_argument("Number of rows in R2 must match number of columns in R1."); }
	if (nc1 == 0 || nr1 == 0){ return 0; }
	if (std::distance(sb,se) < 2 || std::distance(sb,se) % 2 != 0){ throw std::invalid_argument("Pairs of indices must be passed."); }

	int nr = 0; 
	for (size_t i, j; sb != se; sb += 2){
		i = *sb, j = *(sb+1);
		if (i == j){ continue; }
		if (i >= nc1 || j >= nc1){  throw std::invalid_argument("Invalid pairs (i,j) passed ( i or j >= number of columns )."); }
		if (i < j){
			nr += move_right_local(R1, V1, R2, V2, i, j);
		} else {
			nr += move_left_local(R1, V1, R2, V2, i, j);
		}
		f();
		Rcpp::checkUserInterrupt();
	}
	return(nr);
}

// Given a contiguous schedule of integer indices (i,j) via [sb, se), performs a sequence of move operations
// on the iterator range
template< PermutableMatrix Matrix, typename Iter, typename Lambda >
void move_schedule_full(Matrix& R, Matrix& V, Iter sb, const Iter se, Lambda f){
	using entry_t = typename Matrix::entry_t;
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

// ----- END MOVES -----



// [[Rcpp::export]]
int reduce_pspbool(SEXP D_ptr, SEXP V_ptr, bool show_progress=true) {
	BoolMatrix D = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(D_ptr)); 
	BoolMatrix V = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(V_ptr)); 
	auto indices = vector< size_t >(D.dim().second);
	std::iota(indices.begin(), indices.end(), 0);
	reduction_stats[0] = 0; 
	Progress p(indices.size(), show_progress);
	size_t cc = 0; 
	const auto P = [&](){ 
		p.increment(); 
		if ((++cc % 64) == 0){
			D.m.clean(0); V.m.clean(0);
		}
	};
	pHcol(D, V, indices.begin(), indices.end(), P);
	return(reduction_stats[0]);
}

// [[Rcpp::export]]
int reduce_local_pspbool(SEXP D1_ptr, SEXP V1_ptr, SEXP D2_ptr, SEXP V2_ptr, bool clearing, bool show_progress=true) {
	BoolMatrix D1 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(D1_ptr)); 
	BoolMatrix V1 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(V1_ptr)); 
	BoolMatrix D2 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(D2_ptr)); 
	BoolMatrix V2 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(V2_ptr)); 
	auto d1_ind = vector< size_t >(D1.dim().second);
	auto d2_ind = vector< size_t >(D2.dim().second);
	std::iota(d1_ind.begin(), d1_ind.end(), 0);
	std::iota(d2_ind.begin(), d2_ind.end(), 0);
	
	Progress p(d1_ind.size() + d2_ind.size(), show_progress);
	size_t cc = 0; 
	const auto P = [&](){ 
		p.increment(); 
		if ((++cc % 64) == 0){
			D1.m.clean(0); D2.m.clean(0);
			V1.m.clean(0); V2.m.clean(0);
		}
	};
	reduction_stats[0] = 0; 
	if (clearing){
		pHcol_local< true >(D1, V1, D2, V2, d1_ind.begin(), d1_ind.end(), d2_ind.begin(), d2_ind.end(), P);
	} else {
		pHcol_local< false >(D1, V1, D2, V2, d1_ind.begin(), d1_ind.end(), d2_ind.begin(), d2_ind.end(), P);
	}
	return(reduction_stats[0]);
}

//[[Rcpp::export]]
int simulate_vineyard_pspbool(SEXP R_ptr, SEXP V_ptr, IntegerVector schedule, Nullable< Function > f = R_NilValue, bool show_progress=true){
	BoolMatrix R = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(R_ptr)); 
	BoolMatrix V = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(V_ptr)); 
	if (f.isNotNull()){
		Function fun(f);
		size_t cc = 0; 
		transpose_schedule_full(R, V, schedule.begin(), schedule.end(), [&](size_t status){
			fun(status);
			if ((++cc % 1024) == 0){
				R.m.clean(0); V.m.clean(0);
			}
		});
		return(0);
	} else {
		const std::array< size_t, 10 > costs = { 0, 2, 2, 0, 1, 2, 0, 2, 0, 1 };
		size_t nc = 0, cc = 0; 
		Progress p(schedule.size(), show_progress);
		const auto record_costs = [&costs, &nc, &cc, &R, &V, &p](size_t status){ 
			nc += costs[status]; 
			if ((++cc % 1024) == 0){
				R.m.clean(0); V.m.clean(0);
			}
			p.increment();
		};
		transpose_schedule_full(R, V, schedule.begin(), schedule.end(), record_costs);
		return(nc);
	}
}

//[[Rcpp::export]]
int move_schedule_local(SEXP r1, SEXP v1, SEXP r2, SEXP v2, IntegerVector schedule, Nullable< Function > f = R_NilValue){
	BoolMatrix R1 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(r1)); 
	BoolMatrix V1 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(v1)); 
	BoolMatrix R2 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(r2)); 
	BoolMatrix V2 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(v2)); 
	int nr = 0; 
	if (f.isNotNull()){
		Function fun(f);
		nr = move_schedule_local(R1, V1, R2, V2, schedule.begin(), schedule.end(), [&](){
			fun();
		});
	} else {
		const auto do_nothing = [](){ return; };
		nr = move_schedule_local(R1, V1, R2, V2, schedule.begin(), schedule.end(), do_nothing);
	}
	return(nr);
}


// Rcpp Module exports
SEXP as_XPtr(PspBoolMatrix* M){
  Rcpp::XPtr< PspBoolMatrix > p(M, false);
  return(p);
}

void print(PspBoolMatrix* M){
	M->print(Rcout);
}

int low_entry(PspBoolMatrix* M, size_t j){
	auto entry = M->lowest_nonzero(j); 
	return (bool(entry) ? entry->first : -1 );
}

bool at(PspBoolMatrix* M, size_t i, size_t j){
	auto& Mr = *M;
	return(Mr(i,j));
}

IntegerVector find_low(PspBoolMatrix* M, size_t j){
	auto low_j = low_entry(M, j);
	if (low_j != -1){
		for (size_t i = 0; i < j; ++i){
			auto low_i = low_entry(M, i);
			if (low_i != -1 && low_i == low_j){ 
				return(IntegerVector::create(i, low_i));	// Note this is (column index, field value), not (row index, field value)
			}
		}
	}
	return(IntegerVector::create(-1, -1));
}

// Szudziks pairing function. Takes as input two unsigned integral types (a, b), and uniquely
// maps (a, b) to a number c, where c is possibly a different integral type
constexpr inline size_t rosenburg_pair(unsigned int a, unsigned int b){
  return static_cast< size_t >(a >= b ? a * a + a + b : a + b * b);
}

// Assumes a < b!
constexpr inline size_t rosenburg_unpair(unsigned int a, unsigned int b){
  return static_cast< size_t >(a + b * b);
}

IntegerMatrix pivots(PspBoolMatrix* M){
	const size_t nc = M->n_cols();
	vector< size_t > pairs; 
	pairs.reserve(nc);
	for (size_t j = 0; j < nc; ++j){
		auto entry = M->lowest_nonzero(j); 
		if (entry){
			pairs.push_back(entry->first);
			pairs.push_back(j);
		}
	}
	IntegerMatrix pivots(2, pairs.size()/2, pairs.begin());
	return (pivots);
}


//' @export PspBoolMatrix
RCPP_EXPOSED_CLASS_NODECL(PspBoolMatrix)
RCPP_MODULE(PspBoolMatrix) {
	using namespace Rcpp; 
	class_< PspBoolMatrix >( "PspBoolMatrix" )
		.constructor< size_t, size_t >()
  	.method( "as_XPtr", &as_XPtr)
  	.field("nnz", &PspBoolMatrix::nnz )
  	.field("cto", &PspBoolMatrix::cto)
  	.field("otc", &PspBoolMatrix::otc)
  	.property("n_rows", &PspBoolMatrix::n_rows )
  	.property("n_cols", &PspBoolMatrix::n_cols )
  	.method("as.Matrix", &PspBoolMatrix::to_sparseMatrix)
  	.method("swap_rows", &PspBoolMatrix::swap_rows)
  	.method("swap_cols", &PspBoolMatrix::swap_cols)
  	.method("insert", &PspBoolMatrix::insert)
  	.method("remove", &PspBoolMatrix::remove)
  	.method("clean", &PspBoolMatrix::clean)
  	.method("add_columns", &PspBoolMatrix::add_cols)
  	.method("construct", &PspBoolMatrix::construct)
  	.method("print", &print)
  	// .method("permute_rows", &PspIntMatrix::permute_rows)
    // .method("permute_rows", &PspIntMatrix::permute_rows)
    .method("permute_rows", (void (PspBoolMatrix::*)(const vector< size_t >))(&PspBoolMatrix::permute_rows))
    // .method("permute_cols", &PspBoolMatrix::permute_cols)
    .method("permute_cols", (void (PspBoolMatrix::*)(const vector< size_t >))(&PspBoolMatrix::permute_cols))
  	// .method("permute_cols", &PspIntMatrix::permute_cols)
    .method("add_rows", &PspBoolMatrix::add_rows)
    .method("low_entry", low_entry)
    .method("find_low", find_low)
    .method("at", at)
    .method("submatrix", &PspBoolMatrix::submatrix)
    .method("pivots", &pivots)
	;
}

