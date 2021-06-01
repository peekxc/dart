#include <Rcpp.h>
using namespace Rcpp;

#include "reduction.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Dense Matrix class 
struct DenseNumericMatrix {
	using value_type = double;
	NumericMatrix& m;
	
	DenseNumericMatrix(NumericMatrix& d) : m(d){};
	// DenseNumericMatrix(NumericMatrix&& d) : m(std::move(d)){};
	
	void clear_column(size_t j){
		m.column(j) = NumericVector();
	}
	
	static bool equals_zero(double val) {
		return std::abs(val) <= std::numeric_limits<double>::epsilon();
	}
	
	template < typename Iter >
	void assign_column(size_t c, Iter b, const Iter e){
		Rcpp::NumericVector col = Rcpp::NumericVector();
		for (; b != e; ++b){
			col.push_back(b->second);
		}
		m(Rcpp::_, c) = col;
	}
	
	auto dim() const -> pair< size_t, size_t >{ return make_pair(m.nrow(), m.ncol()); }
	
	// Get the lowest non-zero entry at column j
	auto low(const size_t j) const -> optional< pair< size_t, value_type > > {
		NumericVector col_j = m.column(j);
		auto vj = std::vector< double >(col_j.begin(), col_j.end());
		auto entry_it = std::find_if(vj.rbegin(), vj.rend(), [](auto val){ return !equals_zero(val); });
		if (entry_it != vj.rend()){
			const auto row_idx = vj.size() - std::distance(vj.rbegin(), entry_it) - 1;
			return(make_optional(std::make_pair(row_idx, *entry_it)));
		} else {
			return(std::nullopt);
		}
	}
	// Get the lowest non-zero entry at column j
	auto low_index(const size_t j) const -> optional< size_t > {
		auto le = low(j);
		return(le ? std::make_optional(le->first) : std::nullopt);
	}
	auto low_value(const size_t j) const -> optional< value_type > {
		auto le = low(j);
		return(le ? std::make_optional(le->second) : std::nullopt);
	}
	
	void add_cols(size_t i, size_t j, size_t k){
		if (i != k && j != k){ throw std::invalid_argument("i,j must equal k."); }
		const size_t s = i == k ? j : i;
		const NumericMatrix::Column& s_col = m(_, s);
		for (size_t r = 0; r < m.nrow(); ++r){
			m(r,k) += s_col[r];
		}
	}
	
	template < typename Lambda >
	void column(size_t c, Lambda&& f){
		if (c >= m.ncol()){ return; }
		const NumericMatrix::Column& c_col = m(_, c);
		for (size_t i = 0; i < c_col.size(); ++i){
			f(i, c_col[i]);
		} 
	}
	
	bool column_empty(size_t j) const {
		NumericVector col_j = m.column(j);
		return std::all_of(col_j.begin(), col_j.end(), equals_zero);
	}
	
	void scale_col(size_t i, double lambda){
		m.column(i) = lambda*m.column(i);
	}
	void add_scaled_col(size_t i, size_t j, size_t k, double lambda){
		if (i != k && j != k){ throw std::invalid_argument("i or j must equal k."); }
		const size_t s = i == k ? j : i; 
		m.column(k) = m.column(k) + lambda*m.column(s);
	}
	
	// Cancels lowest entry of column j by setting its value = 0 using column i 
	void cancel_lowest(size_t j, size_t i){
		if (column_empty(j)){ return; }
		auto lambda = -(*low_value(j))/(*low_value(i));
		m.column(j) = m.column(j) + lambda*m.column(i);
		// add_scaled_col(i,j, );
	}
	
	// Finds a column index i + its lowest entry low_i such that i < j and low_i == low_j, if it exists
	auto find_low(size_t j, size_t low_j) const -> optional< pair< size_t, value_type > >{
		for (size_t i = 0; i < j; ++i){
			auto low_i = low(i);
			if (low_i && low_i->first == low_j){ 
				return(make_pair(i, low_i->second));	// Note this is (column index, field value), not (row index, field value)
			}
		}
		return(std::nullopt);
	}
};



// [[Rcpp::export]]
List reduce_local_dense(const NumericMatrix& D1, const NumericMatrix& v1, const NumericMatrix& D2, const NumericMatrix& v2, bool show_progress=true) {
	NumericMatrix d1c = clone(D1);
	NumericMatrix v1c = clone(v1);
	NumericMatrix d2c = clone(D2);
	NumericMatrix v2c = clone(v2);
  auto R1 = DenseNumericMatrix(d1c);
	auto V1 = DenseNumericMatrix(v1c);
	auto R2 = DenseNumericMatrix(d2c);
	auto V2 = DenseNumericMatrix(v2c);	
	auto d1_ind = std::vector< size_t >(D1.ncol()-1);
	auto d2_ind = std::vector< size_t >(D2.ncol()-1);
	std::iota(d1_ind.begin(), d1_ind.end(), 1);
	std::iota(d2_ind.begin(), d2_ind.end(), 1);
	Progress p(d1_ind.size() + d2_ind.size(), show_progress);
	const auto P = [&p](){ p.increment(); };
	pHcol_local(R1, V2, R2, V2, d1_ind.begin(), d1_ind.end(), d2_ind.begin(), d2_ind.end(), P);
	return List::create(_["R1"]=R1.m, _["V1"]=V1.m, _["R2"]=R2.m, _["V2"]=V2.m);
}

// [[Rcpp::export]]
List reduce_dense(const NumericMatrix& D, const NumericMatrix& v, bool show_progress=true) {
	NumericMatrix dd = clone(D);
	NumericMatrix vv = clone(v);
  auto R = DenseNumericMatrix(dd);
	auto V = DenseNumericMatrix(vv);	
	auto indices = std::vector< size_t >(D.ncol()-1);
	std::iota(indices.begin(), indices.end(), 1);
	Progress p(indices.size(), show_progress);
	const auto P = [&p](){ p.increment(); };
	pHcol(R, V, indices.begin(), indices.end(), P);
	return List::create(_["R"]=R.m, _["V"]=V.m);
}


