#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "reduction.h"

// Sparse Matrix class using Armadillo
template< typename Value >
struct ReducibleSpMat {
	using value_type = Value;
	arma::SpMat< value_type >& m;
	
	static bool equals_zero(value_type val) {
		if constexpr (std::is_integral_v< value_type >){
			return val == 0; 
		} else {
			return std::abs(val) <= std::numeric_limits< value_type >::epsilon();
		}
	}
	ReducibleSpMat(arma::SpMat< value_type >& d) : m(d){};
	
	void clear_column(size_t j){
		m.col(j) = arma::SpMat< double >(m.n_rows, 1);
	}
	auto dim() const -> pair< size_t, size_t >{ return make_pair(m.n_rows, m.n_cols); }
	
	// Get the lowest non-zero entry at column j
	auto low(const size_t j) const -> optional< pair< size_t, value_type > > {
		auto b = m.begin_col(j);
		auto e = m.end_col(j);
		if (b == e){
			return std::nullopt;
		} else {
			auto entry = std::prev(e);
			return(std::make_optional(std::make_pair(entry.row(), *entry)));
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
		if (i != k && j != k){ throw std::invalid_argument("i or j must equal k."); }
		const size_t s = i == k ? j : i; 
		m.col(k) += m.col(s);
	}
	
	// Adds (lambda * col(i)) + col(j) -> col(j)
	void scale_col(size_t i, double lambda){
		m.col(i) = m.col(i)*lambda;
	}
	
	void add_scaled_col(size_t i, size_t j, size_t k, double lambda){
		if (i != k && j != k){ throw std::invalid_argument("i or j must equal k."); }
		const size_t s = i == k ? j : i; 
		m.col(k) += m.col(s)*lambda;
	}
	
	template< typename Lambda>
	void column(size_t j, Lambda&& f){
		auto b = m.begin_col(j);
		auto e = m.end_col(j);
		for (; b != e; ++b){
			f(b.row(), *b);
		}
	}
	
	template < typename Iter >
	void assign_column(size_t c, Iter b, const Iter e){
		m.col(c) = arma::SpMat< double >(m.n_rows, 1);
		for (; b != e; ++b){
			m.at(b->first, c) = b->second;
		}
	}
	
	bool column_empty(size_t j) const {
		return(m.begin_col(j) == m.end_col(j));
	}

	// Cancels lowest entry of column j by setting its value = 0 using column i
	void cancel_lowest(size_t j, size_t i){
		if (column_empty(j) || column_empty(i)){ return; }
		auto lambda = -(*low_value(j))/(*low_value(i));
		m.col(j) += lambda*m.col(i);
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
Rcpp::List reduce_arma(arma::sp_mat& D, arma::sp_mat& v) {
  auto R = ReducibleSpMat< double >(D);
	auto V = ReducibleSpMat< double >(v);	
	auto indices = std::vector< size_t >(D.n_cols-1);
	std::iota(indices.begin(), indices.end(), 1);
	pHcol(R, V, indices.begin(), indices.end());
	return Rcpp::List::create(Rcpp::_["R"]=Rcpp::wrap(R.m), Rcpp::_["V"]=Rcpp::wrap(V.m));
}

// [[Rcpp::export]]
Rcpp::List reduce_local_arma(arma::sp_mat& D1, arma::sp_mat& v1, arma::sp_mat& D2, arma::sp_mat& v2, bool clearing) {
  auto R1 = ReducibleSpMat< double >(D1);
	auto V1 = ReducibleSpMat< double >(v1);
	auto R2 = ReducibleSpMat< double >(D2);
	auto V2 = ReducibleSpMat< double >(v2);	
	auto d1_ind = std::vector< size_t >(D1.n_cols-1);
	auto d2_ind = std::vector< size_t >(D2.n_cols-1);
	std::iota(d1_ind.begin(), d1_ind.end(), 1);
	std::iota(d2_ind.begin(), d2_ind.end(), 1);
	pHcol_local(R1, V1, R2, V2, d1_ind.begin(), d1_ind.end(), d2_ind.begin(), d2_ind.end());
	auto out = Rcpp::List::create(
		Rcpp::_["R1"]=Rcpp::wrap(R1.m), Rcpp::_["V1"]=Rcpp::wrap(V1.m),
		Rcpp::_["R2"]=Rcpp::wrap(R2.m), Rcpp::_["V2"]=Rcpp::wrap(V2.m)
	);
	return(out);
}
