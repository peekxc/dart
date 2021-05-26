#include <vector>
#include <algorithm>
#include <concepts>
#include <iterator>
#include <cstddef>
#include <array>
#include <optional>
#include <tuple>
#include <sstream>

using std::size_t;
using std::tuple;
using std::pair; 
using std::vector; 
using std::unique_ptr;
using std::array; 
using std::optional;
using std::make_tuple;
using std::make_pair;
using std::make_optional;

#include <Rcpp.h>
using namespace Rcpp;

// Sparse Matrix representation that is robust to row/column permutations 
template< typename T, class BinaryOperation >
class PspMatrix {
public: 
	using entry_t = pair< size_t, T >;
	using column_t = vector< entry_t >;
	using entries_ptr = unique_ptr< vector< entry_t > >;
	using value_type = T;
	
	vector< entries_ptr > columns; // non-zero entries
	vector< size_t > cto;  // current-to-original map
	vector< size_t > otc;	 // original-to-current map
	
	BinaryOperation add; 
	const array< size_t, 2 > size; 
	size_t nnz = 0;
	constexpr size_t n_rows() const { return size.size() < 1 ? 0 : size[0]; };
	constexpr size_t n_cols() const { return size.size() < 2 ? 0 : size[1]; };
	
	PspMatrix(const size_t m, const size_t n) : columns(n), cto(m), otc(m), size({ m, n }){
		for (size_t i = 0; i < n; ++i){
			columns[i] = std::make_unique< vector< entry_t > >();
		}
		std::iota(std::begin(cto), std::end(cto), 0);
		std::iota(std::begin(otc), std::end(otc), 0);
	}; 
	PspMatrix(PspMatrix&&) = default; 									// move construct
	PspMatrix& operator=(PspMatrix&&) = default; 
	PspMatrix(const PspMatrix&) = default;							// copy 
	PspMatrix& operator=(const PspMatrix&) = default; 
	~PspMatrix() = default; 
	
	template< typename TripletIt, typename = decltype(*std::declval<TripletIt&>(), void(), ++std::declval<TripletIt&>(), void()) >
	void initialize(TripletIt b, const TripletIt e){
		using triplet_t = tuple< size_t, size_t, T >;
		using iter_t = typename std::iterator_traits< TripletIt >::value_type;
		static_assert(std::is_same< triplet_t, iter_t >::value, "Invalid triplet specification.");
		T val; 
		for (size_t i, j; b != e; ++b){
			std::tie(i,j,val) = *b;
			this->insert(i,j,val);
		}
	}
	
	// Tests  whether value_type is 0
	static bool equals_zero(value_type val) {
		if constexpr (std::is_integral_v< value_type >){
			return val == 0; 
		} else {
			return std::abs(val) <= std::numeric_limits< value_type >::epsilon();
		}
	}
	
	// Alternative triple-construction method 
	template< typename TripletIt >
	PspMatrix(TripletIt b, const TripletIt e, const size_t m, const size_t n) : PspMatrix(m,n) {
		initialize(b, e);
	}
	
	// Non-iterator construction 
	void construct(vector< size_t > i, vector< size_t > j, vector< T > x){
		// if (i.size() != j.size() || j.size() != x.size()){ Rcpp::stop("Invalid input"); }
		int max_row = *std::max_element(i.begin(), i.end());
		int max_col = *std::max_element(j.begin(), j.end());
		// if (max_row >= m->n_rows() || max_col >= m->n_cols()){ Rcpp::stop("Invalid input"); }
		auto nz_entries = vector< std::tuple< size_t, size_t, T > >();
		for (size_t pos = 0; pos < x.size(); ++pos){
			nz_entries.push_back(std::make_tuple(i[pos], j[pos], x[pos]));
		}
		this->initialize(nz_entries.begin(), nz_entries.end());
	}
	
	// Apply an lambda function to each non-zero element of the matrix
	template< typename Lambda >
	void apply(const Lambda& f){
		using output_t =  typename std::result_of< Lambda(size_t, size_t, T) >::type;
		static constexpr bool is_void = std::is_same< output_t, void >::value;
		static constexpr bool is_T = std::is_same< output_t, T >::value;
		static_assert(is_void || is_T, "Result type must be void or T.");
		
		for (size_t j = 0; j < n_cols(); ++j){
			if (columns[j]){
				vector< entry_t >& v = *columns[j];
				for (auto& e: v){
					if constexpr (is_void){
						f(otc[e.first], j, e.second);
					} else {
						T ne = f(otc[e.first], j, e.second);
						if (ne != e.second){ e.second = ne; }
					}
				}
			}
		}
	}
	
	// Applies given 'f' to non-zero entries of column 'c'. 
	// Calls 'f' with signature f(< row index >, < non-zero entry >)
	template < typename Lambda >
	void column(size_t c, Lambda&& f){
		if (c >= n_cols()){ return; }
		vector< entry_t >& v = *columns[c]; 
		for (entry_t& e: v){
			f(otc[e.first], e.second);
		}
	}
	
	// struct col_iterator {
	// 	
	// };
	
	
	template < typename OutputIter >
	void write_column(size_t c, OutputIter out){
		if (c >= n_cols()){ return; }
		if (!columns[c]){ columns[c] = std::make_unique< vector< entry_t > >(); }; 
		vector< entry_t >& col = *columns[c];
		std::transform(col.begin(), col.end(), out, [&](const entry_t& el){
			return(std::make_pair(otc[el.first], el.second));
		});
	}
	
	template < typename Iter >
	void assign_column(size_t c, Iter b, const Iter e){
		static_assert(std::is_same_v<typename std::iterator_traits< Iter >::value_type, entry_t >, "Incorrect value_type for iterator");
		if (c >= n_cols()){ throw std::invalid_argument("Invalid column index"); }
		
		// Re-map the current indices to the original indices to maintain consistency 
		std::transform(b, e, b, [&](entry_t& el){ return make_pair(cto[el.first], el.second); });
		if (!bool(columns.at(c))){ 
			columns.at(c) = std::make_unique< vector< entry_t > >(b,e); 
		} else {
			columns[c]->clear();
			columns[c]->assign(b,e);
		}
	}
	
	// Returns whether a column is completely zero
	bool column_empty(size_t j) const {
		if (j >= n_cols()){ 
			std::stringstream ss;
			ss << "Invalid column index " << j << " chosen.";
			std::string bad_msg = ss.str();
			throw std::invalid_argument(bad_msg); 
		}
		if (!bool(columns[j]) || columns[j]->size() == 0){ return(true); }
		bool col_empty = std::all_of(columns[j]->begin(), columns[j]->end(), [](auto entr){
			return(equals_zero(entr.second));
		});
		return(col_empty);
	}
	
	// Applies given 'f' to non-zero entries of row 'r'. 
	// Calls 'f' with signature f(< column index >, < non-zero entry >)
	template < typename Lambda >
	void row(size_t r, Lambda&& f){
		if (r >= n_rows()){ return; }
		for (size_t j = 0; j < n_cols(); ++j){
			vector< entry_t >& v = *columns[j];
			// TODO: fix this
			// Search for the a non-zero entry in row 'r' in O(log(n)) time 
			auto el = std::lower_bound(std::begin(v), std::end(v), r, [this](entry_t& e, size_t index){
				return(otc[e.first] < index);
			});
			// If not found, continue to the next column 
			if (el == std::end(v) || otc[el->first] != r){
				continue; 
			} else {
				f(j, el->second);
			}
		}
	}
	
	// Removes entries below a certain threshold from the matrix
	// TODO: Fix this
	void clean(const T threshold){
		auto to_remove = vector< pair< size_t, size_t > > (); 
		apply([threshold, &to_remove](size_t i, size_t j, T val){
			if (val <= threshold){
				to_remove.push_back(std::make_pair(i,j));
			}
		});
		for (auto p: to_remove){ remove(p.first, p.second); }
	}
	
	// Removes the (i,j)th entry of the matrix, if it exists
	// TODO: Fix this
	void remove(size_t i, size_t j){
		if (i >= n_rows() || j >= n_cols()){ return; }
		if (!bool(columns.at(j))){ return; }
		vector< entry_t >& entries = *columns.at(j);
		// auto el = std::lower_bound(std::begin(entries), std::end(entries), i, [this](entry_t& e, size_t index){
		// 	return(otc[e.first] < index);
		// });
		auto el = std::find_if(std::begin(entries), std::end(entries), [this, i](entry_t& e){
			return(otc[e.first] == i);
		});
		if (el == std::end(entries) || otc[el->first] != i){ 
			return; 
		} else {
			entries.erase(el);
			--nnz;
		}
	}
	
	// ( i = i + j ) 
	void add_cols(size_t i, size_t j){
		if (i >= n_cols() || j >= n_cols()){ throw std::invalid_argument("Invalid"); }
		if (!bool(columns.at(i))){ columns.at(i) = std::make_unique< vector< entry_t > >(); }
		if (!bool(columns.at(j))){ columns.at(j) = std::make_unique< vector< entry_t > >(); }
		
		auto first1 = columns[i]->begin(), last1 = columns[i]->end();
		auto first2 = columns[j]->begin(), last2 = columns[j]->end();
		vector< entry_t > to_add; 
		while(true){
			if (first1 == last1){ 
				nnz += std::distance(first2, last2);
				std::copy(first2, last2, std::back_inserter(*columns[i])); 
				break;
			}
			if (first2 == last2){ break; }
			if ((*first1).first == (*first2).first){
				(*first1).second = add((*first1).second, (*first2).second);
				++first1; ++first2;
			} else if ((*first1).first < (*first2).first){
				++first1;
			} else if ((*first2).first < (*first1).first){
				to_add.push_back((*first2));
				++first2;
			} else {
				throw std::logic_error("Invalid case");
			}
			Rcpp::checkUserInterrupt();
		}
		// Add entries retroactively that would've invalidated iterators above
		for (auto e: to_add){
			insert(otc[e.first], i, e.second);
		}
	}	
	
	// Performs the operation: row(i) = row(i) + row(j)  
	void add_rows(size_t i, size_t j){
		if (i >= n_rows() || j >= n_rows()){ throw std::invalid_argument("Invalid"); }
				
		// Collect the rows 
		auto row_i = vector< pair< size_t, T > >();
		auto row_j = vector< pair< size_t, T > >();
		row(i, [&row_i](auto c, auto e){ row_i.push_back(make_pair(c, e)); });
		row(j, [&row_j](auto c, auto e){ row_j.push_back(make_pair(c, e)); });
		
		auto first1 = row_i.begin(), last1 = row_i.end(); 
		auto first2 = row_j.begin(), last2 = row_j.end(); 
		
		auto to_add = vector< tuple< size_t, size_t, T > >(); 
		while(true){
			if (first1 == last1){
				nnz += std::distance(first2, last2);
				for (; first2 != last2; ++first2){
					to_add.push_back(make_tuple(i, (*first2).first, (*first2).second));
				}
				break;
			}
			if (first2 == last2){ break; }
			if ((*first1).first == (*first2).first){
				(*first1).second = add((*first1).second, (*first2).second);
				++first1;
				++first2;
			} else if ((*first1).first < (*first2).first){
				++first1;
			} else if ((*first2).first < (*first1).first){
				to_add.push_back(make_tuple(i, (*first2).first, (*first2).second));
				++first2;
			} else {
				throw std::logic_error("Invalid case");
			}
			Rcpp::checkUserInterrupt();
		}
		
		for (auto te: to_add){
			insert(otc[ std::get<0>(te) ], std::get<1>(te), std::get<2>(te));
		}
		
	}
	
	// Returns an optional containing the lowest-nonzero entry of column j
	optional< pair< size_t, value_type > > lowest_nonzero(const size_t j) const {
		// if column is empty, return nullopt
		if (column_empty(j)){ return(std::nullopt); }
		
		// Get max non-zero entry
		const auto& c = *columns[j];
		auto me = std::max_element(c.begin(), c.end(), [this](auto& e1, auto& e2){
			if (e1.second == 0){ return(true); }
			if (e2.second == 0){ return(false); }
			return otc[e1.first] < otc[e2.first];
		});
		
		// if max element is a (stored) zero, return nullopt
		if ((*me).second == 0){ return(std::nullopt); }
		
		return(std::make_optional(std::make_pair(otc[(*me).first], (*me).second)));
	}
	
	void swap_rows(size_t i, size_t j){
		if (i >= n_rows() || j >= n_cols()){ return; }
		const auto i_idx = cto[i], j_idx = cto[j];
		std::swap(cto[i], cto[j]);
		std::swap(otc[i_idx], otc[j_idx]);
	}
	void swap_cols(size_t i, size_t j){
		if (i >= n_cols() || j >= n_cols()){ return; }
		std::swap(columns[i], columns[j]);
	}
	void swap(size_t i, size_t j){
		swap_rows(i,j);
		swap_cols(i,j);
	}
	
	template< typename Iter >
	void permute_rows(Iter b, const Iter e){
		const size_t n = std::distance(b,e);
		if (n != n_rows()){ throw std::invalid_argument("Permutation must match number of rows."); }
		
		// Prepare permutation + inverse permutation
		vector< size_t > p(b,e), ip(n);
		for (size_t i = 0; i < n; ++i){ ip[p[i]] = i; } //   inverse permutation 
		
		// Map the row permutation to a new vector
		for (size_t i = 0; i < n; ++i){ 
			p[i] = cto[p[i]];    // compute permutation vector 
			ip[i] = cto[ip[i]];	 // compute inverse permutation vector
		}
		
		// Update the row correspondences
		std::move(p.begin(), p.end(), cto.begin());
		std::move(ip.begin(), ip.end(), otc.begin());
	}

	// Convenience method
	void permute_rows(const vector< size_t > p){
		this->permute_rows(p.begin(), p.end());
	}
	
	// Permutes the columns of the matrix according to the permutation given in [b,e)
	template< typename Iter >
	void permute_cols(Iter b, const Iter e){
		const size_t n = std::distance(b,e);
		if (n != n_cols()){ throw std::invalid_argument("Permutation must match number of columns."); }
		
		// Map the columns pointers to a new vector
		vector< entries_ptr > p(n);
		for (size_t i = 0; b != e; ++b, ++i){
			std::swap(p[i], columns[*b]);
		}	
		std::move(p.begin(), p.end(), columns.begin());
	}
	
	template< typename Iter >
	void permute(Iter b, const Iter e){
		if (n_rows() != n_cols()){ throw std::invalid_argument("Permute only possible for square matrices."); }
		permute_cols(b,e);
		permute_rows(b,e);
	}
	
	
	
	// void resize(const size_t n){
	// 	if (n == columns.size()){ return; }
	// 	columns.resize(n);
	// 	if (row_array.size() < n){
	// 		const size_t m = row_array.size();
	// 		row_array.resize(n);
	// 		std::iota(row_array.begin() + m, row_array.end(), m);
	// 	} else {
	// 		row_array.resize(n);
	// 	}
	// }
	
	template < typename OutputStream >
	void print(OutputStream& os){
		std::stringstream ss;
		ss << "Sparse matrix (" << std::to_string(n_rows()) << "," << std::to_string(n_cols()) << ")";
		os << ss.str() << std::endl;
		
		for (size_t c = 0; c < n_cols(); ++c){
			vector< entry_t >& v = *columns[c];
			if (v.size() > 0){
				os << std::to_string(c) << "|";
				for (auto& e: v){
					os << "(" << std::to_string(e.first) << ":" << std::to_string(e.second) << ") "; 
				}
				os << std::endl;
			}
		}
	}
	
	
	// Element access
	T operator()(size_t i, size_t j) const {
		if (i >= n_rows() || j >= n_cols() || !bool(columns.at(j))){ return(static_cast< T >(0)); }
		vector< entry_t >& entries = *columns.at(j);
		auto o_idx = cto[i];
		auto el = std::lower_bound(std::begin(entries), std::end(entries), o_idx, [](entry_t& e, size_t index){
			return(e.first < index);
		});
		if (el == std::end(entries)){ return(static_cast< T >(0)); }
		return(el->first == o_idx ? el->second : 0);
	}
	
	
	
	// T& operator()(size_t i, size_t j) {
	// 	if (i >= size.first || j >= size.second){
	// 		//throw std::invalid_argument("Invalid");
	// 		Rcpp::stop("invalid");
	// 	}
	// 	vector< entry_t >& entries = *columns.at(j);
	// 	auto el = std::lower_bound(std::begin(entries), std::end(entries), i, [](entry_t& e, size_t index){
	// 		return(e.first < index);
	// 	});
	// 	if (el == std::end(entries) || el->first != i){
	// 		auto new_el = entries.insert(el, std::make_pair(i, 0));
	// 		return std::reference_wrapper(entries[std::distance(entries.begin(), new_el)].second);
	// 	}
	// 	return(el->second);
	// }
		
	// Sets an element at position (i,j) to value 'val'
	// TODO: Fix this need to remove lower_bound, otherwise redundant entries exist , or need to guarentee sorting
	void insert(size_t i, size_t j, T val) {
		if (i >= n_rows() || j >= n_cols()){ throw std::invalid_argument("Invalid"); }
		if (!bool(columns.at(j))){ columns.at(j) = std::make_unique< vector< entry_t > >(); }
		vector< entry_t >& entries = *columns.at(j);
		auto o_idx = cto[i]; // original index to search for
		auto el = std::lower_bound(std::begin(entries), std::end(entries), o_idx, [this](entry_t& e, size_t index){
			return(e.first < index);
		});
		if (el == std::end(entries) || el->first != o_idx){
			entries.insert(el, std::make_pair(o_idx, val));
			++nnz;
		} else {
			el->second = val;
		}
	}
	
	
	// Convert sparse matrix into a Matrix sparseMatrix type (dgCMatrix) for compatibility with R
	S4 to_sparseMatrix(){
		const int RTYPE = Rcpp::traits::r_sexptype_traits< T >::rtype;
		IntegerVector m_dim = IntegerVector::create( n_rows(), n_cols() );
		
		Vector< RTYPE > x;
		IntegerVector row_i, col_j;
		apply([&x, &row_i, &col_j](size_t i, size_t j, T val){
			row_i.push_back(i);
			col_j.push_back(j);
			x.push_back(val);
		});
		
		std::string class_name = "dgCMatrix";
		Rcpp::Environment Matrix_pkg("package:Matrix"); 
		Rcpp::Function sm_constructor = Matrix_pkg["sparseMatrix"];    
		
		S4 res = sm_constructor(
			Rcpp::_["i"] = row_i,
			Rcpp::_["j"] = col_j,
			Rcpp::_["x"] = x,
			Rcpp::_["dims"] = m_dim, 
			Rcpp::_["index1"] = false
		);
		return(res);
	}
	S4 submatrix(const size_t i1, const size_t i2, const size_t j1, const size_t j2){
		if (j2 >= n_cols() || i2 >= n_rows() || i1 > i2 || j1 > j2){ throw std::invalid_argument("subscript out of bounds"); }
		const int RTYPE = Rcpp::traits::r_sexptype_traits< T >::rtype;
		IntegerVector m_dim = IntegerVector::create( i2-i1+1, j2-j1+1 );
		Vector< RTYPE > x;
		IntegerVector row_i, col_j;
		apply([&](size_t i, size_t j, T val){
			if (i >= i1 && i <= i2 && j >= j1 && j <= j2){
				row_i.push_back(i-i1);
				col_j.push_back(j-j1);
				x.push_back(val);	
			}
		});
		std::string class_name = "dgCMatrix";
		Rcpp::Environment Matrix_pkg("package:Matrix"); 
		Rcpp::Function sm_constructor = Matrix_pkg["sparseMatrix"];    
		S4 res = sm_constructor(
			Rcpp::_["i"] = row_i,
			Rcpp::_["j"] = col_j,
			Rcpp::_["x"] = x,
			Rcpp::_["dims"] = m_dim, 
			Rcpp::_["index1"] = false
		);
		return(res);
	}
	// S4 submatrix(std::span< size_t > i, std::span< size_t > j){
	// 	const int RTYPE = Rcpp::traits::r_sexptype_traits< T >::rtype;
	// 	
	// 	std::unique
	// 	
	// 	IntegerVector m_dim = IntegerVector::create( i2-i1+1, j2-j1+1 );
	// 	Vector< RTYPE > x;
	// 	IntegerVector row_i, col_j;
	// 	apply([&](size_t i, size_t j, T val){
	// 		if (i >= i1 && i <= i2 && j >= j1 && j <= j2){
	// 			row_i.push_back(i-i1);
	// 			col_j.push_back(j-j1);
	// 			x.push_back(val);	
	// 		}
	// 	});
	// 	std::string class_name = "dgCMatrix";
	// 	Rcpp::Environment Matrix_pkg("package:Matrix"); 
	// 	Rcpp::Function sm_constructor = Matrix_pkg["sparseMatrix"];    
	// 	S4 res = sm_constructor(
	// 		Rcpp::_["i"] = row_i,
	// 		Rcpp::_["j"] = col_j,
	// 		Rcpp::_["x"] = x,
	// 		Rcpp::_["dims"] = m_dim, 
	// 		Rcpp::_["index1"] = false
	// 	);
	// 	return(res);
	// }
	

};


// void transpose_schedule_R(SEXP Rp, SEXP Vp, IntegerVector S, Function f){
// 	XPtr< PspBoolMatrix > rp(Rp);
// 	XPtr< PspBoolMatrix > vp(Vp);
// 	PspBoolMatrix& R = *rp;
// 	PspBoolMatrix& V = *vp;
// 	vector< size_t > schedule(S.begin(), S.end());
// 	transpose_schedule_full(R, V, schedule, [&f](PspBoolMatrix& new_R, PspBoolMatrix& new_V, size_t s){
// 		f(new_R.to_sparseMatrix(), new_V.to_sparseMatrix(), s);
// 	});
// }


// S4 test_spmatrix(IntegerVector i, IntegerVector j, IntegerVector x, const size_t m, const size_t n){
// 	if (i.size() != j.size() || j.size() != x.size()){ Rcpp::stop("Invalid input"); }
// 	auto nz_entries = vector< std::tuple< size_t, size_t, bool > >();
// 	for (size_t pos = 0; pos < x.size(); ++pos){
// 		nz_entries.push_back(std::make_tuple(i[pos], j[pos], bool(x[pos])));
// 	}
// 	PspMatrix< bool, Mod2 > M(nz_entries.begin(), nz_entries.end(), size_t(m), size_t(n));
// 	// int x = M(0,0);
// 	// Rcout << "x: " << x << std::endl;
// 	// 
// 	// M.set(4,5,5);
// 	// M.set(2,4,2);
// 	// M.set(2,3,3);
// 	// M.set(8,9,1);
// 	// M.set(8,9,2);
// 	// // Rcout << "x: " << int(M(4,5)) << std::endl;
// 	// 
// 	// M.print(Rcout);
// 	// 
// 	// // using PspMatrix::entry_t;
// 	// M.apply([](size_t i, size_t j, int val) -> int { return val + 5; });
// 	// M.print(Rcout);
// 	// 
// 	// M.swap_rows(1,2);
// 	// M.swap_rows(2,4);
// 
// 	S4 output = M.to_sparseMatrix();;
// 	return(output);
// }
		// std::transform(b, e, b, [&](entry_t& el){ el.first = cto[el.first]; });
				// vector< entry_t > new_col; 
		// std::transform(b, e, std::back_inserter(new_col), [&](entry_t& el){ 
		// 	cto[el.first]; 
		// });
				// const auto i_idx = otc[i], j_idx = otc[j];
		// std::swap(otc[cto[i_idx]], otc[cto[j_idx]]);
		// std::swap(cto[i_idx], cto[j_idx]);
		// std::swap(cto[otc[i_idx]], cto[otc[j_idx]]);

