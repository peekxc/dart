
#include "PspMatrix.h" // Permutable Sparse Matrix template
#include "reduction.h" // reduction templates

// Define the field operation and specialize
struct Mod2 { 
	constexpr bool operator()(bool b1, bool b2){ return(b1 ^ b2); }
};
typedef PspMatrix< bool, Mod2 > PspBoolMatrix;

// TODO: switch to CRTP-facade like technique like https://vector-of-bool.github.io/2020/06/13/cpp20-iter-facade.html
// or use a type-erasure / duck-typing strategy 
struct BoolMatrix {
	using F = PspBoolMatrix::value_type;
	using value_type = F; 
	using entry_t = PspBoolMatrix::entry_t;
	
	PspBoolMatrix& m;
	BoolMatrix(PspBoolMatrix& pm) : m(pm){ };
	
	// { a.low(size_t(0)) } -> std::same_as< optional< pair< size_t, F > > >;
	auto low(size_t j) { return m.lowest_nonzero(j); }

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
	auto add_scaled_col(size_t i, size_t j, F lambda){
		m.add_cols(j,i);
	}
	
	// Zero's out column j
	auto clear_column(size_t j){
		if (m.columns.at(j)){
			m.columns[j]->clear();
		}
	}
	
	// Use column i to cancel lowest entry of j
	// { a.cancel_lowest(size_t(0), size_t(0)) } -> std::same_as< void >;
	auto cancel_lowest(size_t i, size_t j){
		// auto low_i = low(i), low_j = low(j);
		// auto lambda = low_j->second / low_i->second;
		m.add_cols(j, i);
	}
	
	// TODO: add cancel_lowest(size_t i, size_t j, low_entry& s_low, low_entry& t_low){
	
	// { a.dim() } -> std::same_as< pair< size_t, size_t > >;
	auto dim() -> pair< size_t, size_t > {
		return std::make_pair(m.size[0], m.size[1]);
	} 
	
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
	template <typename ... Args>
	void permute_cols(Args&& ... args){ m.permute_cols(std::forward<Args>(args)...); }
	
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
	void add_cols(Args&& ... args){ m.add_cols(std::forward<Args>(args)...); }
	
	// { a(size_t(0), size_t(0)) } -> std::same_as< F >; 
	template <typename ... Args>
	F operator()(Args&& ... args){ return m(std::forward<Args>(args)...); }
	
	template <typename ... Args>
	void assign_column(Args&& ... args){
		m.assign_column(std::forward<Args>(args)...);
	}
	
};

// [[Rcpp::export]]
void reduce_pspbool(SEXP Ds, SEXP Vs) {
	Rcpp::XPtr< PspBoolMatrix > D_ptr(Ds);
	Rcpp::XPtr< PspBoolMatrix > V_ptr(Vs);
	BoolMatrix D = BoolMatrix(*D_ptr); 
	BoolMatrix V = BoolMatrix(*V_ptr); 
	auto indices = vector< size_t >(D.dim().second);
	std::iota(indices.begin(), indices.end(), 0);
	pHcol(D, V, indices.begin(), indices.end());
}

// [[Rcpp::export]]
void reduce_local_pspbool(SEXP D1s, SEXP V1s, SEXP D2s, SEXP V2s, bool clearing) {
	Rcpp::XPtr< PspBoolMatrix > D1_ptr(D1s);
	Rcpp::XPtr< PspBoolMatrix > V1_ptr(V1s);
	Rcpp::XPtr< PspBoolMatrix > D2_ptr(D2s);
	Rcpp::XPtr< PspBoolMatrix > V2_ptr(V2s);
	BoolMatrix D1 = BoolMatrix(*D1_ptr); 
	BoolMatrix V1 = BoolMatrix(*V1_ptr); 
	BoolMatrix D2 = BoolMatrix(*D2_ptr); 
	BoolMatrix V2 = BoolMatrix(*V2_ptr); 
	auto d1_ind = vector< size_t >(D1.dim().second);
	auto d2_ind = vector< size_t >(D2.dim().second);
	std::iota(d1_ind.begin(), d1_ind.end(), 0);
	std::iota(d2_ind.begin(), d2_ind.end(), 0);
	if (clearing){
		pHcol_local< true >(D1, V1, D2, V2, d1_ind.begin(), d1_ind.end(), d2_ind.begin(), d2_ind.end());
	} else {
		pHcol_local< false >(D1, V1, D2, V2, d1_ind.begin(), d1_ind.end(), d2_ind.begin(), d2_ind.end());
	}
}

//[[Rcpp::export]]
void simulate_vineyard_cpp(SEXP Rs, SEXP Vs, IntegerVector schedule, Nullable< Function > f = R_NilValue){
	Rcpp::XPtr< PspBoolMatrix > R_ptr(Rs);
	Rcpp::XPtr< PspBoolMatrix > V_ptr(Vs);
	BoolMatrix R = BoolMatrix(*R_ptr); 
	BoolMatrix V = BoolMatrix(*V_ptr); 
	if (f.isNotNull()){
		Function fun(f);
		transpose_schedule_full(R, V, schedule.begin(), schedule.end(), [&](size_t status){
			fun(status);
		});
	} else {
		const auto do_nothing = [](size_t status){ return; };
		transpose_schedule_full(R, V, schedule.begin(), schedule.end(), do_nothing);
	}
}

//[[Rcpp::export]]
void move_schedule_local(SEXP r1, SEXP v1, SEXP r2, SEXP v2, IntegerVector schedule, Nullable< Function > f = R_NilValue){
	BoolMatrix R1 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(r1)); 
	BoolMatrix V1 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(v1)); 
	BoolMatrix R2 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(r2)); 
	BoolMatrix V2 = BoolMatrix(*Rcpp::XPtr< PspBoolMatrix >(v2)); 
	if (f.isNotNull()){
		Function fun(f);
		move_schedule_local(R1, V1, R2, V2, schedule.begin(), schedule.end(), [&](){
			// fun(status);
		});
	} else {
		const auto do_nothing = [](size_t status){ return; };
		move_schedule_local(R1, V1, R2, V2, schedule.begin(), schedule.end(), do_nothing);
	}
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
  	// .method("permute_cols", &PspIntMatrix::permute_cols)
    .method("add_rows", &PspBoolMatrix::add_rows)
    .method("low_entry", low_entry)
    .method("at", at)
    .method("submatrix", &PspBoolMatrix::submatrix)
	;
}

