


/*** R
# phtools:::test_matrix()

set.seed(1234)
rips <- pbgrad::r_geometric_complex(4, 0.25, filtered = TRUE)
D <- as.matrix(pbgrad::boundary_matrix(rips))
V <- diag(ncol(D))
wut <- phtools:::test_reduce(D, V)

all(((D %% 2) %*% (wut$V %% 2) %% 2) == (wut$R %% 2))


set.seed(1234)
rips <- pbgrad::r_geometric_complex(16, 0.25, dim = 2, filtered = TRUE)
D <- as.matrix(pbgrad::boundary_matrix(rips))
V <- diag(ncol(D))

idx <- which(D != 0, arr.ind = TRUE)
Dm <- phtools::psp_matrix(dim(D), i = idx[,1]-1L, j = idx[,2]-1L, x = D[idx] %% 2)
idx <- which(V != 0, arr.ind = TRUE)
Vm <- phtools::psp_matrix(dim(V), i = idx[,1]-1L, j = idx[,2]-1L, x = V[idx] %% 2)

phtools:::test_reduce2(Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr())
Dm$as.Matrix()
Vm$as.Matrix()

pbgrad::is_reduced(Dm$as.Matrix())
Matrix::isTriangular(Vm$as.Matrix())
all(((D %*% Vm$as.Matrix()) %% 2) == Dm$as.Matrix())


*/


/*** R
set.seed(1234)
# rips <- pbgrad::r_geometric_complex(4, 0.25, dim = 2, filtered = TRUE)
# D <- as.matrix(pbgrad::boundary_matrix(rips))
# V <- diag(ncol(D))

D <- pbgrad::boundary_matrix(simplextree::simplex_tree(combn(3,2)))
V <- diag(ncol(D))

idx <- which(D != 0, arr.ind = TRUE)
Dm <- phtools::psp_matrix(dim(D), i = idx[,1]-1L, j = idx[,2]-1L, x = D[idx] %% 2)
idx <- which(V != 0, arr.ind = TRUE)
Vm <- phtools::psp_matrix(dim(V), i = idx[,1]-1L, j = idx[,2]-1L, x = V[idx] %% 2)

phtools:::test_reduce2(Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr())
Dm$as.Matrix()
Vm$as.Matrix()

phtools:::simulate_vineyard_cpp(
	Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr(),
	c(4,3,4,4), 
	function(s){ print(s) }
)
Dm$as.Matrix()
Vm$as.Matrix()

## 
set.seed(1234)
rips <- pbgrad::r_geometric_complex(16, 0.30, dim = 2, filtered = TRUE)

D <- pbgrad::boundary_matrix(rips)
V <- diag(ncol(D))

idx <- which(D != 0, arr.ind = TRUE)
Dm <- phtools::psp_matrix(dim(D), i = idx[,1]-1L, j = idx[,2]-1L, x = D[idx] %% 2)
idx <- which(V != 0, arr.ind = TRUE)
Vm <- phtools::psp_matrix(dim(V), i = idx[,1]-1L, j = idx[,2]-1L, x = V[idx] %% 2)

phtools:::test_reduce2(Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr())
Dm$as.Matrix()
Vm$as.Matrix()

s0 <- unlist(simplex_to_str(rips$simplices))
s1 <- phtools::sample_filtration(s0)
S <- apply(schedule(list(s0,s1)), 2, min)-1L
phtools:::simulate_vineyard_cpp(
	Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr(), S, 
	function(s){ print(s) }
)
Dm$as.Matrix()
Vm$as.Matrix()

pbgrad::is_reduced(Dm$as.Matrix())
*/

// restore_right(Matrix& R, Matrix& V, Iter bi, const Iter ei, vector< typename Matrix::entry_t >& dr, vector< typename Matrix::entry_t >& dv){

// // From: https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095
// template<typename T>
// void apply_permutation(std::vector<T>& v, std::vector< size_t >& indices) {
// 	if (v.size() != indices.size()){ throw std::invalid_argument("indices size != v size"); }
// 	using std::swap; // to permit Koenig lookup
// 	for (size_t i = 0; i < indices.size(); i++) {
// 		auto current = i;
// 		while (i != indices[current]) {
// 			auto next = indices[current];
// 			swap(v[current], v[next]);
// 			indices[current] = current;
// 			current = next;
// 		}
// 		indices[current] = current;
// 	}
// }
// 
// 
// // [[Rcpp::export]]
// void move_right_cpp(SEXP Rs, SEXP Vs, IntegerVector schedule, Nullable< Function > f = R_NilValue){
// 	Rcpp::XPtr< PspBoolMatrix > R_ptr(Rs);
// 	Rcpp::XPtr< PspBoolMatrix > V_ptr(Vs);
// 	BoolMatrix R = BoolMatrix(*R_ptr); 
// 	BoolMatrix V = BoolMatrix(*V_ptr); 
// 	if (f.isNotNull()){
// 		Function fun(f);
// 		move_schedule_full(R, V, schedule.begin(), schedule.end(), [&](){
// 			fun(R.m.to_sparseMatrix(), V.m.to_sparseMatrix());
// 		});
// 	} else {
// 		const auto do_nothing = [](){ return; };
// 		move_schedule_full(R, V, schedule.begin(), schedule.end(), do_nothing);
// 	}
// }
// 
// // [[Rcpp::export]]
// void apply_move_permutaton(size_t i, size_t j){
// 	
// }
// 

/*** R
set.seed(1234)
D <- pbgrad::boundary_matrix(simplextree::simplex_tree(combn(3,2)))
V <- diag(ncol(D))

idx <- which(D != 0, arr.ind = TRUE)
Dm <- phtools::psp_matrix(dim(D), i = idx[,1]-1L, j = idx[,2]-1L, x = D[idx] %% 2)
idx <- which(V != 0, arr.ind = TRUE)
Vm <- phtools::psp_matrix(dim(V), i = idx[,1]-1L, j = idx[,2]-1L, x = V[idx] %% 2)

phtools:::test_reduce2(Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr())
Dm$as.Matrix()
Vm$as.Matrix()

phtools:::move_right_cpp(Dm$matrix$as_XPtr(), Vm$matrix$as_XPtr(), c(3,5), function(R,V){
	print("here")
})

Dm$as.Matrix()
Vm$as.Matrix()

*/


	// const size_t m2 = R2.dim().first;
	// for (size_t j; b2 != e2; ++b2){
	// 	j = *b2;
	// 	optional< entry_t > low_i, low_j;
	// 	while((low_j = R2.low(j)) && (low_i = R2.find_low(j, low_j->first))){
	// 		size_t i = low_i->first;
	// 		auto lambda = low_j->second / low_i->second;
	// 		R2.add_scaled_col(i, j, -lambda);
	// 		V2.add_scaled_col(i, j, -lambda);
	// 	}
	// 
	// 	
	// }