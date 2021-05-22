


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


// template< PermutableMatrix Matrix, typename Iter >
// void restore_right(Matrix& R, Matrix& V, Iter bi, const Iter ei, vector< typename Matrix::entry_t >& dr, vector< typename Matrix::entry_t >& dv){
// 	using entry_t = typename Matrix::entry_t;
// 
// 	// Start with donor columns
// 	size_t donor_idx = *bi; 
// 	auto d_low_index = R.low_index(donor_idx);
// 	dr.clear(); dv.clear();
// 	R.write_column(donor_idx, std::back_inserter(dr));
// 	V.write_column(donor_idx, std::back_inserter(dv));
// 	if (std::distance(bi, ei) <= 1){ return; }
// 
// 	// Apply the donor concept
// 	std:advance(bi, 1);
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
// 		R.add_cols(k, donor_idx); V.add_cols(k, donor_idx);
// 		
// 		// Replace the donor columns 
// 		if ((d_low_index && new_low_index) && *d_low_index > *new_low_index){
// 			d_low_index = std::move(new_low_index);
// 			dr.swap(new_dr);
// 			dv.swap(new_dv);
// 		}
// 	}
// }
// 
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
// template< PermutableMatrix Matrix, typename Iter, typename Lambda >
// void move_schedule_full(Matrix& R, Matrix& V, Iter sb, const Iter se, Lambda f){
// 	using entry_t = typename Matrix::entry_t;
// 	using value_t = typename Matrix::value_type;
// 	const size_t nc = R.dim().first;
// 	const size_t nr = R.dim().second;
// 	if (nc == 0 || nr == 0){ return; }
// 	if (nc != nr){ throw std::invalid_argument("R must be square."); }
// 	if (std::distance(sb,se) < 2 || std::distance(sb,se) % 2 != 0){  throw std::invalid_argument("Pairs of indices must be passed."); }
// 
// 	for (size_t i, j; sb != se; sb += 2){
// 		i = *sb, j = *(sb+1);
// 		if (i == j){ continue; }
// 		if (i > j || i >= (nc-1) || j >= nc){  throw std::invalid_argument("Invalid pairs (i,j) passed."); }
// 
// 		// Collect indices I
// 		auto I = vector< size_t >();
// 		V.row(i, [&](auto col_idx, auto v){
// 			if ((col_idx >= i) & (col_idx <= j)){ I.push_back(col_idx); }
// 		});
// 
// 		// Collect indices J
// 		auto J = vector< size_t >();
// 		for (size_t c = 0; c < nc; ++c){
// 			auto low_idx = R.low_index(c);
// 			if (low_idx && (*low_idx) >= i && (low_idx) <= j && R(i,c) != 0){
// 				J.push_back(c);
// 			}
// 		}
// 		
// 		Rcout << "I: "; 
// 		for (auto ii: I){ Rcout << ii << ", "; }
// 		Rcout << std::endl;
// 		
// 		Rcout << "J: "; 
// 		for (auto ii: J){ Rcout << ii << ", "; }
// 		Rcout << std::endl;
// 		 
// 
// 		// Restore invariants
// 		auto dr1 = vector< entry_t >();
// 		auto dv1 = vector< entry_t >();
// 		restore_right(R, V, I.begin(), I.end(), dr1, dv1); // don't wrap in if condition
// 		
// 		// Restore invariants
// 		if (!J.empty()){
// 			auto dr2 = vector< entry_t >();
// 			auto dv2 = vector< entry_t >();
// 			restore_right(R, V, J.begin(), J.end(), dr2, dv2);
// 		}
// 	
// 
// 		// Apply permutations
// 		vector< size_t > p(nc);
// 		std::iota(p.begin(), p.end(), 0);
// 		std::rotate(p.begin()+i, p.begin()+i+1, p.begin()+j);
// 		
// 		for (auto pe: p){ Rcout << pe << " "; }
// 		Rcout << std::endl;
// 		
// 		R.permute(p.begin(), p.end());
// 		V.permute(p.begin(), p.end());
// 		
// 		Rcout << "Donor dr: ";
// 		for (auto pe: dr1){ Rcout << pe.second << " "; }
// 		Rcout << std::endl;
// 		
// 		Rcout << "Donor dv: ";
// 		for (auto pe: dv1){ Rcout << pe.second << " "; }
// 		Rcout << std::endl;
// 		
// 		// Perform donor replacement
// 		apply_permutation(dr1, p); 
// 		apply_permutation(dv1, p);
// 		
// 		Rcout << "Donor dr: ";
// 		for (auto pe: dr1){ Rcout << pe.second << " "; }
// 		Rcout << std::endl;
// 		
// 		Rcout << "Donor dv: ";
// 		for (auto pe: dv1){ Rcout << pe.second << " "; }
// 		Rcout << std::endl;
// 		
// 		R.assign_column(j, dr1.begin(), dr1.end());
// 		V.assign_column(j, dv1.begin(), dv1.end());
// 	}
// }
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