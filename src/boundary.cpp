#include <Rcpp.h>
using namespace Rcpp;

#include "implicit_filtration.h"


// Repeatedly appends whatever is in 'data' (count-1) times
// From: https://stackoverflow.com/questions/31389846/initializing-stdvector-with-a-repeating-pattern
template<typename Container>
void repeat_pattern(Container& data, std::size_t count) {
  auto pattern_size = data.size();
  if(count == 0 or pattern_size == 0) {
    return;
  }
  data.resize(pattern_size * count);
  const auto pbeg = data.begin();
  const auto pend = std::next(pbeg, pattern_size);
  auto it = std::next(data.begin(), pattern_size);
  for(std::size_t k = 1; k < count; ++k) {
    std::copy(pbeg, pend, it);
    std::advance(it, pattern_size);
  }
}

// [[Rcpp::export]]
S4 boundary_matrix_fi_full(SEXP filtration) {
	XPtr< ImplicitFiltration > fi_ptr(filtration);
	ImplicitFiltration& fi = *fi_ptr;
	
	// Prepare sparse matrix inputs 
	const auto nb = std::accumulate(fi.n_simplexes.begin(), fi.n_simplexes.end(), 0);
	vector< size_t > i, j;
	vector< int > x; 
	i.reserve(nb); j.reserve(nb); x.reserve(nb);
	
	// Apply boundary operator to every k-simplex
	auto x_out = std::back_inserter(x);
	auto i_out = std::back_inserter(i);
	auto j_out = std::back_inserter(j);
	for (size_t ii = 0, jj = 0; ii < fi.m; ++ii, ++jj){
		if (fi.dims[ii] > 0){
			const size_t k = fi.dims[ii]; 
			fi.boundary(ii, true, i_out);	
			std::fill_n(j_out, k+1, jj);
			for (size_t c = 0; c < (k+1); ++c){ *x_out++ = c % 2 == 0 ? 1 : -1; }
		}
	}

	// Create the sparse matrix
	Rcpp::Environment Matrix("package:Matrix"); 
  Rcpp::Function constructMatrix = Matrix["sparseMatrix"];    
	auto m_dim = IntegerVector::create(nb, nb);
	S4 m = constructMatrix(_["i"] = i, _["j"] = j, _["x"] = x, _["index1"] = false, _["dims"] = m_dim);
	
	return(m);
}

// [[Rcpp::export]]
S4 boundary_matrix_fi(SEXP filtration, const size_t k) {
	XPtr< ImplicitFiltration > fi_ptr(filtration);
	ImplicitFiltration& fi = *fi_ptr;
	
	// Prepare sparse matrix inputs 
	const auto nb = fi.n_simplexes.size() <= k ? 0 : fi.n_simplexes.at(k)*(k+1);
	vector< size_t > i, j;
	vector< int > x; 
	i.reserve(nb); j.reserve(nb); x.reserve(nb);
	
	// Apply boundary operator to every k-simplex
	if (k > 0){
		auto x_out = std::back_inserter(x);
		auto i_out = std::back_inserter(i);
		auto j_out = std::back_inserter(j);
		for (size_t ii = 0, jj = 0; ii < fi.m; ++ii){
			if (fi.dims[ii] == k){
				fi.boundary(ii, true, i_out);	
				std::fill_n(j_out, k+1, jj++);
			}
		}
		// The boundary operator data
		for (size_t c = 0; c < (k+1); ++c){ *x_out++ = c % 2 == 0 ? 1 : -1; }
		repeat_pattern(x, fi.n_simplexes.at(k));
	}
	
	// Remap row indices to match filtration order 
	if (k > 0){
		// f: (shortlex index) -> (filtration index)
		auto lex_to_full = inverse_permutation(fi.indexes.begin(), fi.indexes.end());	
		const auto offset = k == 1 ? 0 : fi.cum_ns.at(k-2);
		auto index_map = vector< size_t >();
		index_map.reserve(fi.n_simplexes.at(k-1));
		for (size_t i = 0; i < fi.m; ++i){
			if (fi.dims[i] == (k-1)){
				index_map.push_back(lex_to_full[i] - offset);
			}
		}
	}
	
	// Create the sparse matrix
	Rcpp::Environment Matrix("package:Matrix"); 
  Rcpp::Function constructMatrix = Matrix["sparseMatrix"];    
	auto m_dim = IntegerVector::create(fi.n_simplexes.at(k-1), fi.n_simplexes.at(k));
	S4 m = constructMatrix(_["i"] = i, _["j"] = j, _["x"] = x, _["index1"] = false, _["dims"] = m_dim);
	
	return(m);
}

// template< typename RandomIter, typename Size, typename OutputIt>
// void inverse_permutation_n(RandomIter b, Size n, OutputIt out){
//   for (Size i = 0; i != n; ++i){
//   	auto j = *(b+i);
//   	*(out+j) = i; 
//   }
// }


// Given vertex, edge, or simplex weights in dimension-lex order...
// List lex_filtration_order(SEXP stree, const NumericVector weights){
// 	Rcpp::XPtr<SimplexTree> stree_ptr(stree);
// 	SimplexTree& st = *stree_ptr; 
// 	
// 	// Ensure labels are contiguous and start with 0 or 1
// 	// vector< idx_t > v_labels = st.get_vertices();
// 	// if (v_labels.at(0) > 1){ throw std::invalid_argument("Vertex labels must be contiguous starting at 1 or 0."); }
// 	// bool is_1_based = v_labels.at(0) == 1; 
// 	// bool is_contiguous = true; 
// 	// for (size_t i = v_labels.at(0); i < v_labels.size() && is_contiguous; ++i){
// 	// 	is_contiguous &= v_labels[i - is_1_based] == i;
// 	// }
// 	// if (!is_contiguous){ throw std::invalid_argument("Vertex labels must be contiguous starting at 1 or 0."); }
// 	
// 	// TODO: intelligently detect 0:(nv-1) and 1:nv situations and handle directly 
// 	const vector< idx_t > v_labels = st.get_vertices(); // sorted already 
// 	const auto v_index = [&v_labels](idx_t label){
// 		return std::distance(v_labels.begin(), std::lower_bound(v_labels.begin(), v_labels.end(), label));
// 	};
// 	
// 	// Store filtration simplices as ( dim, weight, lexicographical rank )
// 	const size_t n = st.n_simplexes.at(0);
// 	const size_t m = std::accumulate(st.n_simplexes.begin(), st.n_simplexes.end(), 0);
// 	// using fs_t = typename std::tuple< size_t, double, size_t >;
// 	struct simplex_rank {
// 		unsigned int dim; 
// 		size_t rank;
// 	};
// 	auto F = vector< simplex_rank >(); 
// 	F.reserve(m);
// 	
// 	// Encode the filtration simplices 
// 	auto tr = st::level_order< true >(&st);
// 	
// 	// vector< vector< int > > S; 
// 	// for (size_t i = 0; i < filtration.size(); ++i){
// 	// 	auto si = dart::lex_unrank(filtration[i].rank, n, filtration[i].dim+1);
// 	// 	S.push_back(vector< int >(si.begin(), si.end()));
// 	// }
// 	// return(Rcpp::wrap(S));
// 
// 	// Create id vector to sort according to the refinement
//   std::vector< size_t > ids(m);
//   std::iota(ids.begin(), ids.end(), 0);
//   const auto machine_eq = [](double d1, double d2){
//   	return std::abs(d1 - d2) < std::numeric_limits< double >::epsilon();
//   };
// 
//   
//   
//   const auto max_e_weight = [&weights](auto b, auto e) -> double {
//   	for_each_combination(b, b + 2, e, [&weights](auto b, auto e){
//   		lex_rank()
//   	});
//   	return std::reduce(b, e, 0.0, [&weights](auto i, auto j){ return std::max(weights.at(i), weights.at(j)); });
//   };
//   
//   // Lower-star ordering
// 	if (weights.size() == st.n_simplexes.at(0)){
// 		// Computes max weight of simplex using lower-star
// 	  const auto max_v_weight = [&weights](auto b, auto e) -> double {
// 	  	return std::reduce(b, e, 0.0, [&weights](auto i, auto j){ return std::max(weights.at(i), weights.at(j)); });
// 	  };
//   	
//   	// Encode reindexed filtration  
// 		st::traverse(tr, [n, &v_index, &F](node_ptr p, size_t d, simplex_t simplex){
// 			std::transform(simplex.begin(), simplex.end(), simplex.begin(), v_index);
// 			F.push_back({ static_cast< unsigned int >(d - 1), dart::lex_rank(simplex, n) });
// 			return true; 
// 		});
// 		
// 		// Sort using reindexing/refinement map
// 		std::sort(ids.begin(), ids.end(), [m, n, &machine_eq, &max_v_weight, &F](size_t i, size_t j){
// 			auto si = dart::lex_unrank(F[i].rank, n, F[i].dim+1);
// 			auto sj = dart::lex_unrank(F[j].rank, n, F[j].dim+1);
// 			auto wi = max_v_weight(si.begin(), si.end());
// 			auto wj = max_v_weight(sj.begin(), sj.end());
// 	  	if (machine_eq(wi, wj)){ return(F[i].dim != F[j].dim ? F[i].dim < F[j].dim : F[i].rank < F[j].rank); }
// 	  	else if (wi < wj){ return true; }
// 			else { return false; }
// 	  });
// 	} else if (weights.size() == st.n_simplexes.at(1)){
// 		
// 		
// 		std::sort(ids.begin(), ids.end(), [m, n, &machine_eq, &max_weight, &F](size_t i, size_t j){
// 			auto si = dart::lex_unrank(F[i].rank, n, F[i].dim+1);
// 			auto sj = dart::lex_unrank(F[j].rank, n, F[j].dim+1);
// 			auto wi = max_weight(si.begin(), si.end());
// 			auto wj = max_weight(sj.begin(), sj.end());
// 	  	if (machine_eq(wi, wj)){ return(F[i].dim != F[j].dim ? F[i].dim < F[j].dim : F[i].rank < F[j].rank); }
// 	  	else if (wi < wj){ return true; }
// 			else { return false; }
// 	  });
// 	} else {
// 		if (weights.size() != m){ throw std::invalid_argument("Weights array not a valid size."); }
// 	}
// 		return(Rcpp::wrap(v_labels));
// 	// Lower-star ordering 
// 	// if (weights.size() == st.n_simplexes.at(0)){
// 		// auto tr = st::level_order< true >(&st);
// 		// st::traverse(tr, [](st::node_ptr p, size_t d, simplex_t simplex){
// 		// 	for_each_combination(simplex.begin(), simplex.begin() + (d-1), simplex.end(), [](auto b, auto e){
// 		// 
// 		// 	});
// 		// });
// 	// }
// 	
// 	
// 	// if (weights.size() != std::accumulate(st.n_simplexes.begin(), st.n_simplexes.begin(), 0)){
// 	// 	throw std::invalid_argument("weights vector size did not match size of complex.");
// 	// }
// }

// [[Rcpp::export]]
S4 boundary_matrix_st(SEXP stree, const size_t k) {
	Rcpp::XPtr<SimplexTree> stree_ptr(stree);
	SimplexTree& st = *stree_ptr; 
	
	// Map the faces to their combinadics 
	const size_t nv = st.n_simplexes[0];
	auto faces = st::k_simplices< true >(&st, st.root.get(), k-1);
	vector< size_t > face_idx;
	face_idx.reserve(st.n_simplexes.at(k));
	st::traverse(faces, [&face_idx, nv](node_ptr cn, idx_t depth, simplex_t x){
		face_idx.push_back(dart::lex_rank(std::span{x}, nv));
		return true; 
	});
	
	// Prepare containers to store non-zero indices
	auto simplices = st::k_simplices< true >(&st, st.root.get(), k);
	const size_t nb = st.n_simplexes.at(k)*(k+1); 
	vector< size_t > i, j, x;
	i.reserve(nb); j.reserve(nb); x.reserve(nb);
	
	// Apply the boundary operator to the k-simplices
	size_t jj = 0; 
	st::traverse(simplices, [&face_idx, &jj, &i, &j, &x, k, nv](node_ptr cn, idx_t depth, simplex_t s){
		size_t c = 0; 
		for_each_combination(s.begin(), s.begin()+k, s.end(), [&face_idx, &jj, &i, &j, &x, &c, k, nv](auto b, auto e){
			auto sx = simplex_t(b,e);
			auto face_index = dart::lex_rank(std::span{sx}, nv);
			auto it = std::lower_bound(face_idx.begin(), face_idx.end(), face_index);
			i.push_back(std::distance(face_idx.begin(), it));
			j.push_back(jj);
			x.push_back(c++ % 2 == 0 ? 1 : -1);
			return false; 
		});
		++jj;
		return true; 
	});
	
	// Create the sparse matrix
	Rcpp::Environment Matrix("package:Matrix"); 
  Rcpp::Function constructMatrix = Matrix["sparseMatrix"];    
	auto m_dim = IntegerVector::create(st.n_simplexes.at(k-1), st.n_simplexes.at(k));
	S4 m = constructMatrix(_["i"] = i, _["j"] = j, _["x"] = x, _["index1"] = false, _["dims"] = m_dim);
	
	// Return the lexicographically ordered boundary matrix
	return(m);
}


		
		
		
	// 	// Copy ranks for reference
	// 	const auto ranks_lex = ranks;
	// 	
	// 	// Start the id array for all the simplices - vertices 
	// 	auto ids = vector< size_t >(m-n);
	// 	std::iota(ids.begin(), ids.end(), n);
	// 	
	// 	auto s = reindex(simplex);
	// 	double diam = 0.0;
	// 	for_each_combination(s.begin(), s.begin() + 2, s.end(), [this, &diam](auto b, auto e){
	// 		auto it = std::lower_bound(ranks[1].begin(), ranks[1].end(), dart::lex_rank_2(*b, *(b+1), n));
	// 		diam = std::max(diam, grades.at(std::distance(it, ranks[1].end())));
	// 		return false;
	// 	});
	// 	
	// 	// Sort according to weight < dimension < lex order
	// 	std::sort(ids.begin(), ids.end(), [](auto i, auto j){
	// 		
	// 	});
	// 	
	// 	for (size_t ki = 2; ki < d; ++ki){
	// 		
	// 		// auto tr = st::k_simplices< true >(&st, ki);
	// 		// auto k_simplex_weights = vector< double >();
	// 		// k_simplex_weights.reserve(st.n_simplexes.at(1));
	// 		
	// 		// Calculate the diameters
	// 		// st::traverse(tr, [this, &k_simplex_weights](node_ptr p, size_t d, simplex_t simplex){
	// 		// 	auto s = reindex(simplex);
	// 		// 	double diam = 0.0; 
	// 		// 	for_each_combination(s.begin(), s.begin() + 2, s.end(), [this, &diam](auto b, auto e){
	// 		// 		auto it = std::lower_bound(ranks[1].begin(), ranks[1].end(), dart::lex_rank_2(*b, *(b+1), n));
	// 		// 		diam = std::max(diam, grades.at(std::distance(it, ranks[1].end())));
	// 		// 		return false; 
	// 		// 	});
	// 		// 	k_simplex_weights.push_back(diam);
	// 		// });
	// 		
	// 		// std::sort()
	// 	}
	// 	
	// 	// At this point: ranks are sorted lexicographically, weights are not sorted
	// 	std::sort(ranks[2].begin(), ranks[2].end(), [this](auto ri, auto rj){
	// 		auto si = dart::lex_unrank(ri, n, 2);
	// 		auto sj = dart::lex_unrank(ri, n, 2);
	// 	});
	// 	
	// 	// Sort dimension-specific
	// 	
	// 	// Build indices
	// 	auto d_vec = std::vector< double >(d-1);
	// 	for (size_t ki = 1; ki < d; ++ki){
	// 		// auto s = relabel(dart::lex_unrank(ranks[ki].at(0), n, ki+1));
	// 		weights s.begin(), s.end();
	// 		d_vec[ki-1] = ;
	// 	}
	// 	auto me = std::min_element(d_vec.begin(), d_vec.end());
	// 	
	// }
