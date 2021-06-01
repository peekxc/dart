#ifndef IMPLICIT_FILTRATION_H
#define IMPLICIT_FILTRATION_H

// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"
#include "utility/combinations.h"
#include "utility/discrete.h"
#include "combinadic.h"
#include <tuple> 
#include <cstdint>
#include <iterator>
#include <chrono>

// For comparing very close values
template< typename T > [[nodiscard]]
constexpr bool within_machine_prec(T d1, T d2) noexcept {
	static_assert(std::is_floating_point< T >::value, "Must be floating point value.");
	return(std::abs(d1 - d2) < std::numeric_limits< T >::epsilon());
}


// Implicit Filtration class
// Stores a low-dimensional simplexwise filtration of simplices using a light encoding (combinatorial number system). 
// Filtration grades are stored minimally. 
// After construction, allows O(d log m) access to individual simplex < ranks / grades / weights >, 
// or alternatively O(m) access to stream the full filtration
struct ImplicitFiltration {
	const size_t m; 									 // number of total simplices 
	const size_t n; 									 // number of vertices
	const size_t d; 									 // dimension of complex + 1
	const vector< double > grades;		 // filtration grades; size depends on complex
	const vector< idx_t > labels;			 // vertex original labels
	vector< uint_fast64_t > ranks;     // simplex ranks 
	vector< size_t > indexes;          // simplex positions
	vector< uint_least8_t > dims;			 // simplex dimensions
	vector< size_t > n_simplexes; 		 // number of simplices of each dimension
	vector< size_t > cum_ns; 
		
	enum class filt_t { LOWER_STAR, FLAG, GENERIC };
	filt_t type; 
	
	// Return types for () and []
	enum class elem_t { RANK, GRADE, SIMPLEX, ALL };
	
	template< typename Iter >
	ImplicitFiltration(SimplexTree& st, Iter w_begin, const Iter w_end) 
		: m(std::accumulate(st.n_simplexes.begin(), st.n_simplexes.end(), 0)),
    	n(st.n_simplexes.at(0)),
    	d(st.dimension()+1),
    	grades(w_begin, w_end), 
    	labels(st.get_vertices()), 
      n_simplexes(st.n_simplexes.begin(), st.n_simplexes.end()) {
		
		std::partial_sum(n_simplexes.begin(), n_simplexes.begin()+d, std::back_inserter(cum_ns));
		// TODO: check choose(m, d) < 2^64 - 1 for encoding 
		
		// Check the type of filtration to construct
		if (grades.size() == n){
			type = filt_t::LOWER_STAR; 
		} else if (grades.size() == st.n_simplexes.at(1)){
			type = filt_t::FLAG; 
		} else if (grades.size() == m){
			type = filt_t::GENERIC; 
		} else {
			// Rprintf("Weights size: %d, n: %d, m: %d\n", grades.size(), n, m);
			throw std::invalid_argument("Filtration weights must match size of vertices, edges, or all simplices.");
		}
		// auto start1 = std::chrono::steady_clock::now();
    
		// Populate initial ranks + dims + indexes in shortlex order 
		auto tr = st::level_order< true >(&st);
		ranks.reserve(m); 
		dims.reserve(m);
		st::traverse(tr, [this](node_ptr p, size_t d, simplex_t simplex){
			ranks.push_back(dart::lex_rank(reindex(std::span{simplex}), n));
			dims.push_back(static_cast< uint_least8_t >(d-1));
			return true; 
		});
		indexes = vector< size_t >(m);
		std::iota(indexes.begin(), indexes.end(), 0);
	   
// 	  auto end1 = std::chrono::steady_clock::now();
//     std::chrono::duration<double> ms1 = (end1-start1);
//     Rcout << "initial traversal took: " << ms1.count()*1000 << " ms" << std::endl;
    
		// Construct the filtration 
		// Postcondition: indexes array is sorted in filtration order
		// auto start2 = std::chrono::steady_clock::now();
		// if (type != filt_t::FLAG){ throw std::invalid_argument("Filtration weights must match size of vertices, edges, or all simplices."); }
		switch(type){
			case filt_t::LOWER_STAR:{ build_lower_star(); break; }
			case filt_t::FLAG: { build_flag(); break; }
			default: { build_generic(); break; }
		}
// 		auto end2 = std::chrono::steady_clock::now();
// 		std::chrono::duration<double> ms2 = (end2-start2);
//     Rcout << "building complex took: " << ms2.count()*1000 << " ms" << std::endl;

		// Sort ranks and dims vectors according to filtration order
		{
			auto p = indexes;
			apply_permutation(ranks.begin(), ranks.end(), p.begin());
			p = indexes;
			apply_permutation(dims.begin(), dims.end(), p.begin());
		}
// 		auto end3 = std::chrono::steady_clock::now();
// 		std::chrono::duration<double> ms3 = (end3-start3);
//     Rcout << "Applying permutations took: " << ms3.count()*1000 << " ms" << std::endl;

		// Invert indexes to get O(log n) lookup using shortlex order
		// ranks[indexes[i]] := i'th simplex in shortlex order
		// auto start4 = std::chrono::steady_clock::now();
		{
			auto ip = vector< size_t >(m);
			std::iota(ip.begin(), ip.end(), 0);
			apply_reverse_permutation(ip.begin(), ip.end(), indexes.begin());
			indexes = std::move(ip);
		}
// 		auto end4 = std::chrono::steady_clock::now();
// 		std::chrono::duration<double> ms4 = (end4-start4);
//     Rcout << "Inverse permutations took: " << ms4.count()*1000 << " ms" << std::endl;
	}
	
	void build_lower_star(){
		
		// Given a simplex index 'i', compute its grade via its vertices maximum grades
		const auto grade = [this](size_t i){
			if (dims[i] == 0){ return(grades[i]); }
			auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
			return std::reduce(s.begin(), s.end(), 0.0, [this](auto i, auto j){ return std::max(grades.at(i), grades.at(j)); });
		};
	
		// Created the lexicographically-refined ordering
		std::sort(indexes.begin(), indexes.end(), [this, &grade](auto i, auto j){
			auto wi = grade(i), wj = grade(j);
			if (within_machine_prec(wi, wj)){ return(dims[i] == dims[j] ? ranks[i] < ranks[j] : dims[i] < dims[j]); }
			return(wi < wj);
		});
	}
	
	// Created the lexicographically-refined ordering for a generic filtration
	void build_generic(){
		std::sort(indexes.begin(), indexes.end(), [this](auto i, auto j){
			auto wi = grades[i], wj = grades[j];
			if (within_machine_prec(wi, wj)){ return(dims[i] == dims[j] ? ranks[i] < ranks[j] : dims[i] < dims[j]); }
			return(wi < wj);
		});
	}
	
	// Builds a filtration from a flag complex (from edge weights) 
	// Preconditions: 
	// 1. grades.size() == st.n_simplexes[1] 
	// 2. indexes, ranks, and dims vectors are populated in shortlex order
	void build_flag(){
		
		// Extract edge ranks since they define the flag complex 
		// Ranks will be in shortlex order, which is naturally sorted
		auto edge_ranks = vector< uint_fast64_t >(); 
		edge_ranks.reserve(n_simplexes.at(1));
		for (size_t i = 0; i < m; ++i){
			if (dims[i] == 1){ edge_ranks.push_back(ranks[i]); }
		}
		
		// Given an index 'i', obtains the diameter of the simplex at position 'i' in the shortlex order
		const auto diameter = [this, &edge_ranks](size_t i) -> double {
			if (dims[i] == 0){ return 0.0; }
			else if (dims[i] == 1){ 
				auto e_it = std::lower_bound(edge_ranks.begin(), edge_ranks.end(), ranks[i]);
				return(grades[std::distance(edge_ranks.begin(), e_it)]);
			} else {
				auto s = dart::lex_unrank(ranks[i], n, dims[i]+1); 
				double diam = 0.0;
				for_each_combination(s.begin(), s.begin() + 2, s.end(), [this, &diam, &edge_ranks](auto b, auto e){
					auto it = std::lower_bound(edge_ranks.begin(), edge_ranks.end(), dart::lex_rank_2(*b, *(b+1), n));
					diam = std::max(diam, grades[std::distance(edge_ranks.begin(), it)]);
					return false;
				});
				return(diam);
			}
		};
		
		// Cache diameters before sort
		auto diameters = vector< double >(); 
		diameters.reserve(m);
		for (size_t i = 0; i < m; ++i){ diameters.push_back(diameter(i)); }
			
		// Created the lexicographically-refined ordering
		std::sort(indexes.begin(), indexes.end(), [this, &diameters](auto i, auto j){
			auto wi = diameters[i], wj = diameters[j];
			if (within_machine_prec(wi, wj)){ return(dims[i] == dims[j] ? ranks[i] < ranks[j] : dims[i] < dims[j]); }
			return(wi < wj);
		});
	}
	
	// indexes := i -> j := ( short lex i ) -> filtration index j such that ranks[indexes[i]] < ranks[indexes[i+1]]
	double explicit_grade(const size_t i) const {
		if (i >= m){ throw std::invalid_argument("Index larger than filtration"); }
		if (type == filt_t::FLAG){
			if (dims[i] == 0){ return 0.0; }
			else if (dims[i] == 1){
				// Rcout << "Searching for: " << i << " with rank " << ranks[i] << std::endl;
				const auto lb = cum_ns.at(0), ub = cum_ns.at(1);
				const auto it = std::lower_bound(indexes.begin()+lb, indexes.begin()+ub, ranks.at(i), [this](auto j, auto r){
					return(ranks.at(j) < r);
				});
				if (it != indexes.begin()+ub){
					auto gr_idx = std::distance(indexes.begin()+lb, it);
					return(grades.at(gr_idx));
				} else {
					throw std::invalid_argument("Invalid index/dimension given. Edge rank not found.");
				}
			} else {
				auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
				const auto lb = cum_ns.at(0), ub = cum_ns.at(1);
				double diam = 0.0;
				const auto ib = indexes.begin()+lb;
				const auto ie = indexes.begin()+ub;
				for_each_combination(s.begin(), s.begin() + 2, s.end(), [this, &ie, &ib, &diam](auto b, auto e){
					auto it = std::lower_bound(ib, ie, dart::lex_rank_2(*b, *(b+1), n), [this](auto j, auto r){
						return(ranks[j] < r);
					});
					diam = std::max(diam, grades.at(std::distance(ib, it)));
					return false;
				});
				return(diam);
			}
		} else if (type == filt_t::LOWER_STAR){
			// Given a simplex index 'i', compute its grade via its vertices maximum grades
			if (dims[i] == 0){ return(grades[i]); }
			auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
			return std::reduce(s.begin(), s.end(), 0.0, [this](auto i, auto j){ return std::max(grades.at(i), grades.at(j)); });
		} else {
			return(grades[i]);
		}
		return(-1.0);
	}
	vector< size_t > explicit_simplex(const size_t i) const  {
		if (i >= m){ throw std::invalid_argument("Index larger than filtration"); }
		auto s = dart::lex_unrank(ranks[i], n, dims[i]+1);
		// std::span< size_t >(s.begin(), s.end())
		auto ss = relabel(std::span{s});
		// auto sx = relabel(std::span< size_t >(s.begin(), s.end()));
		return(vector< size_t >(ss.begin(), ss.end()));
	}
	
	// Applies the boundary operator to the simplex at position i in the filtration
	template< typename OutputIt >
	void boundary(const size_t i, bool relative, OutputIt out) const {
		if (i >= m){ throw std::invalid_argument("Invalid filtration index given."); }
		const auto d = dims[i];
		if (d == 0){ return; }
		
		// Get the k-simplex at position i in the filtration
		vector< dart::I > s = dart::lex_unrank(ranks[i], n, d+1);
		std::cout << "simplex: " << std::endl;
		for (auto label: s){
			std::cout << label << ",";
		}
		std::cout << std::endl; 
					
		// Range of (k-1)-faces in the filtration (contiguous in indexes order!)
		const auto ib = indexes.begin()+(d == 1 ? 0 : cum_ns.at(d-2));
		const auto ie = indexes.begin()+cum_ns.at(d-1);
		
		// Enumerate faces, outputing relative indices as needed 
		for_each_combination(s.begin(), s.begin() + d, s.end(), [this, &ie, &ib, &out, relative](auto b, auto e){
			
			auto face = vector< dart::I >(b, e);
			std::cout << "face: "; 
			for (auto label: face){
				std::cout << label << ",";
			}
			std::cout << std::endl; 
			
			auto face_rank = dart::lex_rank(std::span{face}, n);
			std::cout << "searching for rank: " << face_rank << std::endl; 
			auto face_it = std::lower_bound(ib, ie, face_rank, [this](auto j, auto r){
				return(ranks[j] < r);
			});
			// Boundary indices relative to dimension 
			if (face_it != ie){
				std::cout << "writing: " << std::distance(ib, face_it) << ";" << std::endl; 
				*out++ = relative ? std::distance(ib, face_it) : *face_it;
			}
			return false;
		});
	}
	
	template< typename T >
	void grading(std::span< T > G) const {
		if (G.size() != m){ throw std::invalid_argument("Output range must match filtration size."); }
		for (size_t i = 0; i < m; ++i){
			G[i] = explicit_grade(i);
		}
	}
	template< typename T >
	void k_grades(std::span< T > G, const size_t k) const {
		if (G.size() != n_simplexes.at(k)){ throw std::invalid_argument("Output range must match filtration size."); }
		for (size_t i = 0, cc = 0; i < m; ++i){
			if (dims[i] == k){
				G[cc++] = explicit_grade(i);
			}
		}
	}
	
	// Returns the k-simplices in filtration order (unranked)
	// template< typename OutputIt >
	// void k_simplices(const size_t k, OutputIt out) const {
	// 	if (k >= d){ return; }
	// 	dart::lex_unrank(ranks[k].begin(), ranks[k].end(), n, k+1, out);
	// } // todo: change to match below
	
	// T must accept a span< size_t > as a constructor
	template< typename T >
	void simplices(std::span< T > S) const {
		if (S.size() != m){ throw std::invalid_argument("Output range must match filtration size."); }
		for (size_t i = 0; i < m; ++i){
			auto sx = dart::lex_unrank(ranks[i], n, dims[i]+1);
			auto ss = relabel(std::span{sx});
			S[i] = T(ss.begin(), ss.end());
		}
	}
	
	// T must accept a span< size_t > as a constructor
	template< typename T >
	void k_simplices(std::span< T > S, const size_t k) const {
		if (S.size() != n_simplexes.at(k)){ throw std::invalid_argument("Output range must match filtration size."); }
		for (size_t i = 0, cc = 0; i < m; ++i){
			if (dims[i] == k){
				auto sx = dart::lex_unrank(ranks[i], n, k+1);
				auto ss = relabel(std::span{sx});
				S[cc++] = T(ss.begin(), ss.end());	
			}
		}
	}
	
	// Looks for a simplex in the filtration. If found, returns the filtration index
	// auto find(std::span< size_t > simplex) const noexcept -> std::optional< size_t > {
	// 	if (simplex.size() > d || simplex.size() == 0){ return(std::nullopt); }
	// 	const size_t ki = simplex.size()-1; 
	// 	auto it = std::lower_bound(ranks[ki].begin(), ranks[ki].end(), dart::lex_rank(simplex, n));
	// 	return(it != ranks[ki].end() ? std::make_optional(std::distance(ranks[ki].begin(), it)) : std::nullopt);
	// }
		
	// Returns the 0-based (contiguous) index of a vertex label ( O(log n) time )
	auto vertex_index(idx_t label) const noexcept -> size_t {
		return std::distance(labels.begin(), std::lower_bound(labels.begin(), labels.end(), label));
	};
	
	template< typename Size_t > [[nodiscard]]
	auto relabel(std::span< Size_t > s) const -> std::span< Size_t > {
		std::transform(s.begin(), s.end(), s.begin(), [this](auto i){ return(this->labels[i]); });
		return(s);
	}
	
	template< typename Size_t > [[nodiscard]]
	auto reindex(std::span< Size_t > s) const -> std::span< Size_t > {
		std::transform(s.begin(), s.end(), s.begin(), [this](auto l){ return(this->vertex_index(l)); });
		return(s);
	}
};

#endif
