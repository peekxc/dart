#include <Rcpp.h>
using namespace Rcpp;


#include "avl_tree.h"
#include "cutset.h"
#include <bit>
#include <bitset>
#include <iostream>

#include <optional>

// [[Rcpp::plugins(cpp20)]]                                        


// [[Rcpp::export]]
void bst(IntegerVector x, int type) {
	
	auto tree = AvlTree< int >(x.begin(), x.end());
	
	tree.print(type);

  return; 
}


using std::vector; 

void dfs_al(vector< vector<int> > &adj, int node, int h = 0) {
	
	
  // for (auto to : adj[node]) {
  //   if (!visited[to]) {
  //       dfs(adj, to, h + 1);
  //       euler.push_back(node);
  //   }
  // }
  
	auto rec_lambda = [](auto&& self, int n){
	    if(n < 2) return 1;               // return type is deduced here
	    else return n * self(/* args */); // this has no impact
	};
  
}


template< typename Iter >
auto lex_to_al(Iter b, const Iter e, const vector< int >& ids, const bool symmetric) -> vector< vector< int > > {
	using value_t = typename std::iterator_traits<Iter>::value_type;
	auto me = std::max_element(b, e);
	vector< vector< int > > al((*me) + 1);
	for (; b != e; b += 2){
		value_t s = *b, t = *(b + 1);
		al[s].push_back(ids[t]);
		if (symmetric){ al[t].push_back(ids[s]); }
	}
	return(al);
}

// Converts an edgelist to an adjacency list
// [[Rcpp::export]]
ListOf< IntegerVector > el_to_al_0(IntegerMatrix el, IntegerVector ids, const bool symmetric) {
	vector< int > ids_copy(ids.begin(), ids.end());
	return(wrap(lex_to_al(el.begin(), el.end(), ids_copy, symmetric)));
}


template<typename Functor>
struct fix_type {
    Functor functor;
    template<typename... Args>
    decltype(auto) operator()(Args&&... args) const&
    { return functor(functor, std::forward<Args>(args)...); }
};

template<typename Functor>
fix_type<typename std::decay<Functor>::type> fix(Functor&& functor)
{ return { std::forward<Functor>(functor) }; }

// [[Rcpp::export]]
IntegerVector euler_tour(ListOf< IntegerVector > al){
	const vector< vector<int> > adj = Rcpp::as< vector< vector< int > > >(al);
	
	int n = adj.size();
	vector< int > height(n), first(n), euler;
	euler.reserve(n * 2);
  vector< bool > visited(n, false);
  
//   auto factorial = fix([](auto&& self, int n) -> int
// 	{ return n < 2 ? 1 : n * self(self, n - 1); });
//   factorial(10);
  auto et_f = fix([&](auto&& et, int node, int h) -> void {
  	visited[node] = true;
		height[node] = h;
		first[node] = euler.size();
	  euler.push_back(node);
		for (auto to : adj[node]) {
		  if (!visited[to]) {
		      et(et, to, h + 1);
		      euler.push_back(node);
		  }
		}
  });
  et_f(0, 0); // execute the tour
  
  return(wrap(euler));
}

using std::pair;
using std::vector;

// Stores and (directed-edge)-based Euler tour that supports re-rooting,  
struct EulerTour {
	using edge_t = pair< size_t, size_t >;
	vector< edge_t > tour;  
	
	EulerTour(const size_t v) : tour(1, edge_t(v,v)){ }
	
	template< typename Iter >
	EulerTour(Iter b, Iter e) : tour(b, e){
		using value_t = typename std::iterator_traits< Iter >::value_type;
		// static_assert(std::is_same< value_t, T >::value, "Iterator must match the parameterized data type.");
		// tour = vector< edge_t >(b, e); // makes the spanning tree
	}
	
	// Concatenates an euler tour in the range (b,e) to the current tour 
	template< typename Iter >
	void merge(Iter b, Iter e){
		tour.insert(tour.end(), b, e);
	}
	
	// Re-roots the euler tour such that (u,u) is at the front
	void reroot(size_t u){
		auto vertex = std::find(tour.begin(), tour.end(), edge_t(u, u));
		if (vertex == tour.begin()){ return; }
		if (vertex != tour.end()){
			std::rotate(tour.begin(), vertex, tour.end());
		}
	}
	
	// Check if (u,v) is a tree edge (by def. of euler tour, (v,u) is also an edge)
	bool is_tree_edge(const size_t u, const size_t v) const {
		auto e = std::find(tour.begin(), tour.end(), edge_t(u, v));
		return(e != tour.end());
	}
	
	// Check if u and v exists in the spanning tree
	bool is_connected(const size_t u, const size_t v) const {
		auto v1 = std::find(tour.begin(), tour.end(), edge_t(u, u));
		if (v1 == tour.end()){ return(false); }
		auto v2 = std::find(tour.begin(), tour.end(), edge_t(v, v));
		if (v2 == tour.end()){ return(false); }
		return(true);
	}
	
	// Returns the set of vertices in the tour
	auto vertices() const -> vector< size_t > {
		vector< size_t > result; 
		for (auto p: tour){ if (p.first == p.second){ result.push_back(p.first); } }
		return(result);
	}
	
};

template < size_t w >
struct DynamicGraph {
	using edge_t = EulerTour::edge_t;
	using v_label = std::bitset< w >;
	// using edge_t = pair< label_t, label_t >;
	
	
                                      
	// const size_t n; 									// number of vertices (fixed)
	vector< v_label > components; 				// cached component index for each vertex
	vector< EulerTour > euler_tours;  // component trees, stored as euler tours
	// vector< edge_t > aux;							// auxiliary edges 
	vector< bool > aux; 									// auxiliary edges *not* in any given tree
	
	constexpr size_t n_vertices() const { return components.size(); }
	constexpr size_t component(const size_t u) const {
		return static_cast< size_t >(components[u].to_ulong());
	}
		
	// template< typename Iter >
	DynamicGraph(size_t nv) : components(nv), aux(nv*(nv-1)/2, false) {
		// assert(w == std::ceil(std::log2(nv)));
		// using value_t = typename std::iterator_traits< Iter >::value_type;
		// // static_assert(std::is_same< value_t, T >::value, "Iterator must match the parameterized data type.");
		// tour = vector< edge_t >(b, e); // makes the spanning tree
		std::iota(components.begin(), components.end(), 0);
		for (auto label: components){
			euler_tours.push_back(EulerTour(label.to_ulong()));
		}
	}
	
	// Checks if an edge exists in the graph
	bool edge_exists(const size_t u, const size_t v) const {
		if (u >= n_vertices() || v >= n_vertices()){ return false ; } 			  // invalid edge
		if (component(u) != component(v)){ return false; }								  // disjoint components => no edge exists
		if (aux[rank_lex_2(u, v, n_vertices())]){ return true; }												// check if edge is in auxillary
		return(euler_tours[component(u)].is_tree_edge(u, v));
		// if (euler_tours[components[u]].is_connected(u, v)){ return true; }	  // see if it's part of the tree
		// return(std::find(aux.begin(), aux.end(), edge_t(u, v)) != aux.end()); // 
	}
	
	// Returns if there is a path from u to v 
	bool is_connected(const size_t u, const size_t v) const {
		if (edge_exists(u,v)){ return(true); }
		return(euler_tours[component(u)].is_connected(u, v));
	}
	
	// Creates a new undirected edge (u,v) using the euler tour representations.
	// The connected component labels are updated. Returns whether a new edge was created or not.  
	bool link(size_t u, size_t v){
		if (u > v){ std::swap(u, v); }
		if (u == v || edge_exists(u, v) || (u >= n_vertices() || v >= n_vertices())){ return false; }
		if (component(u) == component(v)){
			if (euler_tours[component(u)].is_tree_edge(u, v)){ return(false); }
			aux[rank_lex_2(u, v, n_vertices())] = true; 
			return(true);
		} else {
			const auto u_cc = component(u), v_cc = component(v);
			EulerTour& et1 = euler_tours[u_cc];
			EulerTour& et2 = euler_tours[v_cc];
			et1.reroot(u);
			et2.reroot(v);
			et1.tour.push_back(std::make_pair(u, v));
			et1.merge(et2.tour.begin(), et2.tour.end());
			et1.tour.push_back(std::make_pair(v, u));
			std::replace(components.begin(), components.end(), v_label(v_cc), v_label(u_cc));
			et2.tour.clear();
			return(true);
		}
	}
	
	// Removes an undirected edge (u,v) using the euler tour representations. 
	// The connected component labels are updated. Returns whether an edge was removed or not.                                                       
	bool cut(size_t u, size_t v){
		if (u > v){ std::swap(u, v); }
		if (!edge_exists(u, v)){ return false; }
		assert(component(u) == component(v));
		
		// Two cases: either (u,v) is a tree edge, which requires modifying the euler tours, or its an auxillary edge
		EulerTour& et1 = euler_tours.at(component(u));
		if (!et1.is_tree_edge(u, v)){
			Rcout << "not a tree edge" << std::endl;
			aux[rank_lex_2(u, v, n_vertices())] = false; 
			return(true);
		} else {
			// Find/ensure the edges (u,v) and (v,u) exist in the current component
			auto e1 = std::find(et1.tour.begin(), et1.tour.end(), edge_t(u, v));
			auto e2 = std::find(et1.tour.begin(), et1.tour.end(), edge_t(v, u));
			assert(e1 != et1.tour.end() && e2 != et1.tour.end());
			
			// Get the indices of the edges (u,v) and (v,u)
			// TODO: change replacement/source component locations depending on whether (u,v) or (v,u) is first
			size_t e1_idx = std::distance(et1.tour.begin(), e1);
			size_t e2_idx = std::distance(et1.tour.begin(), e2);
			
			auto J = vector< edge_t >(et1.tour.begin(), et1.tour.begin() + e1_idx);
			auto K = vector< edge_t >(et1.tour.begin() + e1_idx, et1.tour.begin() + e2_idx + 1);
			auto L = vector< edge_t >(et1.tour.begin() + e2_idx + 1, et1.tour.end());
			
			J.insert(J.end(), L.begin(), L.end());
			euler_tours.at(component(u)) = EulerTour(J.begin(), J.end());
			euler_tours.at(v) = EulerTour(K.begin() + 1, K.end() - 1);
			
			for (auto e: euler_tours.at(v).tour){
				if (e.first == e.second){ components[e.first] = v; }
			}
			
			// // If T = (J, K, L), left-rotate T such that T* = (J, L, K)
			// auto e2_begin = std::rotate(et1.tour.begin() + e1_idx, et1.tour.begin() + e2_idx + 1, et1.tour.end());
			// 
			// // Make new euler tours
			// euler_tours.at(component(u)) = EulerTour(et1.tour.begin(), e2_begin);
			// euler_tours.at(v) = EulerTour(e2_begin, et1.tour.end());
			// 
			// Rprintf("indices1: %d, %d\n", e1_idx, e2_idx);
			// Rcout << "here0" << std::endl;
			// bool swapped = false; 
			// if (e1_idx > e2_idx){ 
			// 	et1.reroot(u);
			// 	e1_idx = std::distance(et1.tour.begin(), std::find(et1.tour.begin(), et1.tour.end(), edge_t(u, v)));
			// 	e2_idx = std::distance(et1.tour.begin(), std::find(et1.tour.begin(), et1.tour.end(), edge_t(v, u)));
			// 	Rprintf("indices2: %d, %d\n", e1_idx, e2_idx);
			// 	// swapped = true; 
			// }
			// Rcout << "here1 " << euler_tours.size() << std::endl;
			// Rprintf("indices3: %d, %d\n", e1_idx, e2_idx);
			// 
			// // Split the current tour into two disjoint components
			// // TODO: use stable_partition instead
			// assert(euler_tours.at(v).tour.size() == 0);
			// euler_tours.at(v) = EulerTour(et1.tour.begin()+e1_idx+1, et1.tour.begin()+e2_idx);
			// Rcout << "here2" << std::endl;
			// 
			// 
			// for (auto e: euler_tours.at(v).tour){
			// 	Rcout << "(" << e.first << ", " << e.second << ") ";
			// }
			// Rcout << std::endl;
			// 
			// Rcout << "here3" << std::endl;
			// auto L = vector< edge_t >(et1.tour.begin()+e2_idx+1, et1.tour.end());
			// et1.tour.resize(e1_idx);
			// // et1.tour.pop_back();
			// et1.merge(L.begin(), L.end());
			// 
			// // Update vertices in v component
			// Rcout << "here4" << std::endl;
			// for (auto e: euler_tours.at(v).tour){
			// 	if (e.first == e.second){ components[e.first] = v; }
			// }
			// 
			// // At this point: et1 is (u) component, et2 is (v) component, and they are disjoint
			// // Need to find a replacement edge joining the two, if it exists
			// 
			// // If the replacement edge was found, reconnect the components
			// Rcout << "here5" << std::endl;
			// std::optional< edge_t > replacement_edge = search_replacement(component(u), v);
			// if (replacement_edge){
			// 	Rcout << "replacement edge found: ";
			// 	const edge_t e = replacement_edge.value();
			// 	Rprintf("%d,%d\n", e.first, e.second);
			// 	bool relinked = link(e.first, e.second);
			// 	return(!relinked);
			// }
		}
	
		return(true);
	}
	
	// Search for a replacement edge between components
	auto search_replacement(const size_t c1, const size_t c2) const -> std::optional< edge_t > {
		if (c1 >= n_vertices() || c2 >= n_vertices()){ return(std::nullopt); }
		if (c1 == c2){ return(std::nullopt); }
			
		// Extract vertices, take their product, stopping at the first edge in aux that exists 
		auto V1 = euler_tours[c1].vertices();
		auto V2 = euler_tours[c2].vertices();
		std::optional< edge_t > replacement_edge = {}; 
		for (auto v1: V1){
			for (auto v2: V2){
				if (aux[rank_lex_2(v1, v2, n_vertices())]){
					replacement_edge = std::make_optional(edge_t(v1, v2));
				}
			}
		}
		return(replacement_edge);
	}
	
	IntegerVector connected_components(){
		IntegerVector cc = IntegerVector(components.size(), 0);
		for (size_t i = 0; i < components.size(); ++i){
			cc[i] = static_cast< int >(components[i].to_ulong());
		}
		return(cc);
	}
	
};

typedef DynamicGraph< 4 > DynamicGraph4;

void print_tours(DynamicGraph4* G){
	size_t i = 0; 
	for (auto et: G->euler_tours){
		if (et.tour.size() > 0){
			Rcout << i << ": ";
			for (auto e: et.tour){ Rcout << "(" << e.first << "," << e.second << ")"; }
			Rcout << std::endl;
		}
		++i;
	}
}



	
RCPP_EXPOSED_CLASS_NODECL(DynamicGraph4)
RCPP_MODULE(dg4_module) {
	using namespace Rcpp;
	class_< DynamicGraph4 >( "DynamicGraph4" )
	.constructor< int >()
	.method( "components", &DynamicGraph4::connected_components )
	.method( "edge_exists", &DynamicGraph4::edge_exists )
	.method( "is_connected", &DynamicGraph4::is_connected )
	.method( "link", &DynamicGraph4::link )
	.method( "cut", &DynamicGraph4::cut )
	.method( "print_tours", &print_tours)
	;
}



// [[Rcpp::export]]
void test_dc(){
	DynamicGraph< 4 > G = DynamicGraph< 4 >(10);
	for (auto label: G.components){ Rcout << static_cast< size_t >(label.to_ulong()) << ", "; }
	Rcout << std::endl;
	
	G.link(0,1);
	G.link(0,7);
	G.link(0,8);
	G.link(8,9);
		
	G.link(2,3);
	G.link(4,5);
	G.link(3,4);
	
	for (auto et: G.euler_tours){
		if (et.tour.size() > 0){
			for (auto e: et.tour){ Rcout << "(" << e.first << "," << e.second << ")"; }
			Rcout << std::endl;
		}
	}
	
	for (auto label: G.components){ Rcout << static_cast< size_t >(label.to_ulong()) << ", "; }
	Rcout << std::endl;
	
	G.cut(3,4);
	G.cut(0,6);
	G.cut(1,7);
	
	Rcout << std::endl;
	for (auto et: G.euler_tours){
		if (et.tour.size() > 0){
			for (auto e: et.tour){ Rcout << "(" << e.first << "," << e.second << ")"; }
			Rcout << std::endl;
		}
	}
	for (auto label: G.components){ Rcout << static_cast< size_t >(label.to_ulong()) << ", "; }
	Rcout << std::endl;
}

/*** R
# el <- combn(6,2)[,sample(seq(choose(6,2)), size = 10)]
el <- matrix(c(0,1,0,2,0,3,1,5,2,4), nrow = 2)
el <- t(kdtools::lex_sort(t(el)))
al <- phtools:::el_to_al_0(el = el, ids = seq(0, 5), symmetric = TRUE)
et <- phtools:::euler_tour(al)

G <- igraph::graph_from_adj_list(lapply(al, function(x) x + 1), mode = "all")
igraph::vertex_attr(G, "label") <- seq(0, 5)
#al1 <- phtools:::el_to_al_0(el = el, ids = seq(6), symmetric = TRUE)
#G <- igraph::graph_from_adj_list(al1, mode = "all")
phtools:::bst(et)



et <- c("a", "b", "c", "b", "e", "b", "a", "d", "f", "g", "f", "d", "h", "i", "h", "d", "a")

is_euler_tour <- function(p){ return(head(p,1L) == tail(p,1L) && length(et) == (2*length(unique(et)) - 1L)) }

## Re-roots a (degree-type) euler tour path 
reroot_et <- function(p, new_root){
	stopifnot(is_euler_tour(p), new_root %in% p)
	if (p[1L] == new_root){ return(p) }
	nr <- which(p == new_root)
	a_idx <- setdiff(1L:head(nr, 1L), head(nr, 1L))
	A <- p[a_idx]
	B <- p[setdiff(seq_along(p), a_idx)]
	return(c(B, A[-1], new_root))
}

## Links two euler tour paths together at (u,v)
link_et <- function(p1, p2, u, v){
	stopifnot(is_euler_tour(p1), is_euler_tour(p2))
	stopifnot(u %in% p1, v %in% p2)
	p1 <- reroot_et(p1, u)
	p2 <- reroot_et(p2, v)
	return(c(p1, p2, u))
}

## Cuts two euler tour paths 
cut_et <- function(p, u, v){
	stopifnot(is_euler_tour(p), u %in% p, v %in% p)
	uv_idx <- which(et %in% c(u, v))
	split_idx <- uv_idx[which(diff(uv_idx) == 1)]
	J <- p[1L:split_idx[1L]]
	K <- p[(split_idx[1L]+1L):split_idx[2L]]
	L <- p[(split_idx[2L]+1L):length(p)]
	return(list(K, c(J[-length(J)], L)))
}


et <- c("a", "b", "c", "b", "a", "d", "a", "e", "f", "e", "a")

reroot_et_edge <- function(p, new_root){
	r_idx <- head(which(et_edges[1,] == new_root), 1L)
	a_idx <- setdiff(seq(1L, r_idx), r_idx)
	b_idx <- setdiff(seq(r_idx, ncol(p)), r_idx)
	# cbind(p[2:1,b_idx], p[,a_idx])
	# return(cbind(p[2:1,rev(a_idx)], p[2:1,b_idx], matrix(p[2:1,r_idx])))
}
et_edges <- rbind(et[1:(length(et)-1L)], et[2:length(et)])




et1 <- head(et[et %in% c("a", "b", "c", "e")], 7)
et2 <- et[et %in% c("d", "f", "g", "h", "i")]


*/
