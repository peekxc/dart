#include <cmath>
#include <bitset>
#include <bit>
#include <optional>

using std::vector; 
using std::bitset;
using std::pair;


// Lexicographically rank 2-subsets
constexpr auto rank_lex_2(size_t i, size_t j, const size_t n) noexcept {
  if (j < i){ std::swap(i,j); }
  return size_t(n*i - i*(i+1)/2 + j - i - 1);
}

inline std::pair< size_t, size_t > unrank_lex_2(const size_t r, const size_t n) noexcept  {
	auto i = static_cast< size_t >( (n - 2 - floor(sqrt(-8*r + 4*n*(n-1)-7)/2.0 - 0.5)) );
	auto j = static_cast< size_t >( r + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 );
	return std::make_pair(i,j);
}


template< size_t w >
struct cutset {
	using word_t = bitset< w >; // w := ceiling(log2(n))
	using edge_t = pair< word_t, word_t >; 
	const size_t n; 
	// const size_t nb; nb(std::ceil(std::log2(nv))) 
	
	vector< edge_t > vertices; // 
	vector< bool > edges; 		 // adjacency matrix (lower-triangular) for verification
		
	cutset(const size_t nv) : n(nv), vertices(nv), edges(nv*(nv-1)/2, false) {
		static_assert(w == std::ceil(std::log2(nv)));
	}
	
	void insert_edge(const size_t u, const size_t v){
		assert(u < n && v < n);
		vertices[u].first ^= word_t(u);
		vertices[v].second ^= word_t(v);
	}
	
	template< typename Iter >
	void insert(Iter b, const Iter e){
		for(; b != e; b += 2){ insert_edge(*b, *(b+1)); }
	}
	
	// Verifies the existence of an edge (u, v)
	constexpr bool verify(const size_t u, const size_t v) const {
		return edges[rank_lex_2(u,v,n)];
	}
	
	// Searches the set of vertices defined by the range T = [b,e) for a valid edge. 
	template < typename Iter >
	auto search(Iter b, const Iter e) -> std::optional< pair< size_t, size_t > >{
		edge_t result = std::make_pair(word_t(0), word_t(0));
		std::for_each(b, e, [this, &result](auto v){
			result.first ^= vertices[v].first;
			result.second ^= vertices[v].second;
		});
		const bool is_empty = (result.first == word_t(0)) && (result.second == word_t(0));
		return(is_empty ? std::nullopt : std::make_optional(result));
	}
	
};
