#include <Rcpp.h>
using namespace Rcpp;

#include "ref.hpp"
#include "iterator_facade.h"
#include "short_alloc.h"
#include "combinadic.h"

#include <cstdint>
#include <memory>
#include <vector>
#include <deque>

using I = uint_fast64_t;
using std::array;
using std::vector;
using std::unique_ptr;


template<typename Iterator, typename TransformFunction>
class output_transform_iterator {
public:
    using iterator_category = std::output_iterator_tag;
 
    explicit output_transform_iterator(Iterator iterator, TransformFunction transformFunction) : iterator_(iterator), transformFunction_(transformFunction) {}
    output_transform_iterator& operator++(){ ++iterator_; return *this; }
    output_transform_iterator& operator++(int){ ++*this; return *this; }
    output_transform_iterator& operator*(){ return *this; }
    template<typename T>
    output_transform_iterator& operator=(T const& value)  {
      *iterator_ = transformFunction_(value);
      return *this;
    }
private:
    Iterator iterator_;
    TransformFunction transformFunction_;
};
 
template<typename TransformFunction>
class output_transformer {
public:
    explicit output_transformer(TransformFunction transformFunction) : transformFunction_(transformFunction) {}
    template<typename Iterator>
    output_transform_iterator<Iterator, TransformFunction> operator()(Iterator iterator) const {
    	return output_transform_iterator<Iterator, TransformFunction>(iterator, transformFunction_);
    }
    
private:
    TransformFunction transformFunction_;
};
 
template<typename TransformFunction>
output_transformer<TransformFunction> make_output_transformer(TransformFunction transformFunction){
  return output_transformer<TransformFunction>(transformFunction);
}

// #include "easy_iterator.h"
// using namespace easy_iterator;



// TODO: stream lexicographic binary searches into bin sort trick to decode combinadics

// Given a filtration [ sorted vector< pair< T, I > > ],

  // permutation vector transforming strictly increasing indices (from lexicographically-ordered simplex combinadics)
  // to the order given by the actual rips filtration

struct counter_it : neo::iterator_facade< counter_it >{
  ptrdiff_t value = 0;
  ptrdiff_t dereference() const { return value; }
  ptrdiff_t distance_to(counter_it o) const { return *o - value; }
  void advance(ptrdiff_t pos) { value += pos; }
  counter_it(ptrdiff_t d = 0) : value(d) {}
};
struct counter {
  ptrdiff_t max_value = 0;
  counter() = delete;
  counter(ptrdiff_t d) : max_value(d){ }
  auto begin() const { return counter_it(0); }
  auto end() const { return counter_it(max_value); }
};

// int getLIS(vector< size_t > &v) {
// 	size_t N = v.size();
// 	vector< size_t > smallest(N); //index = maxlen-1
// 	smallest[0] = v[0];
// 	size_t mxlen = 1;
// 	for(size_t i = 1; i < N; i++) {
// 		if (v[i] > smallest[mxlen-1]) {  // Change Here for nondecreasing
// 			smallest[mxlen++]=v[i]; 
// 		} else {
// 			size_t k = std::lower_bound(smallest.begin(),smallest.begin()+mxlen,v[i])-smallest.begin(); // Change Here for nondecreasing
// 			assert(k >= 0 && k < mxlen);
// 			smallest[k]=v[i];
// 		}
// 	}
// 	return mxlen;  
// }

template < class T, std::size_t BufSize = 16 >
using SmallVec = vector<T, short_alloc< T, BufSize > >;



	
// SmallVec< I >::allocator_type::arena_type a;
// SmallVec< I > simplex{ a }; 
// 
// SmallVec< I >::allocator_type::arena_type b;
// SmallVec< I > face{ b }; 

// TODO: make a constexpr choice between a log(n) search to make the p-chain vs combinadic 

// template< typename OutputIt >
// void boundary(const I sigma, const size_t n, const size_t k, OutputIt out, vector< I > reference){
// 
// 	lex_rank(sigma, )
// }

struct ReductionBase {
	
	virtual size_t dimension() const = 0;
	virtual ~ReductionBase(){};
	// template< size_t k, typename InputIt,  typename OutputIt  >
	// virtual void boundary_ranks(InputIt s, const InputIt e, OutputIt out) const = 0;
};

// Reduction over coefficients in a field. 
// This class encapsulates an implicitly stored R = DV decomposition of a given filtration. 
// The filtration simplices are stored using the combinatorial number system, and the reduction 
// is performed by forming columns of D and R on-the-fly. Boundary columns are stored non-contiguously
// w/ pointers, and row permutation information is maintained to make permutations of the decomposition 
// more efficient.
template < size_t d, typename Field = bool >
struct Reduction : ReductionBase {
	using F = Field; 
	using boundary_ptr = unique_ptr< vector< size_t > >;
	template < typename T > using k_array = array< T, d >;
	
	size_t n; 
	k_array< vector< I > > simplex_ranks;					// each dimensions simplices
	k_array< vector< Field > > coeff;							// coefficients associated with each boundary chain / column 
	k_array< boundary_ptr > V;     								// upper-triangular portion of V 
	k_array< vector< size_t > > row_p;  					// permutation vector to store the order of the rows
	
	// ~Reduction(){}
	size_t dimension() const override {
		return d; 
	};
	// Performs the multiplication R_i = D * V_i 
	// template < size_t k, typename OutputIt >
	// void Dv(const size_t i, OutputIt ri){
		
		// auto& vi = *V[k]; // the stored boundary chain 
		// vector< >
		// for (auto entry: vi){ // entry of entry_t  
		// 	auto index = entry.first;
		// 	
		// }
	// }
	
	// template< size_t k, typename InputIt,  typename OutputIt  >
	// void boundary_ranks(InputIt s, const InputIt e, OutputIt out) {
	// 	static_assert(k < d);
	// 	
	// 	// Collect ranks of the simplices to take the boundaries of
	// 	SmallVec< I >::allocator_type::arena_type a;
	// 	SmallVec< I > domain_ranks{ a };
	// 	for (; s != e; ++s){
	// 		auto i = *s; 
	// 		domain_ranks.push_back(simplex_ranks[k].at(i));
	// 	}
	// 			
	// 	// Yields the ranks of the faces 
	// 	boundary_ranks(domain_ranks.begin(), domain_ranks.end(), k, n, out);
	// }
		
	// Applies the boundary operator to the simplices at indices pointed to by the input iterators 
	// stores the resulting entries (as pairs of indices + coefficients) in the output iterator
	// template< size_t k, typename InputIt,  typename OutputIt  >
	// void boundary(InputIt s, const InputIt e, OutputIt out) {
	// 	static_assert(k < d);
	// 	auto& k_simplex_ranks = simplex_ranks[k];
	// 	
	// 	// Collect ranks + coefficients of the simplices to take the boundaries of
	// 	SmallVec< I >::allocator_type::arena_type a;
	// 	SmallVec< F >::allocator_type::arena_type b;
	// 	SmallVec< I > domain_ranks{ a };
	// 	SmallVec< F > domain_coeff{ b };
	// 	for (; s != e; ++s){
	// 		auto i = *s; 
	// 		domain_ranks.push_back(simplex_ranks[k].at(i));
	// 		domain_coeff.push_back(coeff[k].at(i));
	// 	}
	// 	
	// 	// Augment the output iterator s.t. both the ranks and the actual coefficients are output 
	// 	bool positive = true; 
	// 	size_t c = 0, j = 0; 
	// 	auto const elementary_chain_f = output_transformer([&c](I face_rank){
	// 		auto coefficient = positive ? domain_coeff.at(j) : -domain_coeff.at(j);
	// 		++c; positive = !positive; 
	// 		if (c % k == 0){ ++j; c = 0; }
	// 		return std::make_pair(face_rank, coefficient);
	// 	});
	// 			
	// 	// Yields the ranks of the faces 
	// 	auto chain_output = elementary_chain_f(out);
	// 	boundary_ranks(domain_ranks.begin(), domain_ranks.end(), k, n, chain_output)
	// }
	
	// void reduce(vector< size_t > indices){
	// 	bool is_inc = std::is_sorted(indices.begin(), indices.end());
	// 	if (!is_inc){ throw std::invalid_argument("Column indices to reducemust be increasing."); }
	// 	for ()
	// }
	// void reduce(size_t offset = 0){
	// 	// vector< size_t > indices ;
	// }

	
};
// auto test_dgm = Reduction< 1, bool >();

// void boundary_matrix(ReductionBase* decomp, const size_t d){
// 	
// 	//void* dgm; 
// 	switch(d){
// 		case 0: {
// 			auto* dgm = dynamic_cast< Reduction< 0, bool >* >(decomp);
// 			break;
// 		}
// 		case 1: {
// 			auto* dgm = dynamic_cast< Reduction< 1, bool >* >(decomp);
// 			vector< bool > entries; 
// 			vector< I > i, j; 
// 			const size_t m = dgm->simplex_ranks[0].size();
// 			auto cc = counter(m);
// 			auto i_ins = std::back_inserter(i);
// 			dgm->boundary_ranks< 1 >(cc.begin(), cc.end(), i_ins);
// 			for (auto jj: cc){
// 				j.push_back(jj);
// 				j.push_back(jj);
// 			}
// 			break; 
// 		}
// 		case 2: 	
// 			break; 
// 		default: 
// 			throw std::invalid_argument("Only d <= 2 supported.");
// 			break; 
// 	}
// 	// auto* dgm = dynamic_cast< Reduction< d, bool >* >(decomp);
// 	// dgm->boundary_ranks();
// }

RCPP_MODULE(rv_module) {
	class_< Reduction< 1, bool > >( "RV" )
	;
}



// [[Rcpp::export]]
IntegerVector test(IntegerVector k_simplices, const size_t n, const size_t k){
  // auto counter = until_7_iter();
  // auto cc = counter(10);
  // for (auto i: cc){
  //   Rcpp::checkUserInterrupt();
  //   Rcout << i << std::endl;
  // }
//   auto face_indices = vector< I >();
// 	auto face_out = std::back_inserter(face_indices);
//   boundary(k_simplices.begin(), k_simplices.end(), k, n, face_out);
//   return(wrap(face_indices));
	return(IntegerVector::create());
}

// IntegerVector faces(IntegerVector v){
// 	vector< I > f;
// 	auto out = std::back_inserter(f);
// 	faces(v.begin(), v.end(), out);
// 	return(wrap(f));
// }

/*** R
phtools:::test()
*/




// Lexicographically-refined Rips filtration (fixed vertices)
struct RipsFiltration {
  // struct entry { I index; Field value; }; // index == combinadic,
  using simplex_ranks = vector< I >;

  size_t n;                             // number of vertices
  vector< size_t > p;                   // permutation vector
  vector< simplex_ranks > filtration;   // weighted simplices

  // template < typename Lamdba >
  // void p_simplices(const size_t d){
  //   if (d == 0){
  //     for (auto v: counter(n)){ f(v); }
  //   } else {
  //     // auto& p_ranks = filtration.at(d);
  //   }
  // }

  // template < typename Lambda >
  // void simplices(Lambda f){
  //   for (auto v: counter(n)){ f(v); }
  //
  //   // Iterate through the edges
  //   // if (filtration.size() > 0){
  //   //   auto& edges = filtration.at(0);
  //   //   for (auto e: edges){
  //   //
  //   //   }
  //   // }
  // }

  // TODO: define an iterator to iterate over a given dimension p

  // TODO: make a faces iterator -- Can the combinadic avoid the binary_search over the dimension to get the boundary operator

  // Writes the positional indices of the boundary of the p-simplex at index i to out
  // template < typename OutputIterator >
  // void boundary(const size_t i, const size_t p, OutputIterator out){
  //
  // }



  // class iterator {
  //
  // };
};
