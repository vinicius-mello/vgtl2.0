#ifndef VGTL_COMBINATIONS_HPP
#define VGTL_COMBINATIONS_HPP

/*! \file
 * \brief Combinations traverser
 */

#include <vgtl/utl/array.hpp>

namespace vgtl {


	template <class NumericIt>
	bool
  next_combination(dim_t m, NumericIt v, NumericIt end) {
    dim_t n=end-v;
		if(v[0]==(m-n)) return false;
    dim_t i;
		for(i=n-1; i>0; --i) if(v[i]!=(m-n+i)) break;
		v[i]+=1;
		for(i=i+1;i<n;++i) v[i]=v[i-1]+1;
		return true;
	}

	template <class NumericIt>
	bool
  prev_combination(dim_t m, NumericIt v, NumericIt end) {
    dim_t n=end-v;
		if(v[n-1]==(n-1)) return false;
    dim_t i;
		for(i=n-1; i>0; --i) if((v[i-1]+1)!=v[i]) break;
		v[i]-=1;
		for(i=i+1;i<n;++i) v[i]=m-n+i;
		return true;
	}
	
	/*! \addtogroup comb 
	 * @{
	 */
	
	//! Combinations traverser
	/*! Traverses the combinations (m,n) in lexicographic order.
	 *
	 * \b Example:
	 * \include combinations.cpp 
	 * \b Output:
	 * \include combinations.out
	 */
  template <dim_t n>
	class combinations_lex_traverser {
		unsigned int m;
    array<dim_t,n> v;
		public:
    typedef array<dim_t,n> result_type;
		combinations_lex_traverser(unsigned int _m) : m(_m) 
		{
			for(unsigned int i=0; i<n; ++i) v[i]=i;
		}
		bool operator++() {
			return next_combination(m,v.begin(),v.end());
		}
    array<dim_t,n>& operator*() {
			return v;
		}
    const array<dim_t,n>& operator*() const {
			return v;
		}
	};
	
	//! Combinations traverser
	/*! Traverses the combinations (m,n) in reverse lexicographic order.
	 */
  template <dim_t n>
	class combinations_lex_rev_traverser {
		unsigned int m;
    array<dim_t,n> v;
		public:
		combinations_lex_rev_traverser(unsigned int _m) : m(_m) 
		{
			for(unsigned int i=(m-n); i<m; ++i) v[i-m+n]=i;
		}
		bool operator++() {
			return prev_combination(m,v.begin(),v.end());
		}
    array<dim_t,n>& operator*() {
			return v;
		}
    const array<dim_t,n>& operator*() const {
			return v;
		}
	};
	
	/*! @} */

}

#endif // VGTL_COMBINATIONS_HPP
