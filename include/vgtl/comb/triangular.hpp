#ifndef VGTL_TRIANGULAR_HPP
#define VGTL_TRIANGULAR_HPP

/*! \file
 * \brief Triangular traverser
 */

#include <vgtl/utl/array.hpp>

namespace vgtl {

  template <class NumericIt>
	bool
  next_triangular(dim_t m, NumericIt v, NumericIt end) {
    const dim_t N=end-v-1;
  	if(v[N]==m) return false; 
    dim_t r;
    for(r=N-1;r>=0;--r) if(v[r]!=0) break;
    v[r]=v[r]-1;
    for(dim_t j=r+1;j<=N;++j) v[j]=0;
    v[r+1]=m;
    for(dim_t j=0;j<=r;++j) v[r+1]-=v[j];
  	return true; 
	}
	
	/*! \addtogroup comb 
	 * @{
	 */
	
	//! Triangular traverser
	/*! Traverses the "triangular" combinations.
	 *
	 * \b Example:
	 * \include triangular.cpp 
	 * \b Output:
	 * \include triangular.out
	 */
  template <dim_t n, class T=dim_t>
	class triangular_traverser {
    dim_t m;
		array<T,n> v;
		public:
    triangular_traverser(dim_t _m) : m(_m)
		{
			v.assign(0);
			v[0]=m;
		}
		bool operator++() {
  		return next_triangular(m,v.begin(),v.end()); 
		}
		array<T,n>& operator*() {
			return v;
		}
		const array<T,n>& operator*() const {
			return v;
		}
	};
	
	/*! @} */

}

#endif // VGTL_TRIANGULAR_HPP
