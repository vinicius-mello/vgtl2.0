#ifndef VGTL_HYPERCUBE_HPP
#define VGTL_HYPERCUBE_HPP

/*! \file
 * \brief Hypercube traverser
 */

#include <vgtl/utl/array.hpp>

namespace vgtl {

			
	/*! \addtogroup comb 
	 * @{
	 */
	
	//! Hypercube traverser
	/*! Traverses the vertices of a hypercube.
	 *
	 * \b Example:
	 * \include hypercube.cpp 
	 * \b Output:
	 * \include hypercube.out
	 */
  template <dim_t n, class T=int>
	class hypercube_traverser {
    dim_t i;
		public:
		typedef array<T,n> result_type;
		hypercube_traverser()
		{
			i=0;
		}
		bool operator++() {
			++i;
			if(i==(1<<n)) return false;
			return true;
		}
		array<T,n> operator*() const {
			array<T,n> v;
      for(dim_t j=0;j<n;++j)
				v[j]=((1<<j)&i) ? 1 : 0;
			return v;
		}
	};
	
	/*! @} */

}

#endif // VGTL_HYPERCUBE_HPP
