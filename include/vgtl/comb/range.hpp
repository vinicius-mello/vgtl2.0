#ifndef VGTL_RANGE_HPP
#define VGTL_RANGE_HPP

/*! \file
 * \brief Range traverser
 */

#include <vgtl/utl/array.hpp>

namespace vgtl {


	/*! \addtogroup comb 
	 * @{
	 */
	
	//! Range traverser
	/*! Traverses a range of values.
	 *
	 * \b Example:
	 * \include range.cpp 
	 * \b Output:
	 * \include range.out
	 */
  template <dim_t n>
	class range_traverser {
		array<int,n> v;
		array<int,n> b;
		public:
		range_traverser(const array<int,n>& _b) : b(_b) {
			for(int j=0;j<n;++j) v[j]=0;
		}
		bool operator++() {
			for(int j=0;j<n;++j) {
				++v[j];
				if(v[j]<b[j]) return true;
				v[j]=0;
			}
			return false;
		}
		array<int,n>& operator*() {
			return v;
		}
		const array<int,n>& operator*() const {
			return v;
		}
	};
	
	/*! @} */

}

#endif // VGTL_RANGE_HPP
