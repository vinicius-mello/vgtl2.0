#ifndef VGTL_ARRAY_INDICES_HPP
#define VGTL_ARRAY_INDICES_HPP

/*! \file
 * \brief Array indices
 */

#include <vgtl/utl/array.hpp>

namespace vgtl {


	/*! \addtogroup utl 
	 * @{
	 */

	//! Array indices
  template <class T, dim_t k, dim_t l>
	void indices(const array<T,l>& as, const array<T,k>& ag, 
               array<dim_t,l> &ao) {
    for(dim_t i=0;i<l;++i) {
      //ao[i]=-1;
      for(dim_t j=0;j<k;++j) {
				if(as[i]==ag[j]) {
					ao[i]=j;
					break;
				}
			}
		}
	}
	
	/*! @} */

}

#endif //VGTL_ARRAY_INDICES_HPP
