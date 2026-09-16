#ifndef VGTL_DO_NOTHING_HPP
#define VGTL_DO_NOTHING_HPP

#include <vgtl/top/sc.hpp>
/*! \file
 * \brief \c do_nothing apply object 
 */

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	//! do nothing apply object
	struct do_nothing {
		template <class T, dim_t k>
		void
		apply(T& t, Simplex(T,k) s) {
		}
	};
	
	/*! @} */

}

#endif // VGTL_DO_NOTHING_HPP
