#ifndef VGTL_DEL_HPP
#define VGTL_DEL_HPP

/*! \file
 * \brief deletes the simplices in a complex_buffer 
 */
 
#include <vgtl/top/complex_buffer.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	struct do_delete {
		template <class T, dim_t k>
		void
		apply(T& t, Simplex(T,k) s) {
			del(t,s);
		}
	};
	
	//! Deletes the simplices of \c t in \c cb
	/*!	\c T must implement the <code>del(t,s)</code> function. 
	 */
	template <class T>
	void
	del(complex_buffer<T>& cb) {
		do_delete dd;
		apply(cb,dd);
	}
	
	/*! @} */
}

#endif // VGTL_DEL_HPP
