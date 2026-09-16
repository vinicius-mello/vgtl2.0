#ifndef VGTL_TOP_CELL_HPP
#define VGTL_TOP_CELL_HPP

/*! \file
 * \brief Top cell
 */

#include <vgtl/top/sc.hpp>

namespace vgtl {

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Returns a cell which contains \c s
	/*! \invariant <code>in(t,top_cell(t,s),s)==true</code>
	 */
	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,Cell(T)>::type
	top_cell(const T& t, Simplex(T,k) s) {
		return top_cell(t,up_simplex(t,s));
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,Cell(T)>::type
	top_cell(const T& t, Simplex(T,k) s) {
		return s;
	}
	
	/*! @} */

}

#endif // VGTL_TOP_CELL_HPP
