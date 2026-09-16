#ifndef VGTL_MARK_HPP
#define VGTL_MARK_HPP

/*! \file
 * \brief Mark functions
 */

#include <vgtl/top/sc.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */
	
	//! 
	template <class T>
	void 
	visit(const T& t, Cell(T) ce) {
		mark_set(t,ce,black_mark);
	}
	
	template <class T>
	void 
	unvisit(const T& t, Cell(T) ce) {
		mark_set(t,ce,white_mark);
	}
	
	template <class T>
	bool 
	visited(const T& t, Cell(T) ce) {
		return (mark(t,ce)!=white_mark);
	}
	
	/*! @} */

}

#endif // VGTL_MARK_HPP
