#ifndef VGTL_INITIAL_HPP
#define VGTL_INITIAL_HPP

/*! \file
 * \brief Change the resolution of a BMT
 */

#include <vgtl/top/fill.hpp>
#include <vgtl/top/link.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */
	
	//! Changes \c t to the initial triangulation
	template <class T>
	void initial(T& t) {
		Cell_it(T) i, end;
		for(simplices(t,i,end); i!=end; ++i) {
			Cell(T) cur=*i;
			if(has_parent(t,cur)) is_current_set(t,cur,false);
			else {
				complex_buffer<T> cb(t);
				fill(cb,cur);
				link(cb);
				is_current_set(t,cur,true);
			}
		}
	}

	//! Changes \c t to the final triangulation
	template <class T>
	void final(T& t) {
		Cell_it(T) i, end;
		for(simplices(t,i,end); i!=end; ++i) {
			Cell(T) cur=*i;
			if(has_children(t,cur)) is_current_set(t,cur,false);
			else {
				complex_buffer<T> cb(t);
				fill(cb,cur);
				link(cb);
				is_current_set(t,cur,true);
			}
		}
	}

	//! Changes \c t to the level \e l triangulation
	template <class T>
	void set_level(T& t, int l) {
		Cell_it(T) i, end;
		for(simplices(t,i,end); i!=end; ++i) {
			Cell(T) cur=*i;
			int lev=level(t,cur);
			if((lev==l)||((lev<l)&&(!has_children(t,cur)))) {
				complex_buffer<T> cb(t);
				fill(cb,cur);
				link(cb);
				is_current_set(t,cur,true);
			}
			else is_current_set(t,cur,false);
		}
	}

	//! Fills the \c initial_cells(t) container with the current top cells
	template <class T>
	void initial_cells_set(T& t) {
		Cell_it(T) i, end;
		for(simplices(t,i,end); i!=end; ++i) {
			Cell(T) cur=*i;
			initial_cells(t).push_back(cur);
		}
	}
	/*! @} */

}

#endif // VGTL_INITIAL_HPP
