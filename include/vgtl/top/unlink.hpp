#ifndef VGTL_UNLINK_HPP
#define VGTL_UNLINK_HPP

/*! \file
 * \brief Generic functions to correct the links of a simplex
 */

#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/top/link.hpp>

namespace vgtl {
	
	using std::pair;

	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T>
	void
	unlink(complex_buffer<T>& cb_in, complex_buffer<T>& cb_notin) {
		const dim_t dim=Dim(T);
		link(cb_notin);
	  for(typename set<Facet(T)>::iterator i=cb_notin.template begin<dim-1>();
				i!=cb_notin.template end<dim-1>(); ++i) {
			Facet(T) ce=*i;
			pair<Cell(T),Cell(T)> p=cells(cb_notin.t,ce);
			bool a=exists(cb_in,p.first);
			bool b=exists(cb_in,p.second);
			if(a&&b) {
				empty_set(p.first);
				empty_set(p.second);
			} else if(a) {
				p.first=p.second;
			} else if(b) {
				p.second=p.first;
			} // else throw exception;
			cells_set(cb_notin.t,ce,p.first,p.second);
		}
	}
	
	/*! @} */

}

#endif // VGTL_UNLINK_HPP
