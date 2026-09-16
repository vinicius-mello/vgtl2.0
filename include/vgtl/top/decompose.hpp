#ifndef VGTL_DECOMPOSE_HPP
#define VGTL_DECOMPOSE_HPP

/*! \file
 * \brief Generic function that decompose the star of a simplex
 */

#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/top/fill.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	//!* Decomposes the star of a simplex
	/*! Decomposes the star of <code>s</code> in two disjoint sets
	 * <code>cb_in</code> and <code>cb_notin</code> such that 
	 * <code>in(t,f,s)==true</code> for each \c f in \c cb_in.
	 */
	template <class T, dim_t k, class StarIterator>
	void
	decompose(T& t,
						Simplex(T,k) s, 
						StarIterator begin, StarIterator end,
						complex_buffer<T>& cb_in, complex_buffer<T>& cb_notin) {
		for(StarIterator i=begin;i!=end; ++i) {
			Cell(T) cur=*i;
			for(dim_t j=0; j<=Dim(T); ++j) {
				Facet(T) f=face_op(t,cur,j);
				if(!in(t,f,s)) {
					fill(cb_notin,f);
				} else {
					fill(cb_in,f,s);
				}
			}
			put(cb_in,cur);
		}
	}
	
	/*! @} */
	
}

#endif // VGTL_DECOMPOSE_HPP
