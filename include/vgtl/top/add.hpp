#ifndef VGTL_ADD_HPP
#define VGTL_ADD_HPP

/*! \file
 * \brief Generic functions to add simplices and set its faces
 */
	
#include <utility>
#include <vgtl/top/sc.hpp>
#include <vgtl/utl/array.hpp>

namespace vgtl {

  using std::pair;
	
	/*! \addtogroup top 
	 * @{
	 */

	template <class T, dim_t kk>
	Simplex(T,kk-1)
	add(T& t, 
			const array<Simplex(T,kk-2),kk>& fs,
			asc_tag) {
		Simplex(T,kk-1) sd;
		add(t,sd);
		for(dim_t i=0; i<kk; ++i) {
			face_op_set(t,sd,i,fs[i]);
		}
		return sd;
	}
	
	template <class T, dim_t kk>
  typename disable_if<cell_c<T,kk-1>,Simplex(T,kk-1)>::type
	add(T& t, 
			const array<Simplex(T,kk-2),kk>& fs,
			apm_tag) {
		Simplex(T,kk-1) sd;
		add(t,sd);
		for(dim_t i=0; i<kk; ++i) {
			face_op_set(t,sd,i,fs[i]);
			if(empty(up_simplex(t,fs[i]))) up_simplex_set(t,fs[i],sd);
		}
		return sd;
	}
	
	template <class T, dim_t kk>
  typename enable_if<cell_c<T,kk-1>,Simplex(T,kk-1)>::type
	add(T& t, 
			const array<Simplex(T,kk-2),kk>& fs,
			apm_tag ) {
		Simplex(T,kk-1) sd;
		add(t,sd);
		for(dim_t i=0; i<kk; ++i) {
			face_op_set(t,sd,i,fs[i]);
			pair<Cell(T),Cell(T)> p=cells(t,fs[i]);
			if(empty(p.first)&&empty(p.second)) cells_set(t,fs[i],sd,sd);
			else if(p.second==p.first) cells_set(t,fs[i],p.first,sd);
			//else throw exception;
		}
		return sd;
	}
						
	//! Adds a simplex and sets its faces
	/*! When \c t models a Mutable ASC, it just 
	 *  adds a new simplex and sets its faces to the 
	 *  simplices in \c fs. When \c t models a Mutable APM, there is
	 *  two cases: if the new simplex is not a cell 
	 *  it sets the up simplex of each <code>fs[i]</code>
	 *  to the new simplex. Otherwise, it sets  
	 *  at least one of the cells of <code>fs[i]</code>
	 *  to the new simplex.
	 * \post <code>face_op(t,add(t,fs),i)==fs[i]</code>
	 */
	template <class T, dim_t kk>
	Simplex(T,kk-1)
	add(T& t, 
			const array<Simplex(T,kk-2),kk>& fs) {
		return add(t,fs,typename T::sc_category());
	}
			
	//! Adds a vertex to \c t 
	template <class T>
	Vertex(T)
	add(T& t) {
		Vertex(T) s;
		return add(t,s);
	}

	/*! @}*/
			
}

#endif // VGTL_ADD_HPP
