#ifndef VGTL_ADD_SIMPLEX_HPP
#define VGTL_ADD_SIMPLEX_HPP

/*! \file
 * \brief Generic function to add a simplex from its vertices
 */

#include <vgtl/top/complex_builder.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	//! Adds a simplex from its vertices
	/*! Adds a simplex corresponding a each subset of <code>vs</code>.
	 *  The <code>complex_builder</code> is updated and 
	 *  checked in order to avoid repetitions.
	 *  The vertices \c vs are allways supposed ordered.
	 * \post <code>vertices(cb.t,add(cb,vs))==vs</code>
	 */
	template <class T, dim_t kk>
	Simplex(T,kk-1) 
	add(complex_builder<T>& cb,
		  const array<Vertex(T),kk>& key) {
		//if(exists(cb,key)) throw error;
		array<Simplex(T, kk-2),kk> fs;
		for(dim_t i=0; i<kk; ++i) {
			array<Vertex(T),kk-1> vfs;
			for(dim_t j=0; j<i; ++j) vfs[j]=key[j];
			for(dim_t j=i+1; j<kk; ++j) vfs[j-1]=key[j];
			if(exists(cb,vfs)) fs[i]=get(cb,vfs);
			else fs[i]=add(cb,vfs);
		}
		Simplex(T,kk-1) s;
		s=add(cb.t,fs);
	  put(cb,key,s);
		return s;
	}
	
	template <class T>
	Edge(T) 
	add(complex_builder<T>& cb,
		  const array<Vertex(T),2>& key) {
		//if(exists(cb,key)) throw error;
		array<Vertex(T),2> fs;
		fs[0]=key[1]; fs[1]=key[0];
		Edge(T) s=add(cb.t,fs);
	  put(cb,key,s);
		return s;
	}

}

#endif // VGTL_ADD_SIMPLEX_HPP
