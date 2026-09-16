#ifndef VGTL_VERTICES_INV_HPP
#define VGTL_VERTICES_INV_HPP

/*! \file
 * \brief Reports the star of a simplex
 */

#include <list>
#include <vgtl/top/star.hpp>
#include <vgtl/top/incidence.hpp>

	
namespace vgtl {
	
	using std::list;

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Reports the simplex from its vertices 
	template <class T, dim_t k>
	typename disable_if<cell_c<T,k-1>,Simplex(T,k-1)>::type
	vertices_inv(const T& t, const array<Vertex(T),k>& vs) {
		list<Cell(T)> st;
		star(t,vs[0],back_inserter(st));
		Simplex(T,k-1) r;
		for(typename list<Cell(T)>::iterator i=st.begin();i!=st.end();++i) {
			Cell(T) cv=*i;
			face_traverser<T,Dim(T),k-1> ft(t,cv);
			do {
				r=*ft;
				array<Vertex(T),k> tvs;
				vertices(t,r,tvs);
				if(tvs==vs) return r;
			} while(++ft);
		}
		empty_set(r);
		return r;
	}

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k-1>,Simplex(T,k-1)>::type
	vertices_inv(const T& t, const array<Vertex(T),k>& vs) {
		list<Cell(T)> st;
		star(t,vs[0],back_inserter(st));
		Simplex(T,k-1) r;
		for(typename list<Cell(T)>::iterator i=st.begin();i!=st.end();++i) {
			Cell(T) cv=*i;
			array<Vertex(T),k> tvs;
			vertices(t,cv,tvs);
			if(tvs==vs) return cv;
		}
		empty_set(r);
		return r;
	}

	/*! @} */

}

#endif // VGTL_VERTICES_INV_HPP
