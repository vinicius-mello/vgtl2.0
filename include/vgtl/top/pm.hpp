#ifndef VGTL_PM_HPP
#define VGTL_PM_HPP

/*! \file
 * \brief Pseudo-manifold functions
 */

#include <utility>
#include <vgtl/top/sc.hpp>

namespace vgtl {
	
	using std::pair;
	using std::make_pair;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	//! 
	template <class T>
	Cell(T) 
	adjacent(const T& t, Cell(T) ce, Facet(T) f) {
		pair<Cell(T),Cell(T)> ces=cells(t,f);
		if(ces.first==ce) return ces.second;
		if(ces.second==ce) return ces.first;
		Cell(T) e;
		empty_set(e);
		return e;
	}
	
	template <class T>
	Cell(T) 
	adjacent(const T& t, Cell(T) ce, int i) {
		pair<Cell(T),Cell(T)> ces=cells(t,face_op(t,ce,i));
		if(ces.first==ce) return ces.second;
		if(ces.second==ce) return ces.first;
		return ces.first; // never, just to avoid warnings
	}
	
	//! Tests if \c ce belongs to the boundary
	template <class T>
	bool
	boundary(const T& t,
		       Facet(T) f) {
		pair<Cell(T),Cell(T)> p=cells(t,f);
		return p.first==p.second;
	}
	
	template <class T>
	pair<int,int> 
	common_facet_ind(const T& t, Cell(T) ce0, Cell(T) ce1) {
		int i,j;
		for(i=0;i<=Dim(T);++i) {
			for(j=0;j<=Dim(T);++j) {
				if(face_op(t,ce0,i)==face_op(t,ce1,j)) return make_pair(i,j);
			}
		}
		//throw exception
	}

	/*! @} */

}

#endif // VGTL_PM_HPP
