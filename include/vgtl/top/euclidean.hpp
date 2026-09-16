#ifndef VGTL_EUCLIDEAN_HPP
#define VGTL_EUCLIDEAN_HPP

/*! \file
 * \brief Euclidean points associated to vertices
 */

#include <vgtl/top/incidence.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	//! Euclidean points associated to vertices
	/*!	\c T must implement the <code>euclidean_point(t,s)</code> function. 
	*/ 
	template <class T, dim_t k, class Point>
	void
	euclidean_points(const T& t,
									 Simplex(T,k) s,
									 array<Point,k+1>& p) {
		array<Vertex(T),k+1> v;
		vertices(t,s,v);
		for(dim_t i=0;i<=k;++i) p[i]=euclidean_point(t,v[i]);
	}
	
	/*! @} */

}

#endif // VGTL_EUCLIDEAN_HPP
