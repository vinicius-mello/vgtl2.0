#ifndef VGTL_RADIAL_HPP
#define VGTL_RADIAL_HPP

/*! \file
 * \brief Reports the incident edges of  a vertex in a oriented 2-pseudo-manifold
 */

#include <list>
#include <vgtl/top/incidence.hpp>
	
namespace vgtl {
	
	using std::list;
	using std::pair;

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Reports the incident edges of  a vertex in a oriented 2-pseudo-manifold
	/*! \c T must implement the <code>orientation(t,cv)</code> function.
	 */
	template <class T>
	typename enable_if<equal_c<Dim(T),2>,void>::type
	radial_edges(const T& t,
							 Vertex(T) v,
							 list<Edge(T)>& le) {
		Cell(T) first=up_simplex(t,up_simplex(t,v));
		Edge(T) initial;
		int initial_sign=orientation(t,first);
		for(int i=0;i<3;++i,initial_sign*=-1) {
			initial=face_op(t,first,i);
			if(!in(t,initial,v)) break;
		}
		Cell(T) cur_face;
		Edge(T) cur_edge;
		for(int j=1;j>=0;--j) {
			int k=(initial_sign>0)?j:(1-j);
			cur_face=first;
			cur_edge=opposite(t,cur_face,face_op(t,initial,k));
			while(true) {
				if(j==1) le.push_back(cur_edge);
				else le.push_front(cur_edge);
				pair<Cell(T),Cell(T)> p=cells(t,cur_edge);
				if(p.first==p.second) break;
				if(p.first==cur_face) cur_face=p.second;
				else cur_face=p.first;
				if(cur_face==first) return;
				cur_edge=opposite(t,cur_face,opposite(t,cur_edge,v));
			} 
		}
	}

	/*! @} */

}
	
#endif // VGTL_RADIAL_HPP
