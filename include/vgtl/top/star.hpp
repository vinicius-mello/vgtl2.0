#ifndef VGTL_STAR_HPP
#define VGTL_STAR_HPP

/*! \file
 * \brief Reports the star of a simplex
 */

#include <queue>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/mark.hpp>
#include <vgtl/top/top_cell.hpp>

	
namespace vgtl {
	
	using std::queue;

	/*! \addtogroup top 
	 * @{
	 */
	

	//! Reports the star of a simplex
	/*! Performs a depth first visit on the Neighborhood Graph
	 *  associated to \c s. Works only on locally completely connected
	 *  APM's.
	 */
	template <class T, dim_t k, class Output>
	typename disable_if<facet_or_cell_c<T,k>,void>::type
	star(const T& t, Simplex(T,k) s, Output o) {
		queue<Cell(T)> st;
		queue<Cell(T)> to_visit;
		to_visit.push(top_cell(t,s));
		while(!to_visit.empty()) {
			Cell(T) c=to_visit.front();
			to_visit.pop();
			if(!visited(t,c)) {
				st.push(c);
				visit(t,c);
        for(dim_t i=0;i<=Dim(T);++i) {
					Facet(T) f=face_op(t,c,i);
					Cell(T) ca=adjacent(t,c,f);
					if((!visited(t,ca))&&in(t,f,s)) to_visit.push(ca);
				}
			}  
		}
		while(!st.empty()) {
			Cell(T) c=st.front();
			st.pop();
			unvisit(t,c);
			*o=c; ++o;
		}
	}
	
	template <class T, dim_t k, class Output>
	typename enable_if<facet_c<T,k>,void>::type
	star(const T& t, Simplex(T,k) s, Output o) {
		std::pair<Cell(T),Cell(T)> p=cells(t,s);
		*o=p.first; ++o;
		if(p.first!=p.second) {
			*o=p.second; ++o;
		}
	}
	
	template <class T, dim_t k, class Output>
	typename enable_if<cell_c<T,k>,void>::type
	star(const T& t, Simplex(T,k) s, Output o) {
		*o=s; ++o;
	}
	

	/*! @} */

}

#endif // VGTL_STAR_HPP
