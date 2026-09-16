#ifndef VGTL_GROW_HPP
#define VGTL_GROW_HPP

/*! \file
 * \brief Reports the star of a simplex
 */

#include <list>
#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/top/star.hpp>

	
namespace vgtl {
	
	using std::list;
	using std::set;

	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T>
	bool
	grow_vc(complex_buffer<T>& cb) {
		set<Vertex(T)> old_vertices;
		set<Cell(T)> new_cells;
		for(typename set<Vertex(T)>::iterator i=cb.template begin<0>();
										i!=cb.template end<0>();++i) {
			Vertex(T) v=*i;
			old_vertices.insert(v);
			list<Cell(T)> st;
			star(cb.t,v,back_inserter(st));
			for(typename list<Cell(T)>::iterator j=st.begin();j!=st.end();++j) {
				Cell(T) ce=*j;
				if(!exists(cb,ce)) {
					put(cb,*j);
					new_cells.insert(ce);
				}
			}
		}
		buffer_set<0,T>(cb).clear();
		for(typename set<Cell(T)>::iterator i=new_cells.begin();
										i!=new_cells.end(); ++i) {
			array<Vertex(T),Dim(T)+1> vs;
			vertices(cb.t,*i,vs);
			for(int j=0;j<=Dim(T);++j)
				if(old_vertices.count(vs[j])==0) 
					put(cb,vs[j]);
		}
		return new_cells.size()!=0;
	}
	

	/*! @} */

}

#endif // VGTL_GROW_HPP
