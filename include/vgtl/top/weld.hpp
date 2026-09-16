#ifndef VGTL_WELD_HPP
#define VGTL_WELD_HPP

#include <list>
#include <set>
#include <vgtl/top/ord_split.hpp>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/del.hpp>

namespace vgtl {

	using std::list;
	using std::set;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	//! Welds a vertex
	/*! Performs a local vertex weld and applies a function to
	 * the simplices turned current again. <code>[begin,end)</code>
	 * must traverse the vertex star.
	 * \pre <code>is_initial(t,s)==false</code>
	 */
	template <class T, class StarIterator, class Apply>
	void
	weld(T& t, Vertex(T) v,
			 StarIterator begin, StarIterator end,
			 Apply& app) {
		SplitSimplex(T) s=split_simplex(t,v);
		//if(empty(s)) throw exception;
		complex_buffer<T> cb_in(t), cb_notin(t), cb(t);
		decompose(t,v,begin,end,cb_in,cb_notin);
		unlink(cb_in,cb_notin);
	  for(StarIterator i=begin;i!=end; ++i) {
			Cell(T) ce=*i;
			Cell(T) par=parent(t,ce);
			fill(cb,par,s);
			is_current_set(t,ce,false);
			is_current_set(t,par,true);
		}
		link(cb);
		apply(cb,app);
	}

	template <class T, class Apply, class Weld>
	void
	recursive_weld(T& t,
								 Vertex(T) v,
								 Weld weld_vertex,
								 Apply& app) {
		set<pair<int,Vertex(T)> > subv;
		list<Cell(T)> st;
		do {
			subv.clear();
			st.clear();
			star(t,v,back_inserter(st));
			for(typename list<Cell(T)>::iterator i=st.begin(); i!=st.end();++i) {
				Cell(T) cur=*i;
				Vertex(T) wv=weld_vertex(t,cur);
				if(wv!=v) subv.insert(make_pair(-level(t,cur),wv));
			}
			for(typename set<pair<int,Vertex(T)> >::iterator j=subv.begin();
					j!=subv.end(); ++j)
				recursive_weld(t,j->second,weld_vertex,app);
		} while(!subv.empty());
		weld(t,v,st.begin(),st.end(),app);
	}
	

	/*! @} */

}

#endif // VGTL_WELD_HPP
