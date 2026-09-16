#ifndef VGTL_MITCHELL_HPP
#define VGTL_MITCHELL_HPP

/*! \file
 * \brief Mitchell Scheme
 */

#include <vgtl/utl/pair_tie.hpp>
#include <vgtl/top/ord_split.hpp>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/del.hpp>
#include <vgtl/top/weld.hpp>

namespace vgtl {

	using std::list;
	using std::set;
	using std::map;
	using std::pair;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T, class ApplyNew>
	void
	mitchell_facet_split(T& t, Facet(T) f, ApplyNew& app, apm_tag) {
		Cell(T) ces[2];
		pair_tie(ces[0],ces[1])=cells(t,f);
		int n=(ces[0]==ces[1]) ? 1 : 2; // n==1 means boundary facet
		complex_buffer<T> cb_in(t), cb(t);
		Vertex(T) v=app.add_facet_vertex(t,f);
		map<Cell(T),Cell(T)> parent;
		ord_split(t,f,ces,ces+n,v,0,cb_in,cb,parent);
		for(typename map<Cell(T),Cell(T)>::iterator i=parent.begin();
										i!=parent.end();++i) {
			app.parent_set(t,i->first,i->second);
		}
		apply(cb,app);
		del(cb_in);
	}
	
	template <class T, class ApplyNew>
	void
	mitchell_facet_split(T& t, Facet(T) f, ApplyNew& app, mt_tag) {
		Cell(T) ces[2];
		pair_tie(ces[0],ces[1])=cells(t,f);
		int n=(ces[0]==ces[1]) ? 1 : 2; // n==1 means boundary facet
		Vertex(T) v=split_vertex(t,f);
		if(empty(v)) {
			complex_buffer<T> cb_in(t), cb(t);
			v=app.add_facet_vertex(t,f);
			ord_split(t,f,ces,ces+n,v,0,cb_in,cb);
			split_vertex_set(t,f,v);
			split_simplex_set(t,v,f);
			apply(cb,app);
		} else {
			complex_buffer<T> cb_in(t), cb_notin(t), cb(t);
			decompose(t,f,ces,ces+n,cb_in,cb_notin);
			unlink(cb_in,cb_notin);
			for(int i=0;i<n; ++i) {
				Cell(T) cur=ces[i];
				Cell(T) p[2];
				p[0]=child(t,cur,0);
				p[1]=child(t,cur,1);
				is_current_set(t,cur,false);
				unvisit(t,cur);
				fill(cb,p[0],v);
				is_current_set(t,p[0],true);
				fill(cb,p[1],v);
				is_current_set(t,p[1],true);
			}
			link(cb);
		}
	}

	template <class T, class ApplyNew>
	void
	mitchell_facet_split(T& t, Facet(T) f, ApplyNew& app) {
		mitchell_facet_split(t,f,app, typename sc_traits<T>::sc_category());
	}
	
	//! 
	template <class T, class ApplyNew>
	void
	mitchell_subdivide(T& t, Cell(T) c, ApplyNew& app)  {
//		Facet(T) f=face_op(t,c,0);
//		Cell(T) ca=adjacent(t,c,f);
//		if(face_op(t,ca,0)!=f) mitchell_subdivide(t,ca,app);
//		mitchell_facet_split(t,f,app);
		Cell(T) ca=adjacent(t,c,0);
		if(level(t,ca)<level(t,c)) mitchell_subdivide(t,ca,app);
		mitchell_facet_split(t,face_op(t,c,0),app);
	}

	template <class T>
	Vertex(T)
	mitchell_weld_vertex(const T& t, Cell(T) ce) {
		return face_ind(t,ce,_<dim_t>(0));
	}
	
	template <class T, class ApplyOld>
	void
	mitchell_weld(T& t, Cell(T) ce, ApplyOld& app)  {
		Vertex(T) v=mitchell_weld_vertex(t,ce);
		recursive_weld(t,v,mitchell_weld_vertex<T>,app);
	}


	/*! @} */

}

#endif // VGTL_MITCHELL_HPP
