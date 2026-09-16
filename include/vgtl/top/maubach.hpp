#ifndef VGTL_MAUBACH_HPP
#define VGTL_MAUBACH_HPP

/*! \file
 * \brief Maubach Scheme
 */

#include <list>
#include <queue>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/ord_split.hpp>
#include <vgtl/top/star.hpp>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/del.hpp>
#include <vgtl/top/weld.hpp>

namespace vgtl {

	using std::list;
	using std::queue;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T, class ApplyNew, class ApplyOld, class StarIt>
	void
	maubach_edge_split(T& t, Edge(T) e, dim_t k,
										 StarIt begin, StarIt end, ApplyNew& app,
										 ApplyOld& appold, apm_tag) {
		complex_buffer<T> cb_in(t), cb(t);
		Vertex(T) v=app.add_edge_vertex(t,e);
		ord_split(t,e,begin,end,v,k,cb_in,cb);
		apply(cb,app);
		apply(cb_in,appold);
		del(cb_in);
	}
	
	template <class T, class ApplyNew, class ApplyOld, class StarIt>
	void
	maubach_edge_split(T& t, Edge(T) e, dim_t k,
										 StarIt begin, StarIt end, ApplyNew& app,
										 ApplyOld& appold, mt_tag) {
		Vertex(T) v=split_vertex(t,e);
		if(empty(v)) {
			complex_buffer<T> cb_in(t), cb(t);
			v=app.add_edge_vertex(t,e);
			ord_split(t,e,begin,end,v,k,cb_in,cb);
			split_vertex_set(t,e,v);
			split_simplex_set(t,v,e);
			apply(cb,app);
			apply(cb_in,appold);
		} else {
			complex_buffer<T> cb_in(t), cb_notin(t), cb(t);
			decompose(t,e,begin,end,cb_in,cb_notin);
			unlink(cb_in,cb_notin);
			for(StarIt i=begin;i!=end; ++i) {
				Cell(T) cur=*i;
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
	
	template <class T, class ApplyNew, class ApplyOld, class StarIt>
	void
	maubach_edge_split(T& t, Edge(T) e, dim_t k,
										 StarIt begin, StarIt end, ApplyNew& app,
										 ApplyOld& appold) {
		maubach_edge_split(t,e,k,begin,end,app,appold,
										   typename sc_traits<T>::sc_category());
	}

	//! 
	template <class T, class ApplyNew, class ApplyOld>
	void
	maubach_subdivide(T& t, Cell(T) ces, ApplyNew& app, ApplyOld& appold)  {
		dim_t k=Dim(T)-level(t,ces)%Dim(T);
		list<Cell(T)> st;
		queue<Cell(T)> to_visit;
		to_visit.push(ces);
		while(!to_visit.empty()) {
			Cell(T) c=to_visit.front();
			to_visit.pop();
			if(!visited(t,c)) {
				st.push_back(c);
				visit(t,c);
        for(dim_t i=0; i<=Dim(T);++i) {
					if((i==0)||(i==k)) continue;
					Cell(T) ca=adjacent(t,c,i);
					if(!visited(t,ca)) {
						if(level(t,ca)<level(t,c)) {
							maubach_subdivide(t,ca,app,appold);
							ca=adjacent(t,c,i);
						}
						to_visit.push(ca);
					}
				}
			}
		}
		Edge(T) e=face_ind(t,ces,_<dim_t>(0,k));
		maubach_edge_split(t,e,k,st.begin(),st.end(),app,appold);
	}

	//! 
	template <class T, class ApplyNew, class ApplyOld>
	void
	maubach_pre_subdivide(T& t, Cell(T) ces, ApplyNew& app, ApplyOld& appold)  {
		dim_t k=Dim(T)-level(t,ces)%Dim(T);
		list<Cell(T)> st;
		queue<Cell(T)> to_visit;
		to_visit.push(ces);
		while(!to_visit.empty()) {
			Cell(T) c=to_visit.front();
			to_visit.pop();
			if(!visited(t,c)) {
				st.push_back(c);
				visit(t,c);
				for(int i=0; i<=Dim(T);++i) {
					if((i==0)||(i==k)) continue;
					Cell(T) ca=adjacent(t,c,i);
					if(!visited(t,ca)) {
						if(level(t,ca)<level(t,c)) {
							maubach_subdivide(t,ca,app,appold);
							ca=adjacent(t,c,i);
						}
						to_visit.push(ca);
					}
				}
			}
		}
		for(typename list<Cell(T)>::iterator i=st.begin();i!=st.end();++i)
			unvisit(t,*i);
	}

	template <class T, class ApplyNew>
	void
	maubach_subdivide(T& t, Cell(T) ces, ApplyNew& app)  {
		do_nothing dn;
		maubach_subdivide(t,ces,app,dn);
	}

	template <class T, class ApplyNew>
	void
	maubach_pre_subdivide(T& t, Cell(T) ces, ApplyNew& app)  {
		do_nothing dn;
		maubach_pre_subdivide(t,ces,app,dn);
	}

	template <class T>
	Vertex(T)
	maubach_weld_vertex(const T& t, Cell(T) ce) {
		dim_t weld_ind=Dim(T)-((level(t,ce)-1)%Dim(T));
		return face_ind(t,ce,_(weld_ind));
	}
	
	template <class T>
	Edge(T)
	maubach_subdivision_edge(const T& t, Cell(T) ce) {
		dim_t k=Dim(T)-level(t,ce)%Dim(T);
		return face_ind(t,ce,_<dim_t>(0,k));
	}
	
	template <class T, class ApplyOld>
	void
	maubach_weld(T& t, Cell(T) ce, ApplyOld& app)  {
		Vertex(T) v=maubach_weld_vertex(t,ce);
		recursive_weld(t,v,maubach_weld_vertex<T>,app);
	}

	template <class T>
	void
	maubach_weld(T& t, Cell(T) ce)  {
		do_nothing dn;
		Vertex(T) v=maubach_weld_vertex(t,ce);
		recursive_weld(t,v,maubach_weld_vertex<T>,dn);
	}

	template <class T>
	int
	maubach_type(const T& t, Cell(T) ce)  {
		return Dim(T)-level(t,ce)%Dim(T);
	}

	/*! @} */

}

#endif // VGTL_MAUBACH_HPP
