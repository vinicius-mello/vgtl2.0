#ifndef VGTL_BARYCENTRIC_HPP
#define VGTL_BARYCENTRIC_HPP

/*! \file
 *  \brief barycentric subdivision
 */

#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/top/star.hpp>
#include <vgtl/top/ord_split.hpp>
#include <vgtl/top/del.hpp>

namespace vgtl {
	

	/*! \addtogroup top 
	 * @{
	 */

	template <class T, class Apply, dim_t k>
	typename enable_if<edge_c<T,k>,void>::type
	barycentric_subdivide(complex_buffer<T>& cb, Apply& app);

	template <class T, class Apply, dim_t k>
	typename disable_if<edge_c<T,k>,void>::type
	barycentric_subdivide(complex_buffer<T>& cb, Apply& app);

	template <class T, class Apply, dim_t k>
	typename enable_if<edge_c<T,k>,void>::type
	barycentric_subdivide(complex_buffer<T>& cba, Apply& app) {
	  for(typename set<Simplex(T,k)>::iterator i=cba.template begin<k>();
				i!=cba.template end<k>(); ++i) {
			Simplex(T,k) s=*i;
			Vertex(T) v=app.barycentric_vertex(cba.t,s);
			list<Cell(T)> st;
			star(cba.t,s,back_inserter(st));
			complex_buffer<T> cb_in(cba.t);
			complex_buffer<T> cb(cba.t);
			map<Cell(T),Cell(T)> parent;
			ord_split(cba.t,s,st.begin(),st.end(),v,k,cb_in,cb,parent);
			apply(cb,app);
			del(cb_in);
		}
	}
	
	template <class T, class Apply, dim_t k>
	typename disable_if<edge_c<T,k>,void>::type
	barycentric_subdivide(complex_buffer<T>& cba, Apply& app) {
	  for(typename set<Simplex(T,k)>::iterator i=cba.template begin<k>();
				i!=cba.template end<k>(); ++i) {
			Simplex(T,k) s=*i;
			Vertex(T) v=app.barycentric_vertex(cba.t,s);
			list<Cell(T)> st;
			star(cba.t,s,back_inserter(st));
			complex_buffer<T> cb_in(cba.t);
			complex_buffer<T> cb(cba.t);
			map<Cell(T),Cell(T)> parent;
			ord_split(cba.t,s,st.begin(),st.end(),v,k,cb_in,cb,parent);
			apply(cb,app);
			del(cb_in);
		}
		barycentric_subdivide<T,Apply,k-1>(cba,app);
	}
	
	//! Performs barycentric subdivision of \c t. Class \c app must implement
	//\c barycentric_vertex(t,s) that returns the new vertex associated to \c s
	template <class T, class Apply>
	void barycentric_subdivision(T& t, Apply& app) {
		complex_buffer<T> cb(t);
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			fill(cb,c);
		}
		barycentric_subdivide<T,Apply,Dim(T)>(cb,app);
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			level_set(t,c,0);
		}
	}

	template <class T, class Apply, dim_t k>
	typename enable_if<edge_c<T,k>,void>::type
	barycentric_subdivide2(T& t, set<Vertex(T)>& sv, Apply& app);

	template <class T, class Apply, dim_t k>
	typename disable_if<edge_c<T,k>,void>::type
	barycentric_subdivide2(T& t, set<Vertex(T)>& sv, Apply& app);

	template <class T, class Apply, dim_t k>
	typename enable_if<edge_c<T,k>,void>::type
	barycentric_subdivide2(T& t, set<Vertex(T)>& sv, Apply& app) {
		Simplex_it(T,k) si,send;
		list<Simplex(T,k)> ls;
	  for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			array<Vertex(T),k+1> vs;
			vertices(t,s,vs);
			bool flag=true;
			for(int i=0;i<=k;++i) {
				if(sv.count(vs[i])==0) {
					flag=false;
					break;
				}
			}
			if(flag) ls.push_back(s);
		}
		for(typename list<Simplex(T,k)>::iterator i=ls.begin();i!=ls.end();++i) {
			Simplex(T,k) s=*i;
			Vertex(T) v=app.barycentric_vertex(t,s);
			list<Cell(T)> st;
			star(t,s,back_inserter(st));
			complex_buffer<T> cb_in(t);
			complex_buffer<T> cb(t);
			map<Cell(T),Cell(T)> parent;
			ord_split(t,s,st.begin(),st.end(),v,k,cb_in,cb,parent);
			apply(cb,app);
			del(cb_in);
		}
	}
	
	template <class T, class Apply, dim_t k>
	typename disable_if<edge_c<T,k>,void>::type
	barycentric_subdivide2(T& t, set<Vertex(T)>& sv, Apply& app) {
		Simplex_it(T,k) si,send;
		list<Simplex(T,k)> ls;
	  for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			array<Vertex(T),k+1> vs;
			vertices(t,s,vs);
			bool flag=true;
			for(int i=0;i<=k;++i) {
				if(sv.count(vs[i])==0) {
					flag=false;
					break;
				}
			}
			if(flag) ls.push_back(s);
		}
		for(typename list<Simplex(T,k)>::iterator i=ls.begin();i!=ls.end();++i) {
			Simplex(T,k) s=*i;
			Vertex(T) v=app.barycentric_vertex(t,s);
			list<Cell(T)> st;
			star(t,s,back_inserter(st));
			complex_buffer<T> cb_in(t);
			complex_buffer<T> cb(t);
			map<Cell(T),Cell(T)> parent;
			ord_split(t,s,st.begin(),st.end(),v,k,cb_in,cb,parent);
			apply(cb,app);
			del(cb_in);
		}
		barycentric_subdivide2<T,Apply,k-1>(t,sv,app);
	}
	
	template <class T, class Apply>
	void barycentric_subdivision2(T& t, Apply& app) {
		set<Vertex(T)> sv;
		Vertex_it(T) vi,vend;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			Vertex(T) v=*vi;
			if(!is_current(t,v)) continue;
			sv.insert(v);
		}
		barycentric_subdivide2<T,Apply,Dim(T)>(t,sv,app);
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			level_set(t,c,0);
		}
	}

	/*! @} */
	
}

#endif // VGTL_BARYCENTRIC_HPP
