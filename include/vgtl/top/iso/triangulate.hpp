#ifndef VGTL_TRIANGULATE_HPP
#define VGTL_TRIANGULATE_HPP

#include <map>
#include <vgtl/utl/array.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/ord_split.hpp>
#include <vgtl/top/del.hpp>
#include <vgtl/top/iso/iso.hpp>

namespace vgtl {

	namespace iso {

		template <class Iso, class T, class Apply/*, class ApplyOld*/>
		void triangulate(const Iso& ti, T& t, Apply& app, /*ApplyOld& appold, */
										bool lin=false) {
			map<Edge(Iso),Vertex(T)> ev;
			Edge_it(Iso) ei, eend;
			for(simplices(ti,ei,eend);ei!=eend;++ei) {
				Edge(Iso) e=*ei;
				if(!is_current(ti,e)) continue;
				array<Vertex(Iso),2> vs; 
				vertices(ti,e,vs);
				if(signal(ti,vs[0])*signal(ti,vs[1])>=0) continue;
				Vertex(T) v=add(t);
				ev[e]=v;
				app.edge_set(t,v,e);
			}
			complex_builder<T> cb(t);
			Cell_it(Iso) ci, cend;
			for(simplices(ti,ci,cend);ci!=cend;++ci) {
				Cell(Iso) ce=*ci;
				if(!is_current(ti,ce)) continue;
				int m,p;
				pair_tie(p,m)=dim(ti,ce);
				if(m==-1 || p==-1) continue;
				list<array<Edge(Iso),Dim(Iso)> > ls;
				list<int> lori;
				triangulate(ti,ce,ls,lori);
				list<Cell(T)> lce;
				typename list<int>::iterator lorii=lori.begin();
				for(typename list<array<Edge(Iso),Dim(Iso)> >::iterator lsit=ls.begin();
						lsit!=ls.end();++lsit,++lorii) {
					array<Vertex(T),Dim(T)+1> vs;
					for(int i=0;i<=Dim(T);++i) vs[i]=ev[(*lsit)[i]];
					Cell(T) c=add(cb,vs);
					app.isocell_set(t,c,ce);
					orientation_set(t,c,*lorii);
					lce.push_back(c);
				}
				if(!lin) {
					complex_buffer<T> cb_in(t), cb_new(t);
					Cell(T) c=*lce.begin();
					if(m==0 || p==0) {
						Vertex(T) v=app.add_cell_vertex(t,c,ce);
						ord_split(t,c,lce.begin(),lce.end(),v,0,cb_in,cb_new);
					} else {
						Edge(T) e=face_ind(t,c,_<dim_t>(0,Dim(T)));
						Vertex(T) v=app.add_edge_vertex(t,e,ce);
						ord_split(t,e,lce.begin(),lce.end(),v,0,cb_in,cb_new);
					}
					for(typename set<Cell(T)>::iterator i=cb_new.template begin<Dim(T)>();
											i!=cb_new.template end<Dim(T)>(); ++i) {
						app.isocell_set(t,*i,ce);
					}
//					apply(cb_new,app);
			//		apply(cb_in,appold);
					del(cb_in);
				}
			}
			Vertex_it(T) vi, vend;
			for(simplices(t,vi,vend);vi!=vend;++vi) {
				Vertex(T) v=*vi;
				app.apply(t,v);
			}
			apply(cb,app);
		}
	}

}

#endif // VGTL_TRIANGULATE_HPP
