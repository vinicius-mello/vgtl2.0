#ifndef VGTL_ISOSIMPLICIAL_TOPOLOGY_ESTIMATION_HPP
#define VGTL_ISOSIMPLICIAL_TOPOLOGY_ESTIMATION_HPP

#include <vgtl/top/star.hpp>
#include <vgtl/comb/combinations.hpp>

namespace vgtl {

	using std::list;			

	template <class T>
	double isosimplicial_energy(const T& t) {
		const int d=Dim(T);
		Cell_it(T) i, end;
		double ge=0;
		for(simplices(t,i,end); i!=end; ++i) {
			Cell(T) cur=*i;
			if(!is_current(t,cur)) continue;
			array<int,d+1> ls;
			array<Vertex(T),d+1> vs;
			vertices(t,cur,vs);
			local_signals(t,cur,ls);
			double le=local_confidence(t,cur);
			combinations_lex_traverser<2> cb(d+1);
			do {
				array<unsigned int,2> j=*cb;
				ge+=-ls[j[0]]*ls[j[1]]*signal(t,vs[j[0]])*signal(t,vs[j[1]])*le;		
			} while(++cb);
		}
		return ge;
	}

/*	template <class T>
	double isosimplicial_energy2(const T& t) {
		const int d=Dim(T);
		Edge_it(T) ei, eend;
		double ge=0;
		for(simplices(t,ei,eend); ei!=eend; ++ei) {
			Edge(T) e=*ei;
			if(!is_current(t,e)) continue;
			array<Vertex(T),2> ve;
			vertices(t,e,ve);
			double w=0;
			list<Cell(T)> st;
			star(t,e,back_inserter(st));
			for(typename list<Cell(T)>::iterator i=st.begin(); i!=st.end(); ++i) {
				Cell(T) cur=*i;
				array<int,d+1> ls;
				array<Vertex(T),d+1> vs;
				vertices(t,cur,vs);
				local_signals(t,cur,ls);
				double le=local_confidence(t,cur);
				double de=1;
				for(int j=0;j<=d;++j) 
					if((ve[0]==vs[j])||(ve[1]==vs[j])) de*=ls[j];
				w+=-de*le;
			}
			ge+=signal(t,ve[0])*signal(t,ve[1])*w;
		}
		return ge;
	}*/

	template <class T>
	double isosimplicial_change_energy(const T& t, Vertex(T) v) {
		const int d=Dim(T);
		double ce=0;
		list<Cell(T)> st;
		star(t,v,back_inserter(st));
		for(typename list<Cell(T)>::iterator i=st.begin(); i!=st.end(); ++i) {
			Cell(T) cur=*i;
			array<int,d+1> ls;
			array<Vertex(T),d+1> vs;
			vertices(t,cur,vs);
			local_signals(t,cur,ls);
			double le=local_confidence(t,cur);
			combinations_lex_traverser<2> cb(d+1);
			do {
				array<unsigned int,2> j=*cb;
				if(v==vs[j[0]] || v==vs[j[1]]) 
					ce+=ls[j[0]]*ls[j[1]]*signal(t,vs[j[0]])*signal(t,vs[j[1]])*le;		
			} while(++cb);
		}
		return 2*ce; // each edge must be counted twice
	}
	
	template <class T>
	class isosimplicial_cg {
		Vertex(T) rv;
		T& t;
		public:
		isosimplicial_cg(T& _t) : t(_t) {}
		double energy() {
			return isosimplicial_energy(t);
		}
		double change_configuration() {
			rv=random_vertex(t);
			return isosimplicial_change_energy(t,rv);
		}
		void save_configuration() {
			Vertex_it(T) vi,vend;
			for(simplices(t,vi,vend); vi!=vend; ++vi) {
				Vertex(T) v=*vi;
				if(!is_current(t,v)) continue;
				best_signal_set(t,v,signal(t, v));
			}
			best_signal_flip(t,rv);
		}
		void restore_configuration() {
			Vertex_it(T) vi,vend;
			for(simplices(t,vi,vend); vi!=vend; ++vi) {
				Vertex(T) v=*vi;
				if(!is_current(t,v)) continue;
   			signal_set(t,v,best_signal(t,v));
  		} 
		}
		void accept_configuration() {
			signal_flip(t,rv);
		}
	};

}

#endif //  VGTL_ISOSIMPLICIAL_TOPOLOGY_ESTIMATION_HPP
