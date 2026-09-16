#ifndef VGTL_ISO_HPP
#define VGTL_ISO_HPP

#include <list>
#include <vgtl/utl/pair_tie.hpp>
#include <vgtl/utl/array_indices.hpp>
#include <vgtl/alg/vec.hpp>
#include <vgtl/top/incidence.hpp>

namespace vgtl {

	namespace iso {

		using std::pair;
		using std::list;
		using std::make_pair;
		
		template <class T, dim_t kk>
		pair<int,int> dim(const T& t, const array<Vertex(T),kk>& vs) {
			const dim_t k=kk-1;
			int p=0;
      for(dim_t i=0;i<=k;++i)
				if(signal(t,vs[i])>0) ++p;
			if((p==0)||(p==(k+1))) return make_pair(-1,-1);
			return make_pair(p-1,k-p);
		}

		template <class T, dim_t k>
		pair<int,int> dim(const T& t, Simplex(T,k) s) {
			int p=0;
			face_traverser<T,k,0> ft(t,s);
			do {
				if(signal(t,*ft)>0) ++p;
			} while(++ft);
			if((p==0)||(p==(k+1))) return make_pair(-1,-1);
			return make_pair(p-1,k-p);
		}
	
		template <class T, dim_t k>
		Simplex(T,k-1)
		face_op_m(const T& t, Simplex(T,k) s, int i) {
			array<Vertex(T),k+1> vs;
			vertices(t,s,vs);
			int j=0;
			while(true) {
				if(signal(t,vs[j])<0) {
					if(i==0) return face_op(t,s,j);
					--i;	
				}
				++j;
			}
		}
		
		template <class T, dim_t k>
		Simplex(T,k-1)
		face_op_p(const T& t, Simplex(T,k) s, int i) {
			array<Vertex(T),k+1> vs;
			vertices(t,s,vs);
			int j=0;
			while(true) {
				if(signal(t,vs[j])>0) {
					if(i==0) return face_op(t,s,j);
					--i;	
				}
				++j;
			}
		}
		
		template <class T>
		void
		triangulate(const T& t, Edge(T) s, 
								list<array<Edge(T),1> >& ls) {
			array<Edge(T),1> a;
			a[0]=s;
			ls.push_back(a);
		}
		
		template <class T, dim_t k>
		void
		triangulate(const T& t, Simplex(T,k) s, 
								list<array<Edge(T),k> >& ls) {
			int m,p;
			pair_tie(p,m)=dim(t,s);
			if(m==0 || p==0) {
				face_traverser<T,k,1> ft(t,s);
				array<Edge(T),k> a;
				int i=0;
				do {
					Edge(T) e=*ft;
					array<Vertex(T),2> v;
					vertices(t,e,v);
					if(signal(t,v[0])*signal(t,v[1])<0) {
						a[i]=e;++i;
					}
				} while(++ft);
				ls.push_back(a);
			} else {
				list<array<Edge(T),k-1> > lts;
				triangulate(t,face_op_p(t,s,0),lts);
				triangulate(t,face_op_m(t,s,0),lts);
				face_traverser<T,k,1> ft(t,s);
				Edge(T) e;
				do {
					e=*ft;
					array<Vertex(T),2> v;
					vertices(t,e,v);
					if(signal(t,v[0])*signal(t,v[1])<0) {
						break;
					}
				} while(++ft);
				for(typename list<array<Edge(T),k-1> >::iterator i=lts.begin();
						i!=lts.end(); ++i) {
					array<Edge(T),k> a;
          for(dim_t j=0;j<(k-1); ++j) a[j+1]=(*i)[j];
					a[0]=e;
					ls.push_back(a);
				}
			}
		}
		
		template <class T>
		void edge_indices(const T& t,
											const array<Vertex(T),Dim(T)+1>& vs,
											const array<Edge(T),Dim(T)>& a,
                      array<array<dim_t,2>,Dim(T)>& ai) {
			const dim_t d=Dim(T);
      for(dim_t i=0; i<d; ++i) {
				array<Vertex(T),2> e;
				vertices(t,a[i],e);
				indices(e,vs,ai[i]);
			}
		}
		
		template <class T>
		int orientation(const T& t,
										const array<Vertex(T),Dim(T)+1>& vs,
                    const	array<array<dim_t,2>,Dim(T)>& ai) {
			const dim_t d=Dim(T);
      array<vec<d+1,dim_t>,d+1> m;
			m[0][0]=1;
      for(dim_t i=0; i<d; ++i) {
				m[i+1][ai[i][0]]=1;
				m[i+1][ai[i][1]]=1;
			}
			return -signal(t,vs[0])*det(m);
		}
		
		template <class T>
		void basis(const T& t,
							 Cell(T) ce,
							 const array<Vertex(T),Dim(T)+1>& vs,
               array<array<dim_t,2>,Dim(T)>& ai) {
      array<dim_t,2> ae;
			ae[0]=0;ae[1]=1;
      dim_t i=0;
			do {
				if(signal(t,vs[ae[0]])*signal(t,vs[ae[1]])<0) {
					ai[i]=ae;
					++i;
				}
			} while((i<Dim(T))&&next_combination(Dim(T)+1,ae.begin(),ae.end()));
			if(orientation(t,ce)*orientation(t,vs,ai)<0) {
				swap(ai[Dim(T)-1],ai[Dim(T)-2]);
			}
		}
		
		template <class T, dim_t k>
		void basis(const T& t,
							 const array<Vertex(T),k+1>& vs,
               array<array<dim_t,2>,k>& ai) {
      array<dim_t,2> ae;
			ae[0]=0;ae[1]=1;
      dim_t i=0;
			do {
				if(signal(t,vs[ae[0]])*signal(t,vs[ae[1]])<0) {
					ai[i]=ae;
					++i;
				}
			} while((i<k)&&next_combination(k+1,ae.begin(),ae.end()));
		}

		template <class T>
		int orientation(const T& t, Cell(T) s,
										const array<Edge(T),Dim(T)>& a) {
			const dim_t d=Dim(T);
			array<Vertex(T),d+1> vs;
			vertices(t,s,vs);
      array<array<dim_t,2>,d> ai;
			edge_indices(t,vs,a,ai);
			return orientation(t,vs,ai);
		}
		
		template <class T, dim_t k>
		void
		triangulate(const T& t, Simplex(T,k) s, 
								list<array<Edge(T),k> >& ls, list<int>& ori) {
			triangulate(t,s,ls);
			array<Vertex(T),Dim(T)+1> vs;
			vertices(t,s,vs);
      array<array<dim_t,2>,Dim(T)> ai;
			for(typename list<array<Edge(T),k> >::iterator li=ls.begin();
					li!=ls.end(); ++li) {
				edge_indices(t,vs,*li,ai);
				ori.push_back(orientation(t,s)*orientation(t,vs,ai));
			}
		}

	}

}

#endif // VGTL_ISO_HPP
