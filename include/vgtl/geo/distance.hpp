#ifndef VGTL_GEO_DISTANCE_HPP
#define VGTL_GEO_DISTANCE_HPP

#include <vector>
#include <list>
#include <algorithm>
#include <vgtl/alg/point.hpp>
#include <vgtl/top/star.hpp>
#include <vgtl/top/iso/tubular.hpp>

/*! \file
 * \brief 
 */


namespace vgtl {

	using namespace std;
	
	/*! \addtogroup geo 
	 * @{
	 */

	template <class T, typename Scalar>
	void distance_points(T& t, Scalar m) {
		const dim_t d=Dim(T);
		Vertex_it(T) vi, vend;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			Vertex(T) v=*vi;
			if(!is_current(t,v)) continue;
			Scalar dis=abs(scalar_value(t,v));
			list<Cell(T)> ls;
			star(t,v,back_inserter(ls));
			for(typename list<Cell(T)>::iterator i=ls.begin();i!=ls.end();++i) {
				for(typename vector<point<d,Scalar> >::iterator pi=points_begin(t,*i);
					pi!=points_end(t,*i);++pi) {
					point<d,Scalar> p=*pi;
					vec<d,Scalar> ve=euclidean_point(t,v)-p;
					dis=min(dis,sqrt(dot(ve,ve)));
				}
			}
			dis=max(dis,m);
			scalar_value_set(t,v,signal(t,v)*dis);
		}
	}

	
	template <class T>
	void compute_section(T& t) {
		Edge_it(T) ei,eend;
		for(simplices(t,ei,eend);ei!=eend;++ei) {
			Edge(T) e=*ei;
			if(!is_current(t,e)) continue;
			array<Vertex(T),2> vs;
			vertices(t,e,vs);
			if(!iso::cross(t,vs)) {
				section_set(t,e,-1);
			} else {
				array<point<Dim(T),double>,2> ps;
				ps[0]=euclidean_point(t,vs[0]);
				ps[1]=euclidean_point(t,vs[1]);
				vec<Dim(T),double> v=ps[1]-ps[0];
				double len=sqrt(dot(v,v));
				normalize(v);
				double m=1.0e10;
				list<Cell(T)> st;
				star(t,e,back_inserter(st));
				section_set(t,e,0.5);
				for(typename list<Cell(T)>::iterator i=st.begin();i!=st.end();++i) {
					Cell(T) c=*i;
					for(typename vector<point<Dim(T),double> >::iterator j=points_begin(t,c);
													j!=points_end(t,c);++j) {
						point<Dim(T),double> p=*j;
						double pr=dot(p-ps[0],v);
						vec<3,double> l=pr*v-(p-ps[0]);
						double ll=dot(l,l);
						if(ll<m) {
							m=ll;
							double s=1.0-pr/len;
							s=min(s,0.99);
							s=max(s,0.01);
							section_set(t,e,s);
						}
					}
				}
			}
		}
	}

	/*! @} */

}

#endif // VGTL_GEO_DISTANCE_HPP
