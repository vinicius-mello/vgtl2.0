#ifndef VGTL_GEO_UNIFORM_TRIANGLE_HPP
#define VGTL_GEO_UNIFORM_TRIANGLE_HPP

#include <map>
#include <vgtl/top/complex_builder.hpp>
#include <vgtl/alg/point.hpp>
#include <vgtl/comb/triangular.hpp>

namespace vgtl {
		
	using namespace std;

	template <class Tri>
	void uniform_triangle_set(Tri& tri, int r) {
		map<point<3,int>,Vertex(Tri)> mp;
		{
			triangular_traverser<3,int> tt(r);
			do {
				point<3,int> p=*tt;
				Vertex(Tri) v2=add(tri);
				coord_set(tri,v2,p);
				mp[p]=v2;
			} while(++tt);
		}
		complex_builder<Tri> cb2(tri);
		{
			array<vec<3,int>,3> e;
			e[0][0]=1;
			e[1][1]=1;
			e[2][2]=1;
			triangular_traverser<3,int> tt(r-1);
			do {
				point<3,int> b=*tt;
				array<Vertex(Tri),3> vs;
				for(int j=0;j<3;++j) {
					point<3,int> be=(b+e[j]);
					vs[j]=mp[be];
				}
				orientation_set(tri,add(cb2,vs),+1);
			} while(++tt);
		}
		{
			array<vec<3,int>,3> e;
			e[0][0]=1; e[0][1]=1;
			e[1][0]=1; e[1][2]=1;
			e[2][1]=1; e[2][2]=1;
			triangular_traverser<3,int> tt(r-2);
			do {
				point<3,int> b=*tt;
				array<Vertex(Tri),3> vs;
				for(int j=0;j<3;++j) {
					point<3,int> be=(b+e[j]);
					vs[j]=mp[be];
				}
				orientation_set(tri,add(cb2,vs),-1);
			} while(++tt);
		}
	}

}

#endif // VGTL_GEO_UNIFORM_TRIANGLE_HPP
