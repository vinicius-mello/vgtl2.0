#include <iostream>
#include <vgtl/top/model/lc.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/radial.hpp>

using namespace std;
using vgtl::dim_t;
using vgtl::array;

typedef vgtl::lc<2> T;

T t;

namespace vgtl {
	template <>
	struct extra_data<0> {
		int id;
	};

	template <>
	struct extra_data<2> {
		int o;
	};

	int id(const T& t, Vertex(T) v) {
		return attr(t,v)->id;
	}
	
	void id_set(T& t, Vertex(T) v, int i) {
		attr(t,v)->id=i;
	}

	int orientation(const T& t, Cell(T) f) {
		return attr(t,f)->o;
	}
	
	void orientation_set(T& t, Cell(T) f, int o) {
		attr(t,f)->o=o;
	}

}

namespace std {
	template <class T, dim_t k>
	ostream& operator<<(ostream& o, Simplex(T,k) s) {
		if(empty(s)) {
			o<<"()";
		} else {
			vgtl::array<Vertex(T),k+1> vs;
			vertices(t,s,vs);
			o<<"(";
      for(dim_t i=0; i<k; ++i) o<<id(t,vs[i])<<" ";
			o<<id(t,vs[k])<<")";
		}
		return o;
	}
}

int main(int argc, char * argv[]) {
	Vertex(T) v[6];
	for(int i=0; i<6; ++i) {
		v[i]=add(t);
		id_set(t,v[i],i);
	}
	{
		/*
		       3---1
					/ \ / \
				 5---2---0
		      \	 |  /
					 \ | /
					 	\|/
						 4
		 */
		vgtl::complex_builder<T> cb(t);
		orientation_set(t,add(cb,_(v[1],v[2],v[3])),-1);
		orientation_set(t,add(cb,_(v[0],v[2],v[4])),+1);
		orientation_set(t,add(cb,_(v[0],v[1],v[2])),+1);
		orientation_set(t,add(cb,_(v[2],v[3],v[5])),+1);
		orientation_set(t,add(cb,_(v[2],v[4],v[5])),-1);
	}
	cout<<endl;
	list<Edge(T)> re;
	radial_edges(t,v[2],re);	
	for(list<Edge(T)>::iterator i=re.begin(); i!=re.end(); ++i) {
		cout<<*i<<endl;
	}
}
