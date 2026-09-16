#include <iostream>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/utl/pair_tie.hpp>
//#include <vgtl/top/complex_io.hpp>
#include <vgtl/top/iso/iso.hpp>
#include <vgtl/top/model/lc.hpp>

using namespace std;
using namespace vgtl;

typedef vgtl::lc<3> T;

T t;

namespace vgtl {
	template <>
	struct extra_data<0> {
		int id;
		int signal;
	};

	int id(const T& t, Vertex(T) v) {
		return attr(t,v)->id;
	}
	
	void id_set(T& t, Vertex(T) v, int i) {
		attr(t,v)->id=i;
	}

	int signal(const T& t, Vertex(T) v) {
		return attr(t,v)->signal;
	}
	
	void signal_set(T& t, Vertex(T) v, int i) {
		attr(t,v)->signal=i;
	}

}

namespace std {
	template <class T, dim_t k>
	ostream& operator<<(ostream& o,
											Simplex(T,k) s) {
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
	Vertex(T) v[4];
	for(int i=0; i<4; ++i) {
		v[i]=add(t);
		id_set(t,v[i],i);
	}
	Cell(T) s;
	{
		vgtl::complex_builder<T> cb(t);
		s=add(cb,_(v[0],v[1],v[2],v[3]));
	}
	signal_set(t,v[0],+1);
	signal_set(t,v[1],-1);
	signal_set(t,v[2],-1);
	signal_set(t,v[3],+1);
	list<vgtl::array<Edge(T),3> > ls;
	iso::triangulate(t,s,ls);
	for(list<vgtl::array<Edge(T),3> >::iterator i=ls.begin(); i!=ls.end(); ++i) {
		cout<<*i<<endl;
		cout<<iso::orientation(t,s,*i)<<endl;
	}
}
