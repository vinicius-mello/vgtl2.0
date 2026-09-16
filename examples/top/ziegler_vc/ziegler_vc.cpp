#include <iostream>
#include <vgtl/top/model/vc.hpp>
#include <list>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/star.hpp>

using namespace std;
using vgtl::dim_t;

typedef vgtl::vc<3> T;

T t;

namespace vgtl {
	template <>
	struct extra_data<0> {
		int id;
	};

	int id(const T& t, Vertex(T) v) {
		return attr(t,v)->id;
	}
	
	void id_set(T& t, Vertex(T) v, int i) {
		attr(t,v)->id=i;
	}

}

namespace vgtl {
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
	Vertex(T) v[11];
	for(int i=1; i<=10; ++i) {
		v[i]=add(t);
		id_set(t,v[i],i);
	}
	{
		vgtl::complex_builder<T> cb(t);
		add(cb,_(v[1],v[2],v[3],v[4]));
		add(cb,_(v[1],v[2],v[5],v[6]));
		add(cb,_(v[2],v[3],v[6],v[7]));
		add(cb,_(v[3],v[4],v[7],v[8]));
		add(cb,_(v[1],v[4],v[5],v[8]));
		add(cb,_(v[1],v[5],v[6],v[9]));
		add(cb,_(v[1],v[2],v[6],v[9]));
		add(cb,_(v[1],v[2],v[4],v[9]));
		add(cb,_(v[1],v[4],v[8],v[9]));
		add(cb,_(v[1],v[5],v[8],v[9]));
		add(cb,_(v[2],v[5],v[6],v[10]));
		add(cb,_(v[2],v[6],v[7],v[10]));
		add(cb,_(v[2],v[3],v[7],v[10]));
		add(cb,_(v[1],v[2],v[3],v[10]));
		add(cb,_(v[1],v[2],v[5],v[10]));
		add(cb,_(v[3],v[6],v[7],v[8]));
		add(cb,_(v[2],v[3],v[4],v[8]));
		add(cb,_(v[2],v[3],v[6],v[8]));
		add(cb,_(v[4],v[5],v[7],v[8]));
		add(cb,_(v[1],v[3],v[4],v[7]));
		add(cb,_(v[1],v[4],v[5],v[7]));
	}
	cout<<endl;
	list<Cell(T)> st;
	star(t,v[1],back_inserter(st));	
	for(list<Cell(T)>::iterator i=st.begin(); i!=st.end(); ++i) {
		cout<<*i<<endl;
	}
}
