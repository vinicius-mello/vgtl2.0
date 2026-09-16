#include <iostream>
#include <vgtl/top/model/lc.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/star.hpp>
#include <vgtl/top/decompose.hpp>
#include <vgtl/top/unlink.hpp>
#include <vgtl/top/del.hpp>
#include <vgtl/alg/point.hpp>
#include <vgtl/geo/cube.hpp>

using namespace std;
using namespace vgtl;

#define NN 3

typedef vgtl::lc<NN> T;

T t;

namespace vgtl {
	template <>
	struct extra_data<0> {
		point<NN,float> p;
	};

	template <>
	struct extra_data<NN> {
		vgtl::array<unsigned int,NN> p;
	};

	inline
	point<NN,float> euclidean_point(const T& t, Vertex(T) v) {
		return attr(t,v)->p;
	}
	
	inline
	void euclidean_point_set(T& t, Vertex(T) v, const vgtl::array<int,NN>& a) {
		point<NN,float> p(a);
		attr(t,v)->p=p;
	}

	inline
	vgtl::array<unsigned int,NN> perm(const T& t, Cell(T) cv) {
		return attr(t,cv)->p;
	}
	
	inline
	void perm_set(T& t, Cell(T) cv, const vgtl::array<unsigned int,NN>& p) {
		attr(t,cv)->p=p;
	}

	int id(const T& t, Vertex(T) v) {
		point<NN,float> p=euclidean_point(t,v);
		int r=0;
		for(int i=0;i<NN;++i) 
			r+=(int)((1<<i)*p[i]);
		return r;
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
	do_nothing dn;
	Edge(T) e=set_cube(t,dn);
	permutations_traverser<NN> p;
	Cell_it(T) cvit, cvend;
	simplices(t,cvit,cvend);
	do {
		perm_set(t,*cvit,*p);
		++cvit;
	} while(++p);
	list<Cell(T)> st;
	star(t,e,back_inserter(st));	
	for(list<Cell(T)>::iterator i=st.begin(); i!=st.end(); ++i) {
		Cell(T) cur=*i;
		cout<<opposite(t,cur,e)<<" : "<<perm(t,cur)<<endl;
	}
}
