#include <iostream>
#include <string>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/top/subdivide_until.hpp>
#include <vgtl/top/cube.hpp>
#include <vgtl/top/radial.hpp>
#include <vgtl/top/initial.hpp>

#define N 2

using namespace std;
using namespace vgtl;

typedef vgtl::nmt<N> T;

T t;

namespace vgtl {
	template <>
	struct extra_data<0> {
		string id;
	};

	string id(const T& t, Vertex(T) v) {
		return attr(t,v)->id;
	}
	
	void id_set(T& t, Vertex(T) v, string i) {
		attr(t,v)->id=i;
	}

	struct new_simplices : public do_nothing {
		Vertex(T) add_edge_vertex(T& t, Edge(T) e) {
			Vertex(T) v=add(t);
			vgtl::array<Vertex(T),2> vs;
			vertices(t,e,vs);
			string ii="["+id(t,vs[0])+","+id(t,vs[1])+"]";
			id_set(t,v,ii);
			return v;
		}
	};
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
			for(int i=0; i<k; ++i) o<<id(t,vs[i])<<" ";
			o<<id(t,vs[k])<<")";
		}
		return o;
	}
}

void cvshow(const T& t) {
	Cell_it(T) cvi, cvend;
	for(simplices(t,cvi,cvend);cvi!=cvend;++cvi) {
		Cell(T) cv=*cvi;
		cout<<(orientation(t,cv)>0 ? "+" : "-")<<level(t,cv)<<cv;
	}
	cout<<endl;
}

Cell(T) current_cell(const T& t) {
	Cell_it(T) cvi, cvend;
	simplices(t,cvi,cvend);
	Cell(T) cv;
	do {
		cv=*cvi;
		if(is_current(t,cv)) return cv;
		++cvi;
	} while(1);
}

int main(int argc, char * argv[]) {
	vgtl::array<Vertex(T),4> vs;
	char idbuf[2];
	for(int i=0;i<4;++i) {
		vs[i]=add(t);
		idbuf[0]='a'+i;idbuf[1]=0;
		id_set(t,vs[i],idbuf);
	}
	Cell(T) cv;
	{
		vgtl::complex_builder<T> cb(t);
		add_cube(cb,vs);
	}
	cvshow(t);
	new_simplices ns;

	maubach_subdivide_to_level(t,3,ns);

	cvshow(t);
	Vertex_it(T) vi, vend;
  for(simplices(t,vi,vend);vi!=vend;++vi) {
	  list<Edge(T)> re;
		cout<<*vi<<"-->"<<endl;
  	radial_edges(t,*vi,re);
  	for(list<Edge(T)>::iterator i=re.begin(); i!=re.end(); ++i) {
   		cout<<*i<<endl;
  	}
	}
	cout<<"--------------------"<<endl;
	initial(t);
  for(simplices(t,vi,vend);vi!=vend;++vi) {
		if(!is_current(t,*vi)) continue;
	  list<Edge(T)> re;
		cout<<*vi<<"-->"<<endl;
  	radial_edges(t,*vi,re);
  	for(list<Edge(T)>::iterator i=re.begin(); i!=re.end(); ++i) {
   		cout<<*i<<endl;
  	}
	}

}
