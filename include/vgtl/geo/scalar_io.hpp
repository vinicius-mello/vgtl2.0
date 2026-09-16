#ifndef VGTL_SCALAR_IO_HPP
#define VGTL_SCALAR_IO_HPP

/*! \file
 * \brief Mesh IO
 */

#include <fstream>
#include <vgtl/top/sc.hpp>

namespace vgtl {
	
	using std::vector;
	using std::map;
	using std::ifstream;
	using std::ofstream;
	using std::ios;
	using std::endl;

	/*! \addtogroup top 
	 * @{
	 */

	template <class T>
	void 
	read_scalar(ifstream& in, T& t, vector<Vertex(T)>& vv) {
		Vertex_it(T) si,send;
		for(simplices(t,si,send);si!=send;++si) {
			Vertex(T) s=*si;
			double d;
			in>>d;
			scalar_value_set(t,s,d);
			vv.push_back(s);
		}
	}

	template <class T>
	void 
	write_scalar(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		Vertex_it(T) si,send;
		int i=0;
		for(simplices(t,si,send);si!=send;++si,++i) {
			Vertex(T) s=*si;
			out<<scalar_value(t,s)<<endl;
			mv[s]=i;
		}
	}
	
	template <class T>
	void 
	read_scalar(ifstream& in, T& t) {
		Vertex_it(T) si,send;
		for(simplices(t,si,send);si!=send;++si) {
			Vertex(T) s=*si;
			double d;
			in>>d;
			scalar_value_set(t,s,d);
		}
	}

	template <class T>
	void 
	write_scalar(ofstream& out, const T& t) {
		Vertex_it(T) si,send;
		int i=0;
		for(simplices(t,si,send);si!=send;++si,++i) {
			Vertex(T) s=*si;
			out<<scalar_value(t,s)<<endl;
		}
	}
	
}

#endif // VGTL_SCALAR_IO_HPP
