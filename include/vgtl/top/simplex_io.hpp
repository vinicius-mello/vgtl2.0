#ifndef VGTL_SIMPLEX_IO_HPP
#define VGTL_SIMPLEX_IO_HPP

/*! \file
 * \brief Simplex IO operations
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <vgtl/top/vertices_inv.hpp>

namespace vgtl {
	
  using std::vector;
	using std::map;
	using std::ifstream;
	using std::ofstream;

	template <class T, dim_t k>
	typename enable_if<vertex_c<T,k>,void>::type
	read_simplex(ifstream& in, T& t, Simplex(T,k)& s, vector<Vertex(T)>& vv) {
		int i;
		in>>i;
		s=vv[i];
	}

	template <class T, dim_t k>
	typename disable_if<vertex_c<T,k>,void>::type
	read_simplex(ifstream& in, T& t, Simplex(T,k)& s, vector<Vertex(T)>& vv) {
		array<Vertex(T),k+1> vs;
		for(int i=0;i<=k;++i) {
			int j;
			in>>j;
			vs[i]=vv[j];
		}
		s=vertices_inv(t,vs);
	}

	template <class T, dim_t k>
	typename enable_if<vertex_c<T,k>,void>::type
	write_simplex(ofstream& out, const T& t, Simplex(T,k) s, map<Vertex(T),int>& mv) {
		out<<mv[s]<<" ";
	}

	template <class T, dim_t k>
	typename disable_if<vertex_c<T,k>,void>::type
	write_simplex(ofstream& out, const T& t, Simplex(T,k) s, map<Vertex(T),int>& mv) {
		array<Vertex(T),k+1> vs;
		vertices(t,s,vs);
		for(int i=0;i<=k;++i) {
			out<<mv[vs[i]]<<" ";
		}
	}

}

#endif // VGTL_SIMPLEX_IO_HPP
