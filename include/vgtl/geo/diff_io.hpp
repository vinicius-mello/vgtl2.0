#ifndef VGTL_DIFF_IO_HPP
#define VGTL_DIFF_IO_HPP

/*! \file
 * \brief Mesh IO
 */

#include <vgtl/geo/mesh_io.hpp>
#include <vgtl/top/simplex_io.hpp>

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

	template <class T, dim_t k>
	void 
	read_diff(ifstream& in, T& t, vector<Vertex(T)>& vv) {
		Simplex_it(T,k) si,send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			read_simplex(in,t,s,vv);
			point<k+1,double> d;
			in>>d;
			diffpars_set(t,s,d);
		}
	}
	
	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	read_diffs(ifstream& in, T& t, vector<Vertex(T)>& vv);
					
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	read_diffs(ifstream& in, T& t, vector<Vertex(T)>& vv) {
		read_diff<T,k>(in,t,vv);
	}
	
	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	read_diffs(ifstream& in, T& t, vector<Vertex(T)>& vv) {
		read_diff<T,k>(in,t,vv);
		read_diffs<T,k+1>(in,t,vv);
	}
	
	template <class T>
	void
	read_diffs(ifstream& in, T& t, vector<Vertex(T)>& vv) {
		read_diffs<T,1>(in,t,vv);
	}
	
	template <class T, class Apply>
	void 
	read_mesh_diffs(ifstream& in, T& t, Apply& app) {
		vector<Vertex(T)> vv;
		read_mesh_iso(in,t,app,vv);
		read_diffs<T,1>(in,t,vv);
	}

	template <class T, dim_t k>
	void 
	write_diff(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		Simplex_it(T,k) si,send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			write_simplex(out,t,s,mv);
			out<<diffpars(t,s)<<endl;
		}
	}
	
	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	write_diffs(ofstream& out, const T& t, map<Vertex(T),int>& mv);
					
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	write_diffs(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		write_diff<T,k>(out,t,mv);
	}

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	write_diffs(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		write_diff<T,k>(out,t,mv);
		write_diffs<T,k+1>(out,t,mv);
	}
	
	template <class T>
	void 
	write_diffs(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		write_diffs<T,1>(out,t,mv);
	}
	
	template <class T>
	void 
	write_mesh_diffs(ofstream& out, T& t) {
		map<Vertex(T),int> mv;
		write_mesh_iso(out,t,mv);
		write_diffs<T,1>(out,t,mv);
	}


}

#endif // VGTL_DIFF_IO_HPP
