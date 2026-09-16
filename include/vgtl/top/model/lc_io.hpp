#ifndef VGTL_MODEL_LC_IO_HPP
#define VGTL_MODEL_LC_IO_HPP

/*! \file
 * \brief 
 */

#include <map>
#include <fstream>
#include <vgtl/top/model/lc.hpp>
#include <vgtl/top/incidence.hpp>

namespace vgtl {

	using std::ifstream;
	using std::ofstream;
	using std::ios;
	using std::endl;
	using std::map;

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Writes a ply file 
	template <template <dim_t> class Data>
	void write_ply(ofstream& out, const lc<2,Data>& t) {
		out<<"ply"<<endl; 
		out<<"format ascii 1.0"<<endl;
		int c=simplex_container<0,2,Data>(t).size();
		out<<"element vertex "<<c<<endl;
		out<<"property float x"<<endl;
		out<<"property float y"<<endl;
		out<<"property float z"<<endl; 
		c=simplex_container<2,2,Data>(t).size();
		out<<"element face "<<c<<endl;
		out<<"property list uchar int vertex_indices"<<endl;
		out<<"end_header"<<endl;
		map<simplex_descriptor<lc<2,Data>,0>,int> desc;
		simplex_iterator<lc<2,Data>,0> vi,vend;
		int i=0;
		for(simplices(t,vi,vend);vi!=vend;++vi,++i) {
			simplex_descriptor<lc<2,Data>,0> v=*vi;
			out<<euclidean_point(t,v)<<endl;
			desc[v]=i;
		}
		simplex_iterator<lc<2,Data>,2> ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			simplex_descriptor<lc<2,Data>,2> c=*ci;
			array<simplex_descriptor<lc<2,Data>,0>,3> vs;
			vertices(t,c,vs);
			out<<"3 ";
			out<<desc[vs[0]]<<" ";
			if(orientation(t,c)<0) {
				out<<desc[vs[2]]<<" ";
				out<<desc[vs[1]]<<endl;
			} else {
				out<<desc[vs[1]]<<" ";
				out<<desc[vs[2]]<<endl;
			}
		}
		
	}

	//! Writes a ply file with normals
	template <template <dim_t> class Data>
	void write_ply_with_normal(ofstream& out, const lc<2,Data>& t) {
		out<<"ply"<<endl; 
		out<<"format ascii 1.0"<<endl;
		int c=simplex_container<0,2,Data>(t).size();
		out<<"element vertex "<<c<<endl;
		out<<"property float x"<<endl;
		out<<"property float y"<<endl;
		out<<"property float z"<<endl; 
		out<<"property float nx"<<endl;
		out<<"property float ny"<<endl;
		out<<"property float nz"<<endl; 
		c=simplex_container<2,2,Data>(t).size();
		out<<"element face "<<c<<endl;
		out<<"property list uchar int vertex_indices"<<endl;
		out<<"end_header"<<endl;
		map<simplex_descriptor<lc<2,Data>,0>,int> desc;
		simplex_iterator<lc<2,Data>,0> vi,vend;
		int i=0;
		for(simplices(t,vi,vend);vi!=vend;++vi,++i) {
			simplex_descriptor<lc<2,Data>,0> v=*vi;
			out<<euclidean_point(t,v)<<endl;
			out<<normal(t,v)<<endl;
			desc[v]=i;
		}
		simplex_iterator<lc<2,Data>,2> ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			simplex_descriptor<lc<2,Data>,2> c=*ci;
			array<simplex_descriptor<lc<2,Data>,0>,3> vs;
			vertices(t,c,vs);
			out<<"3 ";
			out<<desc[vs[0]]<<" ";
			if(orientation(t,c)<0) {
				out<<desc[vs[2]]<<" ";
				out<<desc[vs[1]]<<endl;
			} else {
				out<<desc[vs[1]]<<" ";
				out<<desc[vs[2]]<<endl;
			}
		}
		
	}


	/*! @} */

}

#endif // VGTL_MODEL_LC_IO_HPP
