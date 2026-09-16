#ifndef VGTL_MODEL_VC_IO_HPP
#define VGTL_MODEL_VC_IO_HPP

/*! \file
 * \brief 
 */

#include <vgtl/top/model/vc.hpp>
#include <vgtl/top/incidence.hpp>
#include <fstream>

namespace vgtl {

	using std::ifstream;
	using std::ofstream;
	using std::ios;
	using std::endl;

	/*! \addtogroup top 
	 * @{
	 */
	
	template <template <dim_t> class Data>
	void write_ply(ofstream& out, const vc<2,Data>& t) {
		out<<"ply"<<endl; 
		out<<"format ascii 1.0"<<endl;
		int i=simplex_container<0,2,Data>(t).size();
		out<<"element vertex "<<i<<endl;
		out<<"property float x"<<endl;
		out<<"property float y"<<endl;
		out<<"property float z"<<endl; 
		simplex_iterator<vc<2,Data>,2> ci,cend;
		i=0;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			simplex_descriptor<vc<2,Data>,2> c=*ci;
			if(is_current(t,c)) i++;
		}
		out<<"element face "<<i<<endl;
		out<<"property list uchar int vertex_indices"<<endl;
		out<<"end_header"<<endl;
		simplex_iterator<vc<2,Data>,0> vi,vend;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			simplex_descriptor<vc<2,Data>,0> v=*vi;
			out<<euclidean_point(t,v)<<endl;
		}
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			simplex_descriptor<vc<2,Data>,2> c=*ci;
			if(!is_current(t,c)) continue;
			array<simplex_descriptor<vc<2,Data>,0>,3> vs;
			vertices(t,c,vs);
			out<<"3 ";
			out<<vs[0].desc<<" ";
			if(orientation(t,c)<0) {
				out<<vs[2].desc<<" ";
				out<<vs[1].desc<<endl;
			} else {
				out<<vs[1].desc<<" ";
				out<<vs[2].desc<<endl;
			}
		}
		
	}

	/*! @} */

}

#endif // VGTL_MODEL_VC_IO_HPP
