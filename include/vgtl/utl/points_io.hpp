#ifndef VGTL_POINTS_IO_HPP
#define VGTL_POINTS_IO_HPP

/*! \file
 * \brief Points IO operations
 */

#include <vgtl/alg/point.hpp>
#include <iostream>
#include <vector>

namespace vgtl {

	using std::vector;
	using std::istream;
	using std::ostream;
	using std::endl;
	
	/*! \addtogroup utl 
	 * @{
	 */

	//! Reads an array from a stream
  template <dim_t k, class Scalar>
	void read(istream& in, vector<point<k,Scalar> >& vp) {
		int s;
		in>>s;
		for(int i=0;i<s;++i) {
			point<k,Scalar> p;
			for(int j=0;j<k;++j) in>>p[j];
			vp.push_back(p);
		} 
	}
	
	//! Reads an array from a stream
  template <dim_t k, class Scalar>
	void write(ostream& out, vector<point<k,Scalar> >& vp) {
		int s=vp.size();
		out<<s<<endl;
		for(int i=0;i<s;++i) {
			point<k,Scalar> p=vp[i];
			for(int j=0;j<k;++j) out<<p[j]<<" ";
			out<<endl;
		} 
	}

	template <class Scalar>
	void write_ply(ostream& out, vector<point<3,Scalar> >& vp) {
		int s=vp.size();
		out<<"ply"<<endl;
		out<<"format ascii 1.0"<<endl;
		out<<"element vertex ";
		out<<s<<endl;
		out<<"property float x"<<endl;
		out<<"property float y"<<endl;
		out<<"property float z"<<endl;
		out<<"end_header"<<endl;
		for(int i=0;i<s;++i) {
			point<3,Scalar> p=vp[i];
			for(int j=0;j<3;++j) out<<p[j]<<" ";
			out<<endl;
		} 
	}

	/*! @} */

}

#endif //VGTL_POINTS_IO_HPP
