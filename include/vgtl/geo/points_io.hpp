#ifndef VGTL_GEO_POINTS_IO_HPP
#define VGTL_GEO_POINTS_IO_HPP

/*! \file
 * \brief Mesh IO
 */

#include <algorithm>
#include <vgtl/alg/point.hpp>
#include <vgtl/top/simplex_io.hpp>

namespace vgtl {
	
	using std::vector;
	using std::map;
	using std::ifstream;
	using std::ofstream;
	using std::min;
	using std::max;
	using std::endl;

	/*! \addtogroup top 
	 * @{
	 */
	template <class T>
	void 
	read_points(ifstream& in, T& t, vector<Vertex(T)>& vv) {
		typename vector<point<Dim(T),double> >::iterator nu=points(t).begin();
		typename vector<point<Dim(T),double> >::iterator nv=points(t).begin()+points(t).size();
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			points_begin_set(t,c,nv);
			points_end_set(t,c,nu);
		}
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			read_simplex(in,t,c,vv);
			int ii,iend;
			in>>ii;
			in>>iend;
			points_begin_set(t,c,points(t).begin()+ii);
			points_end_set(t,c,points(t).begin()+iend);
		}
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(has_children(t,c)) continue;
			Cell(T) par=c;
			while(has_parent(t,par)) {
				par=parent(t,par);
				points_begin_set(t,par,min(points_begin(t,par),points_begin(t,c)));
				points_end_set(t,par,max(points_end(t,par),points_end(t,c)));
			}
		}
	}

	template <class T>
	void 
	write_points(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			write_simplex(out,t,c,mv);
			out<<(points_begin(t,c)-points(t).begin())<<" ";
			out<<(points_end(t,c)-points(t).begin())<<endl;
		}
	}

}

#endif // VGTL_GEO_POINTS_IO_HPP
