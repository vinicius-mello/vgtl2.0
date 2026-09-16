#ifndef VGTL_GEO_TOPOLOGICAL_SORT_HPP
#define VGTL_GEO_TOPOLOGICAL_SORT_HPP

#include <utility>
#include <vgtl/alg/vec.hpp>
#include <vgtl/top/pm.hpp>

/*! \file
 * \brief sort functions
 */

namespace vgtl {
	
	using std::pair;

	/*! \addtogroup geo 
	 * @{
	 */

	template <class T> 
	class normal_arrow {
		vec<Dim(T),double> dir;
		public:
		normal_arrow(const vec<Dim(T),double>& _dir) : dir(_dir) {}
		Cell(T) operator()(const T& t, Cell(T) c, int i) const {
			Facet(T) f=face_op(t,c,i);
			Cell(T) ca=adjacent(t,c,f);
			if(c==ca) return c;
			double o=((i%2)?-1:1)*orientation(t,c);
			if(dot(dir,o*normal(t,f))>0) return ca;
			else return c;
		}
	};
				
	template <class T, class OutputIt, class Arrow>
	void
	ts_dfs(const T& t, Cell(T) c, OutputIt out, Arrow& arrow) {
		mark_set(t,c,gray_mark);
		for(int i=0;i<=Dim(T);++i) {
			Cell(T) ca=arrow(t,c,i);
			if(ca!=c) {
				mark_type mk=mark(t,ca);
				if(mk==white_mark) {
					ts_dfs(t,ca,out,arrow);
				} else if(mk==gray_mark) {
					//throw exception
				}
			}
		}
		mark_set(t,c,black_mark);
		*out=c;++out;
	}
	
	template <class T, class OutputIt, class Arrow>
	void
	topological_sort(const T& t, OutputIt out, Arrow& arrow) {
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			if(mark(t,c)==white_mark) {
				ts_dfs(t,c,out,arrow);
			}
		}
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			mark_set(t,c,white_mark);
		}
	}

	/*! @} */

}

#endif // VGTL_GEO_TOPOLOGICAL_SORT_HPP
