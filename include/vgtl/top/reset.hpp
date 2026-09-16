#ifndef VGTL_RESET_HPP
#define VGTL_RESET_HPP

/*! \file
 *  \brief reset subdivision
 */

#include <vgtl/top/sc.hpp>

namespace vgtl {
	

	/*! \addtogroup top 
	 * @{
	 */

	template <class T, dim_t k>
	typename disable_if<facet_c<T,k>,void>::type
	reset_subdivision(T& t);
	
	template <class T, dim_t k>
	typename enable_if<facet_c<T,k>,void>::type
	reset_subdivision(T& t);
	
	template <class T, dim_t k>
	typename disable_if<facet_c<T,k>,void>::type
	reset_subdivision(T& t) {
		Simplex_it(T,k) si,send;
		list<Simplex(T,k)> ls;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) ls.push_back(s);
		}
		for(typename list<Simplex(T,k)>::iterator i=ls.begin();i!=ls.end();++i) {
			del(t,*i);
		}
		reset_subdivision<T,k+1>(t);
	}
	
	template <class T, dim_t k>
	typename enable_if<facet_c<T,k>,void>::type
	reset_subdivision(T& t) {
		Simplex_it(T,k) si,send;
		list<Simplex(T,k)> ls;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) ls.push_back(s);
		}
		for(typename list<Simplex(T,k)>::iterator i=ls.begin();i!=ls.end();++i) {
			del(t,*i);
		}
	}
	
	template <class T>
	void reset_subdivision(T& t) {
		reset_subdivision<T,0>(t);
		Cell_it(T) ci,cend;
		list<Cell(T)> lc;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(is_current(t,c)) {
				Cell(T) empt;
				empty_set(empt);
				parent_set(t,c,empt);
				for(int i=0;i<=SplitDim(T);++i) child_set(t,c,i,empt);
				level_set(t,c,0);
			} else lc.push_back(c);
		}
		for(typename list<Cell(T)>::iterator i=lc.begin();i!=lc.end();++i) {
			del(t,*i);
		}
	}

	/*! @} */
	
}

#endif // VGTL_RESET_HPP
