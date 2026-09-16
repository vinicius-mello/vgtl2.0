#ifndef VGTL_LINK_HPP
#define VGTL_LINK_HPP

/*! \file
 * \brief Generic functions to correct the links of a simplex
 */

#include <utility>
#include <vgtl/top/complex_buffer.hpp>

namespace vgtl {
	
	using std::set;
	using std::pair;

	template <class T>
	inline
	void
	link(complex_buffer<T>& cb);

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	link(const complex_buffer<T> &cb);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	link(const complex_buffer<T> &cb);

	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T, dim_t k>
	void
	link(T& t,
			 Simplex(T,k) s) {
		for(dim_t i=0; i<=k; ++i)
			up_simplex_set(t,face_op(t,s,i),s);
	}

 	//!Corrects the links of a simplex
	template <class T>
	inline
	void
	link(complex_buffer<T>& cb) {
		link<T,1>(cb);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	link(const complex_buffer<T> &cb)
	{
		for(typename set<Simplex(T,k)>::const_iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			Cell(T) cv=*i;
      for(dim_t j=0;j<=Dim(T);++j) {
				Facet(T) ce=face_op(cb.t,cv,j);
				pair<Cell(T),Cell(T)> q=cells(cb.t,ce);
				if(exists(cb,ce)) {
					if(q.first==q.second) { //Begin nice trick
						q.first=q.second=cv;
					} else if(empty(q.second)) {
						q.second=cv;
					} else {
						q.first=cv;
						empty_set(q.second);
					}	// End nice trick
					cells_set(cb.t,ce,q.first,q.second);
				} else {
					if(empty(q.first)&&empty(q.second)) cells_set(cb.t,ce,cv,cv);
					else if(q.first==q.second) cells_set(cb.t,ce,q.first,cv);
					// else throw exception
				}
			}
		}
	}

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	link(const complex_buffer<T> &cb)
	{
		for(typename set<Simplex(T,k)>::const_iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			link(cb.t,*i);
		}
		link<T,k+1>(cb);
	}
	
	/*! @} */

}

#endif // VGTL_LINK_HPP
