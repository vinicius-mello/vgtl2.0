#ifndef VGTL_MS2_HPP
#define VGTL_MS2_HPP

/*! \file
 * \brief BMT simplex storage class
 */

#include <vgtl/top/sc.hpp>
#include <vgtl/top/model/bs.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T, dim_t k, bool t=(k!=Dim(T))> 
	struct ms_hier {
	};
				
	template <class T, dim_t k> 
	struct ms_hier<T,k,false> {
		Cell(T) child[SplitDim(T)+1];
		Cell(T) parent;
		ms_hier() {}
	};
				
	template <class T, dim_t k, bool t=(k!=SplitDim(T))> 
	struct ms_split {
	};
				
	template <class T, dim_t k> 
	struct ms_split<T,k,false> {
		Vertex(T) split;
		ms_split() {}
	};
				
	template <class T, dim_t k> 
	struct ms_vertex {
	};
				
	template <class T> 
	struct ms_vertex<T,0> {
		SplitSimplex(T) split;
		ms_vertex() {}
	};

	template <class T, dim_t k,
	          template <dim_t> class Data=extra_data>
	struct ms : bs<T,k,Data>, ms_hier<T,k>, ms_split<T,k>, ms_vertex<T,k> {
	};
	
	/*! @} */

}

#endif // VGTL_MS2_HPP
