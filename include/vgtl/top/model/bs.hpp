#ifndef VGTL_BS2_HPP
#define VGTL_BS2_HPP

/*! \file 
 * \brief ASC simplex storage class
 */

#include <vgtl/top/sc.hpp>
#include <vgtl/top/model/extra_data.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	template <class T, dim_t k> 
	struct bs_faces {
		Simplex(T,k-1) fop[k+1];
	};
				
	template <class T> 
	struct bs_faces<T,0> {
	};
				
	template <class T, dim_t k, bool t=(k<(Dim(T)-1))> 
	struct bs_up {
		Simplex(T,k+1) up;
	};

	template <class T, dim_t k> 
	struct bs_up<T,k,false> {
	};
	
	template <class T, dim_t k, bool t=(k!=Dim(T))> 
	struct bs_cell {
	};
	
	template <class T, dim_t k> 
	struct bs_cell<T,k,false> {
		bool cur : 1;
		bool ori : 1;
		mark_type mk : 2;
		short lvl: 16;
    bs_cell() : cur(true), ori(false), mk(white_mark), lvl(0) {}
	};
				
	template <class T, dim_t k, bool t=(k!=(Dim(T)-1))> 
	struct bs_facet {
	};
				
	template <class T, dim_t k> 
	struct bs_facet<T,k,false> {
		Cell(T) cov[2];
		bs_facet() {};
	};

	template <class T, dim_t k,
	          template <dim_t> class Data=extra_data>
	struct bs : Data<k>, bs_faces<T,k>,
			bs_cell<T,k>, bs_facet<T,k>, bs_up<T,k> {
	};

	/*! @} */

}

#endif // VGTL_BS2_HPP
