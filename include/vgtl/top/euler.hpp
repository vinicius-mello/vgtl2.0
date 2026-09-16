#ifndef VGTL_EULER_HPP
#define VGTL_EULER_HPP

/*! \file
 * \brief Computes the Euler characteristic
 */

#include <vgtl/top/sc.hpp>

namespace vgtl {
	
	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, mt_tag); 

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, mt_tag);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, asc_tag);

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, asc_tag);

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Computes the Euler characteristic
	template <class T>
	int
	euler_characteristic(const T& t) 
	{
		return euler_characteristic<T,0>(t, typename sc_traits<T>::sc_category());
	}

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, mt_tag)
	{
		Simplex_it(T,k) i,end;
		int j=0;
		for(simplices(t,i,end);i!=end;++i) 
			if(is_current(t,*i)) ++j;
		return ((k%2)?-1:1)*j+euler_characteristic<T,k+1>(t,mt_tag());
	}

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, mt_tag)
	{
		Simplex_it(T,k) i,end;
		int j=0;
		for(simplices(t,i,end);i!=end;++i) 
			if(is_current(t,*i)) ++j;
		return ((k%2)?-1:1)*j;
	}

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, asc_tag)
	{
		Simplex_it(T,k) i,end;
		int j=0;
		for(simplices(t,i,end);i!=end;++i) 
			++j;
		return ((k%2)?-1:1)*j+euler_characteristic<T,k+1>(t,asc_tag());
	}

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,int>::type
	euler_characteristic(const T& t, asc_tag)
	{
		Simplex_it(T,k) i,end;
		int j=0;
		for(simplices(t,i,end);i!=end;++i) 
			++j;
		return ((k%2)?-1:1)*j;
	}

	/*! @} */

}

#endif // VGTL_EULER_HPP
