#ifndef VGTL_STATS_HPP
#define VGTL_STATS_HPP

/*! \file
 * \brief Reports simplicial complex statistics 
 */

#include <vgtl/top/sc.hpp>

namespace vgtl {
	
	using std::cout;
	using std::endl;

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	stats(const T& t, mt_tag);

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	stats(const T& t, mt_tag);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	stats(const T& t, asc_tag);

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	stats(const T& t, asc_tag);

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Shows some informations about a simplicial complex
	/*! The correct information is displayed accordingly
	 *  the <code>sc_traits\<T\>::sc_category</code> type.
	 */
	template <class T>
	void
	stats(const T& t) 
	{
		cout<<"stats:"<<endl;
		stats<T,0>(t, typename sc_traits<T>::sc_category());
	}

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	stats(const T& t, mt_tag)
	{
		Simplex_it(T,k) i,end;
		int cur=0;
		int all=0;
		for(simplices(t,i,end);i!=end;++i,++all) 
			if(is_current(t,*i)) ++cur;
		cout<<"--"<<k<<"--"<<endl;
		cout<<"current :"<<cur<<endl;
		cout<<"total   :"<<all<<endl;
		cout<<"ratio   :"<<((float)cur)/all<<endl;
		stats<T,k+1>(t,mt_tag());
	}

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	stats(const T& t, mt_tag)
	{
		Simplex_it(T,k) i,end;
		int cur=0;
		int all=0;
		for(simplices(t,i,end);i!=end;++i,++all) 
			if(is_current(t,*i)) ++cur;
		cout<<"--"<<k<<"--"<<endl;
		cout<<"current :"<<cur<<endl;
		cout<<"total   :"<<all<<endl;
		cout<<"ratio   :"<<((float)cur)/all<<endl;
	}

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	stats(const T& t, asc_tag)
	{
		Simplex_it(T,k) i,end;
		int all=0;
		for(simplices(t,i,end);i!=end;++i,++all) 
			;
		cout<<"--"<<k<<"--"<<endl;
		cout<<"total   :"<<all<<endl;
		stats<T,k+1>(t,asc_tag());
	}

	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	stats(const T& t, asc_tag)
	{
		Simplex_it(T,k) i,end;
		int all=0;
		for(simplices(t,i,end);i!=end;++i,++all) 
			;
		cout<<"--"<<k<<"--"<<endl;
		cout<<"total   :"<<all<<endl;
	}

	/*! @} */

}

#endif // VGTL_STATS_HPP
