#ifndef VGTL_COMPLEX_BUFFER_IO_HPP
#define VGTL_COMPLEX_BUFFER_IO_HPP

/*! \file
 *  \brief complex_buffer IO operations
 */

#include <vgtl/top/complex_buffer.hpp>

namespace vgtl {
	
	using std::cout;

	/*! \addtogroup top 
	 * @{
	 */

	template <class T, dim_t k>
	typename disable_if<vertex_c<T,k>,void>::type
	show(const complex_buffer<T> &cb)
	{
		cout<<"--"<<k<<"--"<<endl;
	  for(typename set<Simplex(T,k)>::iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			cout<<*i<<endl;
		}
		show<T,k-1>(cb);
	}

	template <class T, dim_t k>
	typename enable_if<vertex_c<T,k>,void>::type
	show(const complex_buffer<T> &cb)
	{
		cout<<"--"<<k<<"--"<<endl;
	  for(typename set<Simplex(T,k)>::iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			cout<<*i<<endl;
		}
		std::cout<<"-----"<<endl;
	}

	//! Shows the contents of a complex buffer
	template <class T>
	void
	show(const complex_buffer<T> &cb) 
	{
		show<T,Dim(T)>(cb);
	}

	/*! @} */
	
}

#endif // VGTL_COMPLEX_BUFFER_IO_HPP
