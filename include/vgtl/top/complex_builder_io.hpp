#ifndef VGTL_COMPLEX_BUILDER_IO_HPP
#define VGTL_COMPLEX_BUILDER_IO_HPP

/*! \file
 * \brief complex_builder IO operations
 */

#include <vgtl/top/complex_builder.hpp>

namespace vgtl {
	
	using std::map;

	/*! \addtogroup top 
	 * @{
	 */

	//! Shows the contents of a complex builder
	template <class T>
	void
	show(const complex_builder<T> &cb) 
	{
		show<T,Dim(T)>(cb);
	}

	template <class T, dim_t k>
	typename disable_if<edge_c<T,k>,void>::type
	show(const complex_builder<T> &cb)
	{
		typedef map<array<Vertex(T), k+1>, Simplex(T,k)> maps;
		std::cout<<"--"<<k<<"--"<<endl;
	  for(typename maps::const_iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			std::cout<<i->second<<endl;
		}
		show<T,k-1>(cb);
	}

	template <class T, dim_t k>
	typename enable_if<edge_c<T,k>,void>::type
	show(const complex_builder<T> &cb)
	{
		typedef map<array<Vertex(T), k+1>, Simplex(T,k)> maps;
		std::cout<<"--"<<k<<"--"<<endl;
	  for(typename maps::const_iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			std::cout<<i->second<<endl;
		}
		std::cout<<"-----"<<endl;
	}

	/*! @} */

}

#endif // VGTL_COMPLEX_BUILDER_IO_HPP
