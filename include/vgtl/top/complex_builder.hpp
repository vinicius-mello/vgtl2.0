#ifndef VGTL_COMPLEX_BUILDER_HPP
#define VGTL_COMPLEX_BUILDER_HPP

/*! \file
 * \brief Maps vertex sets to simplices
 */

#include <map>
#include <vgtl/utl/array.hpp>
#include <vgtl/top/sc.hpp>
#include <vgtl/top/incidence.hpp>
#include <vgtl/top/do_nothing.hpp>
#include <vgtl/top/complex_buffer.hpp>

namespace vgtl {
	
	using std::map;
	
	/*! \addtogroup top 
	 * @{
	 */

	template <class T, dim_t k>
	struct complex_builder_rep : complex_builder_rep<T,k-1> {
		map<array<Vertex(T),k+1>,Simplex(T,k)> m;
		complex_builder_rep<T,k>(T& _t) : complex_builder_rep<T,k-1>(_t) {}
	};

	template <class T>
	struct complex_builder_rep<T,0> {
		T& t;
		complex_builder_rep<T,0>(T& _t) : t(_t) {}
	};


	//! Maps vertex sets to simplices
	/*! A complex builder its just a set of <code>sc_traits\<T\>::dim</code>
	 * maps between vertex sets and simplices, one for each dimension above zero.
	 * It works like a heterogeneus map between vertex sets and
	 * <code>simplex_descriptor</code>.*/
	template <class T>
	struct complex_builder : complex_builder_rep<T,Dim(T)> {
		complex_builder(T& _t) : complex_builder_rep<T,Dim(T)>(_t) {}
		template <dim_t k>
		typename map<array<Vertex(T),k+1>,Simplex(T,k)>::const_iterator
		begin() const {
			return (*this).complex_builder_rep<T,k>::m.begin(); 
		}
		template <dim_t k>
		typename map<array<Vertex(T),k+1>,Simplex(T,k)>::const_iterator
		end() const {
			return (*this).complex_builder_rep<T,k>::m.end(); 
		}
	};

	//! Returns the map of \c k dimensional simplices
	template <dim_t k, class T>
	map<array<Vertex(T),k+1>,Simplex(T,k)>&
	builder_map(complex_builder<T> &cb) {
	  return cb.complex_builder_rep<T,k>::m;
	}

	//! Returns the map of \c k dimensional simplices
	template <dim_t k, class T>
	const map<array<Vertex(T),k+1>,Simplex(T,k)>&
	builder_map(const complex_builder<T> &cb) {
	  return cb.complex_builder_rep<T,k>::m;
	}

	//! Tests if a vertex set belongs to a complex builder
	template <class T, dim_t kk>
	bool
	exists(const complex_builder<T> &cb, 
		     const array<Vertex(T),kk>& key)
	{
	  if(builder_map<kk-1>(cb).count(key)==0) return false;
		return true;
	}

	//! Returns the simplex corresponding to a vertex set 
	template <class T, dim_t kk>
	Simplex(T,kk-1) 
	get(const complex_builder<T>& cb, 
		  const array<Vertex(T),kk>& key)
	{
	  return (*builder_map<kk-1>(cb).find(key)).second;
	}

	//! Associates a simplex to a vertex set 
	template <class T, dim_t kk>
	void
	put(complex_builder<T>& cb, 
		  const array<Vertex(T),kk>& key,
			Simplex(T,kk-1) s)
	{
	  builder_map<kk-1>(cb)[key]=s;
	}
	
	template <class T, dim_t k>
	void
	put(complex_builder<T>& cb, 
			Simplex(T,k) s)
	{
		array<Vertex(T),k+1> key;
		vertices(cb.t,s,key);
	  builder_map<k>(cb)[key]=s;
	}
	
	template <class T, class Apply, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	apply(complex_builder<T> &cb,
				Apply& app) {
		typedef map<array<Vertex(T),k+1>,Simplex(T,k)> maps;
	  for(typename maps::const_iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			app.apply(cb.t,i->second);
		}
	}
	
	template <class T, class Apply, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	apply(complex_builder<T> &cb,
				Apply& app) {
		typedef map<array<Vertex(T),k+1>,Simplex(T,k)> maps;
	  for(typename maps::const_iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			app.apply(cb.t,i->second);
		}
		apply<T,Apply,k+1>(cb,app);
	}
	
	//! Applies a function to all simplices in a complex builder
	template <class T, class Apply>
	inline
	void
	apply(complex_builder<T> &cb, Apply& app) {
		apply<T,Apply,1>(cb,app);
	}
	
	template <class T>
	struct build_apply {
		complex_builder<T>& cb;
		build_apply(complex_builder<T>& _cb) : cb(_cb) {}
		template <dim_t k>
		void
		apply(T& t, Simplex(T,k) s) {
			put(cb,s);
		}
		void
		apply(T& t, Simplex(T,0) s) {
		}
	};

	template <class T>
	void buffer_to_build(complex_buffer<T>& cb1,
		complex_builder<T>& cb2) {
		build_apply<T> ba(cb2);
		apply(cb1,ba);
	}
	
	/*! @} */

}

#endif // VGTL_COMPLEX_BUILDER_HPP
