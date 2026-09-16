#ifndef VGTL_COMPLEX_BUFFER_HPP
#define VGTL_COMPLEX_BUFFER_HPP

/*! \file
 * \brief A buffer to store simplices of different dimensions
 */

#include <set>
#include <vgtl/top/sc.hpp>
#include <vgtl/top/do_nothing.hpp>

namespace vgtl {
	
	using std::set;
	
	/*! \addtogroup top 
	 * @{
	 */

	template <class T, dim_t k>
	struct complex_buffer_rep : complex_buffer_rep<T,k-1> {
		set<Simplex(T,k)> s;
		complex_buffer_rep<T,k>(T& _t) : complex_buffer_rep<T,k-1>(_t) {}
	};

	template <class T>
	struct complex_buffer_rep<T,0> {
		set<Vertex(T)> s;
		T& t;
		complex_buffer_rep<T,0>(T& _t) : t(_t) {}
	};


	//! A buffer to store simplices of different dimensions
	/*! A complex buffer its just a set of <code>sc_traits\<T\>::dim+1</code>
	 * sets of simplices, one for each dimension. It works like a heterogeneus
	 * set of <code>simplex_descriptor</code>.*/
	template <class T>
	struct complex_buffer : complex_buffer_rep<T,Dim(T)> {
		complex_buffer(T& _t) : complex_buffer_rep<T,Dim(T)>(_t) {}
		template <dim_t k>
		typename set<Simplex(T,k)>::const_iterator 
		begin() const {
			return (*this).complex_buffer_rep<T,k>::s.begin(); 
		}
		template <dim_t k>
		typename set<Simplex(T,k)>::const_iterator 
		end() const {
			return (*this).complex_buffer_rep<T,k>::s.end(); 
		}
		template <dim_t k>
		typename set<Simplex(T,k)>::iterator 
		begin() {
			return (*this).complex_buffer_rep<T,k>::s.begin(); 
		}
		template <dim_t k>
		typename set<Simplex(T,k)>::iterator 
		end() {
			return (*this).complex_buffer_rep<T,k>::s.end(); 
		}

	};

	//! Returns the set of \c k dimensional simplices
	template <dim_t k, class T>
	set<Simplex(T,k)>&
	buffer_set(complex_buffer<T> &cb) {
	  return cb.complex_buffer_rep<T,k>::s;
	}

	//! Returns the set of \c k dimensional simplices
	template <dim_t k, class T>
	const set<Simplex(T,k)>&
	buffer_set(const complex_buffer<T> &cb) {
	  return cb.complex_buffer_rep<T,k>::s;
	}

	//! Tests if a simplex belongs to a complex buffer
	template <class T, dim_t k>
	bool
	exists(const complex_buffer<T> &cb, 
		     Simplex(T,k) key)
	{
	  if(buffer_set<k>(cb).count(key)==0) return false;
		return true;
	}

	//! Puts a simplex into a complex buffer
	template <class T, dim_t k>
	void
	put(complex_buffer<T> &cb, 
		  Simplex(T,k) s)
	{
	  buffer_set<k>(cb).insert(s);
	}

	template <class T, class Apply, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	apply(complex_buffer<T> &cb, Apply& app) {
	  for(typename set<Simplex(T,k)>::iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			app.apply(cb.t,*i);
		}
	}
	
	template <class T, class Apply, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	apply(complex_buffer<T> &cb, Apply& app) {
	  for(typename set<Simplex(T,k)>::iterator i=cb.template begin<k>();
				i!=cb.template end<k>(); ++i) {
			app.apply(cb.t,*i);
		}
		apply<T,Apply,k+1>(cb,app);
	}
	
	//! Applies a function to all simplices in a complex buffer
	template <class T, class Apply>
	void
	apply(complex_buffer<T> &cb, Apply& app) {
		apply<T,Apply,0>(cb,app);
	}
	
	template <class T>
	inline
	void
	apply(complex_buffer<T> &cb, do_nothing& app) {
	}
	
	/*! @} */

}

#endif // VGTL_COMPLEX_BUFFER_HPP
