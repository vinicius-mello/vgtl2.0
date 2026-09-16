#ifndef VGTL_FILL_HPP
#define VGTL_FILL_HPP

/*! \file
 * \brief Recursevely fills a complex_builder/complex_buffer
 * with subsimplices of a simplex
 */

#include <vgtl/top/complex_builder.hpp>
#include <vgtl/top/complex_buffer.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T, dim_t k>
	void
	fill(complex_builder<T>& cb, Simplex(T,k) s) {
		array<Vertex(T),k+1> vs;
		vertices(cb.t,s,vs);
		if(!exists(cb,vs)) {
			put(cb,vs,s);
			for(dim_t i=0; i<=k; ++i) {
				Simplex(T,k-1) fs;
				fs=face_op(cb.t,s,i);
				fill(cb,fs);
			}
		}
	}
	
	template <class T>
	void
	fill(complex_builder<T>& cb, Edge(T) s) {
		array<Vertex(T),2> vs;
		vertices(cb.t,s,vs);
		if(!exists(cb,vs)) put(cb,vs,s);
	}
	
	template <class T, dim_t k>
	void
	fill(complex_buffer<T>& cb, Simplex(T,k) s) {
		if(!exists(cb,s)) {
			put(cb,s);
			for(dim_t i=0; i<=k; ++i) {
				Simplex(T,k-1) fs;
				fs=face_op(cb.t,s,i);
				fill(cb,fs);
			}
		}
	}
	
	template <class T>
	void
	fill(complex_buffer<T>& cb, Edge(T) s) {
		if(!exists(cb,s)) put(cb,s);
	}
	
	template <class T>
	void
	fill_vc(complex_buffer<T>& cb, Cell(T) ce) {
		put(cb,ce);
		array<Vertex(T),Dim(T)+1> vs;
		vertices(cb.t,ce,vs);
		for(dim_t i=0; i<=Dim(T); ++i) {
			put(cb,vs[i]);
		}
	}
	
	template <class T, dim_t k, dim_t l>
	typename enable_if<greater_c<k,l>,void>::type
	fill(complex_buffer<T>& cb,
		   Simplex(T,k) s,
		   Simplex(T,l) e) {
		if(!exists(cb,s)) {
			put(cb,s);
			for(dim_t i=0; i<=k; ++i) {
				Simplex(T,k-1) fs;
				fs=face_op(cb.t,s,i);
				if(in(cb.t,fs,e)) fill(cb,fs,e);
			}
		}
	}
	
	template <class T, dim_t k, dim_t l>
	typename enable_if<equal_c<k,l>,void>::type
	fill(complex_buffer<T>& cb,
		   Simplex(T,k) s,
		   Simplex(T,l) e) {
		if(!exists(cb,s)) put(cb,s);
	}
	
	template <class T, dim_t k, dim_t l>
	typename enable_if<lesser_c<k,l>,void>::type
	fill(complex_buffer<T>& cb,
		   Simplex(T,k) s,
		   Simplex(T,l) e) {
	}

	/*! @} */
	
}

#endif // VGTL_FILL_HPP
