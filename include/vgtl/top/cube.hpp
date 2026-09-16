#ifndef VGTL_CUBE_HPP
#define VGTL_CUBE_HPP

/*! \file
 * \brief Generic function to add a cube from its vertices
 */

#include <vgtl/comb/permutations.hpp>
#include <vgtl/top/complex_builder.hpp>
#include <vgtl/top/add_simplex.hpp>

namespace vgtl {
	
	/*! \addtogroup top 
	 * @{
	 */

	//! Adds a topological hypercube
	/*! \c vs must point to a sequence of \f$2^n\f$ vertices,
	 * with \f$n=\f$<code>sc_traits\<T\>::dim</code>.
	 */
	template <class T, class ArrayIter>
	void
	add_cube(complex_builder<T>& cb, ArrayIter vs) {
		const dim_t dim=Dim(T);
		permutations_traverser<dim> p;
		do {
			unsigned int i=0;
			array<Vertex(T),dim+1> a;
			array<unsigned int,dim> per=*p;
			a[0]=vs[0];
			for(dim_t j=0; j<dim; ++j) {
				i|=1<<(per[j]);
				a[j+1]=vs[i];
			}
			orientation_set(cb.t,add(cb,a),signature(per.begin(),per.end()));
		} while(++p);
	}

	/*! @} */

}

#endif // VGTL_CUBE_HPP
