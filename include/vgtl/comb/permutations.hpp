#ifndef VGTL_PERMUTATIONS_HPP
#define VGTL_PERMUTATIONS_HPP

/*! \file
 * \brief Permutations traverser
 */

#include <algorithm>
#include <vgtl/utl/array.hpp>

namespace vgtl {

  using std::next_permutation;
			
	template <class ArrayIt>
	int signature(ArrayIt begin, ArrayIt end) {
		int inversions=0;
		for(ArrayIt i=begin;i!=end;++i) {
			for(ArrayIt j=i+1;j!=end;++j) {
				if(*i>*j) inversions++;
			}
		}
		return (inversions%2) ? -1 : 1;
	}
	
	/*! \addtogroup comb 
	 * @{
	 */
	
	//! Permutations traverser
	/*! Traverses the permutations on \c n elements.
	 *
	 * \b Example:
	 * \include permutations.cpp 
	 * \b Output:
	 * \include permutations.out
	 */
  template <dim_t n>
	class permutations_traverser {
		array<unsigned int,n> v;
		public:
		typedef array<unsigned int,n> result_type;
		permutations_traverser()
		{
			for(unsigned int i=0; i<n; ++i) v[i]=i;
		}
		bool operator++() {
			return next_permutation(&v[0],&v[0]+n);
		}
		array<unsigned int,n>& operator*() {
			return v;
		}
		const array<unsigned int,n>& operator*() const {
			return v;
		}
	};


	/*! @} */

}

#endif // VGTL_PERMUTATIONS_HPP
