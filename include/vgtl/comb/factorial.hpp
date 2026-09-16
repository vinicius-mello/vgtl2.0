#ifndef VGTL_FACTORIAL_HPP
#define VGTL_FACTORIAL_HPP

/*! \file
 * \brief Factorial function
 */

#include <vgtl/utl/array.hpp>

namespace vgtl {

	/*! \addtogroup comb 
	 * @{
	 */
	
	//! Factorial function
	inline
	int factorial(int n) {
	static int lut[14]={
		1,
		1,
		2,
		6,
		24,
		120,
		720,
		5040,
		40320,
		362880,
		3628800,
		39916800,
		479001600,
		1932053504};
		return lut[n];
	}
	
	/*! @} */

}

#endif // VGTL_FACTORIAL_HPP
