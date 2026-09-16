#ifndef VGTL_CONSTRAINT_HPP
#define VGTL_CONSTRAINT_HPP

/*! \file
 * \brief Constraint definitions
 */	

#include <vgtl/utl/types.hpp>
#include <vgtl/utl/enable_if.hpp>


namespace vgtl {

				
	template <dim_t k, dim_t d>
	struct equal_c {
		static const bool value=(k==d);
	};

	//! Different dimension constraint
	template <dim_t k, dim_t d>
	struct different_c {
		static const bool value=(k!=d);
	};
	
	//! Lesser dimension constraint
	template <dim_t k, dim_t d>
	struct lesser_c {
		static const bool value=(k<d);
	};

	//! Greater dimension constraint
	template <dim_t k, dim_t d>
	struct greater_c {
		static const bool value=(k>d);
	};

}

#endif // VGTL_CONSTRAINT_HPP
