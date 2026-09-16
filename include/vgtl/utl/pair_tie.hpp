#ifndef VGTL_PAIR_TIE_HPP
#define VGTL_PAIR_TIE_HPP

/*! \file
 * \brief Pair tie
 */

#include <utility>

namespace vgtl {

	/*! \addtogroup utl 
	 * @{
	 */

	template <class T, class U>
	struct rpair {
		T& first;
		U& second;
		rpair(T& a, U& b) : first(a), second(b) {}
		rpair& operator=(const std::pair<T,U>& p) {first=p.first;second=p.second; return *this;}
	};
	
	//! Great for a two variables assignment
	template <class T, class U>
	rpair<T,U> pair_tie(T& a,U& b) {
		return rpair<T,U>(a,b);
	}

	/*! @} */

}

#endif //VGTL_PAIR_TIE_HPP
