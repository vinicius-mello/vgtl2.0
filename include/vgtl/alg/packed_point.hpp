#ifndef VGTL_PACKED_POINT_HPP
#define VGTL_PACKED_POINT_HPP

/*! \file
 * \brief A packed point class usefull in dimensions 2 and 3
 */

#include <vgtl/alg/point.hpp>

namespace vgtl {

	/*! \addtogroup alg 
	 * @{
	 */

	//! A packed point class
	/*! Usefull if you only need dyadic points in \f$[0,1]^2\f$ or \f$[0,1]^3\f$.
	 *
	 * \b Example:
	 * \include packed_point.cpp 
	 * \b Output:
	 * \include packed_point.out
	 */
  template <dim_t n, class Scalar>
	class packed_point {};
	
	template <class Scalar>
	struct packed_point<2,Scalar> {
		packed_point() : x(0),y(0) {};
		packed_point(const packed_point& pp) {
			x=pp.x;
			y=pp.y;
		}
		packed_point(const point<2,Scalar>& p) {
			x=static_cast<unsigned int>(p[0]*32768);
			y=static_cast<unsigned int>(p[1]*32768);
		}
		operator point<2,Scalar>() const {
			point<2,Scalar> temp;
			temp[0]=x/32768.0f;
			temp[1]=y/32768.0f;
			return temp;
		}
		unsigned int x : 16;
		unsigned int y : 16;
	};


	template <class Scalar>
	inline
	packed_point<2,Scalar> mediator(const packed_point<2,Scalar>& p0,
																	const packed_point<2,Scalar>& p1) {
		packed_point<2,Scalar> temp;
		temp.x=p0.x+(p1.x-p0.x)/2;
		temp.y=p0.y+(p1.y-p0.y)/2;
		return temp;
	}

	template <class Scalar>
	struct packed_point<3,Scalar> {
		packed_point() : x(0),y(0),z(0),fill(0) {};
		packed_point(const packed_point& pp) {
			x=pp.x;
			y=pp.y;
			z=pp.z;
			fill=0;
		}
		packed_point(const point<3,Scalar>& p) {
			x=static_cast<unsigned int>(p[0]*32768);
			y=static_cast<unsigned int>(p[1]*32768);
			z=static_cast<unsigned int>(p[2]*32768);
			fill=0;
		}
		operator point<3,Scalar>() const {
			point<3,Scalar> temp;
			temp[0]=x/32768.0f;
			temp[1]=y/32768.0f;
			temp[2]=z/32768.0f;
			return temp;
		}
		unsigned int fill : 16;
		unsigned int x : 16;
		unsigned int y : 16;
		unsigned int z : 16;
	};

	template <class Scalar>
	inline
	packed_point<3,Scalar> mediator(const packed_point<3,Scalar>& p0,
																	const packed_point<3,Scalar>& p1) {
		packed_point<3,Scalar> temp;
		temp.x=p0.x+(p1.x-p0.x)/2;
		temp.y=p0.y+(p1.y-p0.y)/2;
		temp.z=p0.z+(p1.z-p0.z)/2;
		return temp;
	}
	
	template <class Scalar>
	struct packed_point<4,Scalar> {
		packed_point() : x(0),y(0),z(0),w(0) {};
		packed_point(const packed_point& pp) {
			x=pp.x;
			y=pp.y;
			z=pp.z;
			w=pp.w;
		}
		packed_point(const point<4,Scalar>& p) {
			x=static_cast<unsigned int>(p[0]*256);
			y=static_cast<unsigned int>(p[1]*256);
			z=static_cast<unsigned int>(p[2]*256);
			w=static_cast<unsigned int>(p[3]*256);
		}
		operator point<4,Scalar>() const {
			point<4,Scalar> temp;
			temp[0]=x/256.0f;
			temp[1]=y/256.0f;
			temp[2]=z/256.0f;
			temp[3]=w/256.0f;
			return temp;
		}
		unsigned int x : 8;
		unsigned int y : 8;
		unsigned int z : 8;
		unsigned int w : 8;
	};

	template <class Scalar>
	inline
	packed_point<4,Scalar> mediator(const packed_point<4,Scalar>& p0,
																	const packed_point<4,Scalar>& p1) {
		packed_point<4,Scalar> temp;
		temp.x=p0.x+(p1.x-p0.x)/2;
		temp.y=p0.y+(p1.y-p0.y)/2;
		temp.z=p0.z+(p1.z-p0.z)/2;
		temp.w=p0.w+(p1.w-p0.w)/2;
		return temp;
	}

	/*! @} */

}

#endif //VGTL_PACKED_POINT_HPP
