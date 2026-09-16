#ifndef VGTL_ALG_PLUCKER
#define VGTL_ALG_PLUCKER

/*! \file
 * \brief Models a pluckerernion Space
 */

#include <vgtl/alg/point.hpp>

namespace vgtl {


  /*! \addtogroup alg
   * @{
   */

  //! Plucker Coordinates
  /*! Provides the usual Plucker coordinates operations.
   *
   */
  template <typename Scalar>
  class plucker {
    vec<3,Scalar> d_;
    vec<3,Scalar> c_;
    public:
    plucker() {}
    plucker(const plucker<Scalar>& p) : d_(p.d_), c_(p.c_) {}
    plucker(const point<3,Scalar>& p, const vec<3,Scalar>& dir) {
      d_=dir;
      point<3,Scalar> zero;
      c_=cross(dir,p-zero);
    }
    plucker(const point<3,Scalar>& p, const point<3,Scalar>& q) {
      d_=q-p;
      point<3,Scalar> zero;
      c_=cross(q-zero,p-zero);
    }
    vec<3,Scalar> d() const {
      return d_;
    }
    vec<3,Scalar> c() const {
      return c_;
    }
  };

  template <typename Scalar>
  double operator*(const plucker<Scalar>& p, const plucker<Scalar>& q) {
    return dot(p.d(),q.c())+dot(q.d(),p.c());
  }

  /*! @} */

}

#endif //VGTL_ALG_PLUCKER
