#ifndef VGTL_BARYGRAD_HPP
#define VGTL_BARYGRAD_HPP

/*! \file
	\brief Computes the gradient of the linear function represented by vector \c f in point \c p
 */

#include <vgtl/alg/point.hpp>
#include <vgtl/alg/linalg.hpp>

namespace vgtl {

	/*! \addtogroup alg 
	 * @{
	 */

	//! Computes the gradient of the linear function represented by vector \c f in point \c p
  template <dim_t n>
  vec<n,double> barygrad(const array<point<n,double>,n+1>& p,
												 vec<n+1,double>& f) {
		array<vec<n+1,double>,n+1> Pt;
		vec<n,double> fg;
		vec<n+1,double> t;
    array<dim_t,n+1> per;
    for(dim_t i=0;i<=n; ++i) {
			Pt[i][0]=1;
		}
    for(dim_t i=0;i<n; ++i)
      for(dim_t j=0;j<=n; ++j) Pt[j][i+1]=p[j][i];
		LU_decomp(Pt,per);
		LU_solve(Pt,per,f,t);
		t[0]=0;
    for(dim_t i=0;i<n;++i) fg[i]=t[i+1];
		return fg;
	}

	/*! @} */

}

#endif // VGTL_BARYGRAD_HPP
