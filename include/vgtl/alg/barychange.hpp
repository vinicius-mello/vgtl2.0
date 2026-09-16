#ifndef VGTL_BARYCHANGE_HPP
#define VGTL_BARYCHANGE_HPP

/*! \file
 * \brief Change from cartesian to barycentric coordinates with respect a given simplex
 */

#include <vgtl/alg/point.hpp>
#include <vgtl/alg/linalg.hpp>

namespace vgtl {

	/*! \addtogroup alg 
	 * @{
	 */

	//! Change from cartesian to barycentric coordinates with respect a given simplex
  template <dim_t n, class Scalar>
	class barychange {
		array<vec<n+1,Scalar>,n+1> a;
    array<dim_t,n+1> p;
		public:
		//! Constructs the change of coordinates from ps
		barychange(const array<point<n,Scalar>,n+1>& ps) {
      for(dim_t i=0;i<=n;++i) a[0][i]=1;
      for(dim_t j=1;j<=n;++j)
      for(dim_t i=0;i<=n;++i) a[j][i]=ps[i][j-1];
			LU_decomp(a,p);
		}
		//! Converts affine point to barycentric point wrt ps
		point<n+1,Scalar> operator()(const point<n,Scalar>& q) {
			vec<n+1,Scalar> b;
			vec<n+1,Scalar> x;
			b[0]=1;
      for(dim_t i=1;i<=n;++i) b[i]=q[i-1];
			LU_solve(a,p,b,x);
			point<n+1,Scalar> r;
      for(dim_t i=0;i<=n;++i) r[i]=x[i];
			return r;
		}
		//! Converts affine vector to barycentric vector wrt ps
		vec<n+1,Scalar> operator()(const vec<n,Scalar>& v) {
			vec<n+1,Scalar> b;
			vec<n+1,Scalar> x;
			b[0]=0;
      for(dim_t i=1;i<=n;++i) b[i]=v[i-1];
			LU_solve(a,p,b,x);
			vec<n+1,Scalar> r;
      for(dim_t i=0;i<=n;++i) r[i]=x[i];
			return r;
		}
	};


	/*! @} */

}

#endif // VGTL_BARYCHANGE_HPP
