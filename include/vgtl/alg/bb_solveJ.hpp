#ifndef VGTL_BB_SOLVEJ_HPP
#define VGTL_BB_SOLVEJ_HPP

/*! \file
 * \brief BB solve
 */

#include <vgtl/alg/bb.hpp>
#include <vgtl/alg/linalg.hpp>

namespace vgtl {


	/*! \addtogroup alg 
	 * @{
	 */

  template <dim_t n, class Scalar>
	point<n,Scalar> solveJ(const bbform<n,Scalar>& bb,
												 const point<n,Scalar>& p,
												 array<vec<n,Scalar>,n>& J,
											 		int max_iter=20, double eps=1e-6) {
		array<vec<n-1,Scalar>,n-1> A;
    array<dim_t,n-1> per;
		point<n,Scalar> q=p;
		vec<n,Scalar> f;
		f=jacobian(bb,q,J)-p;
		vec<n,Scalar> dw;

		int iter=0;
		bool conv;
		do {
			iter++;
      for(dim_t i=0;i<(n-1);++i)
      for(dim_t j=0;j<(n-1);++j) A[i][j]=J[i][j]-J[i][n-1];
			LU_decomp(A,per);
			LU_solve(A,per,reinterpret_cast<vec<n-1,Scalar>& >(f),
				reinterpret_cast<vec<n-1,Scalar>& >(dw));
			q[n-1]=1;
			dw[n-1]=0;
      for(dim_t i=0;i<(n-1);++i) {
				Scalar e=dw[i];
				dw[i]=-e;
				q[i]=q[i]-e;
				q[n-1]-=q[i];
				dw[n-1]-=dw[i];
			}
			f=jacobian(bb,q,J)-p;
			conv=true;
      for(dim_t i=0;i<n;++i)
				if(abs(f[i])>eps) {
					conv=false;
					break;
				}
		} while((!conv) && iter<max_iter);
		return q;
	}

	/*! @} */

}

#endif // VGTL_BB_SOLVEJ_HPP
