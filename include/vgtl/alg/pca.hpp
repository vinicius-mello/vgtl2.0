#ifndef VGTL_PCA_HPP
#define VGTL_PCA_HPP

/*! \file
 * \brief Covariance matrix and Principal Component Analysis
 */

#include <vgtl/alg/point.hpp>
#include <vgtl/alg/linalg.hpp>

namespace vgtl {


	/*! \addtogroup alg 
	 * @{
	 */

	//! Builds the covariance matrix from a set of points
  template <typename Iterator, dim_t n, typename Scalar>
  void covariance(Iterator begin, Iterator end, array<vec<n,Scalar>,n>& cov) {
    for(dim_t i=0; i<n; ++i) {
      for(dim_t j=0; j<=i; ++j) {
				Scalar txy=0;
				Scalar tx=0;	
				Scalar ty=0;
				int l=0;
				for(Iterator k=begin; k!=end; ++k, ++l) {
	  			point<n, Scalar> p=*k;
	  			txy+=p[i]*p[j];
	  			tx+=p[i];
	  			ty+=p[j];
				}
				cov[i][j]=(txy-tx*ty/l)/(l-1);
			}
		}
    for(dim_t j=1; j<n; ++j) {
      for(dim_t i=0; i<j; ++i) {
				cov[i][j]=cov[j][i];
			}
		}
	}

	//! Performs PCA and returns the ellipsoid axis and radii
  template <typename Iterator, dim_t n, typename Scalar>
	void pca(Iterator begin, Iterator end, point<n,Scalar>& bc,
	    array<vec<n,Scalar>,n>& axes,
	    vec<n,Scalar>& a) {
    barycenter(begin,end,bc);
    array<vec<n,Scalar>,n> cov;
    covariance(begin,end,cov);
    eigenvalues(cov,axes,a);
  }

	/*! @} */
}

#endif // VGTL_PCA_HPP
