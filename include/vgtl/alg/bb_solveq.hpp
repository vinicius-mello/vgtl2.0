#ifndef VGTL_BB_SOLVEQ_HPP
#define VGTL_BB_SOLVEQ_HPP

/*! \file
 * \brief BB solveq
 */

#include <vgtl/alg/bb.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

namespace vgtl {

	/*! \addtogroup alg 
	 * @{
	 */

  template <dim_t n, class Scalar>
	struct bbq_params {
		const bbform<n+1,Scalar>& bb;
		array<vec<n,Scalar>,n>& Jf;
		const point<n,Scalar>& initial;
		bbq_params(const bbform<n+1,Scalar>& _bb, const point<n,Scalar>& _initial,
							array<vec<n,Scalar>,n>& _Jf) 
			: bb(_bb), initial(_initial), Jf(_Jf) {};
	};

  template <dim_t n, class Scalar>
	int bbq_f(const gsl_vector * x, void * params,
					 gsl_vector * f) {
		bbq_params<n,Scalar> * par=(bbq_params<n,Scalar> *)params;
		point<n,Scalar> p;
		p[n-1]=1;
    for(dim_t i=0; i<(n-1); ++i) {
			p[i]=(Scalar)gsl_vector_get(x,i);
			p[n-1]-=p[i];
		}
		point<n,Scalar> q=evalq(par->bb,p);
    for(dim_t i=0; i<(n-1); ++i) {
			gsl_vector_set(f,i,q[i]-par->initial[i]);
		}
		return GSL_SUCCESS;
	}
	
  template <dim_t n, class Scalar>
	int bbq_df(const gsl_vector * x, void * params,
				    gsl_matrix * J) {
		bbq_params<n,Scalar> * par=(bbq_params<n,Scalar> *)params;
		point<n,Scalar> p;
		p[n-1]=1;
    for(dim_t i=0; i<(n-1); ++i) {
			p[i]=(Scalar)gsl_vector_get(x,i);
			p[n-1]-=p[i];
		}
//		array<vec<n,Scalar>,n> Jf;
		point<n,Scalar> q=jacobianq(par->bb,p,par->Jf);
    for(dim_t i=0; i<(n-1); ++i) {
      for(dim_t j=0; j<(n-1); ++j) {
				gsl_matrix_set(J,i,j,par->Jf[i][j]-par->Jf[i][n-1]);
			}
		}
		return GSL_SUCCESS;
	}
	
  template <dim_t n, class Scalar>
	int bbq_fdf(const gsl_vector * x, void * params,
						 gsl_vector * f, gsl_matrix * J) {
		bbq_params<n,Scalar> * par=(bbq_params<n,Scalar> *)params;
		point<n,Scalar> p;
		p[n-1]=1;
    for(dim_t i=0; i<(n-1); ++i) {
			p[i]=(Scalar)gsl_vector_get(x,i);
			p[n-1]-=p[i];
		}
//		array<vec<n,Scalar>,n> Jf;
		point<n,Scalar> q=jacobianq(par->bb,p,par->Jf);
    for(dim_t i=0; i<(n-1); ++i) {
			gsl_vector_set(f,i,q[i]-par->initial[i]);
		}
    for(dim_t i=0; i<(n-1); ++i) {
      for(dim_t j=0; j<(n-1); ++j) {
				gsl_matrix_set(J,i,j,par->Jf[i][j]-par->Jf[i][n-1]);
			}
		}
		return GSL_SUCCESS;
	}
	
  template <dim_t n, class Scalar>
	point<n,Scalar> solveqJ(const bbform<n+1,Scalar>& bb,
											 	 const point<n,Scalar>& p,
												 array<vec<n,Scalar>,n>& Jf,
												 int max_iter=20, double eps=1e-6) {
		bbq_params<n,Scalar> par(bb,p,Jf);
		gsl_multiroot_fdfsolver * s;
		gsl_multiroot_function_fdf f={bbq_f<n,Scalar>,
																	bbq_df<n,Scalar>,
																	bbq_fdf<n,Scalar>,n-1,&par};
		gsl_vector * x=gsl_vector_alloc(n-1);
    for(dim_t i=0; i<(n-1); ++i) {
			gsl_vector_set(x,i,p[i]);
		}
		s=gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj,n-1);
		//s=gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridj,n-1);
		//s=gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_newton,n-1);
		//s=gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_gnewton,n-1);
		gsl_multiroot_fdfsolver_set(s,&f,x);
		int iter=0;
		int status;
		do {
			iter++;
			status=gsl_multiroot_fdfsolver_iterate(s);
			if(status) break;
			status=gsl_multiroot_test_residual(s->f,eps);
		} while(status==GSL_CONTINUE && iter<max_iter);
		point<n,Scalar> q;
		q[n-1]=1;
    for(dim_t i=0; i<(n-1); ++i) {
			q[i]=(Scalar)gsl_vector_get(s->x,i);
			q[n-1]-=q[i];
		}
		gsl_multiroot_fdfsolver_free(s);
		gsl_vector_free(x);
		return q;
	}

  template <dim_t n, class Scalar>
	point<n,Scalar> solveq(const bbform<n+1,Scalar>& bb,
											 const point<n,Scalar>& p,
											 int max_iter=20, double eps=1e-6) {
		array<vec<n,Scalar>,n> Jf;
		return solveqJ(bb,p,Jf,max_iter,eps);
	}

	/*! @} */

}

#endif // VGTL_BB_SOLVEQ_HPP
