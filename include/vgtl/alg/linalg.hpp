#ifndef VGTL_LINALG_HPP
#define VGTL_LINALG_HPP

/*! \file
 * \brief Wrappers to GSL linear algebra functions
 */

#include <vector>
#include <vgtl/alg/vec.hpp>
#include <vgtl/alg/point.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

namespace vgtl {

	using std::vector;

	/*! \addtogroup alg 
	 * @{
	 */

	//! Converts vgtl::vec to gsl_vector
  template <dim_t n>
	void convert(vec<n,double>& a, gsl_vector& gv) {
		gv.size=n;
		gv.data=&a[0];
		gv.stride=1;
		gv.owner=0;
	}

	//! Converts stl::vector to gsl_vector
	inline
	void convert(vector<double>& a, gsl_vector& gv) {
		gv.size=a.size();
		gv.data=&a[0];
		gv.stride=1;
		gv.owner=0;
	}

	//! Converts vgtl::array<vec<n>,m> to gsl_matrix
  template <dim_t n, dim_t m>
	void convert(array<vec<n,double>,m>& a, gsl_matrix& gm) {
		gm.size1=m;
		gm.size2=n;
		gm.tda=n;
		gm.data=&a[0][0];
		gm.owner=0;
	}

	//! Converts std::vector to gsl_matrix
	inline
  void convert(dim_t m, dim_t n, vector<double>& a, gsl_matrix& gm) {
		gm.size1=m;
		gm.size2=n;
		gm.tda=n;
		gm.data=&a[0];
		gm.owner=0;
	}

  //! Converts vgtl::array<dim_t,n> to gsl_permutation
  /*! gsl_permutation::data is a size_t*, not an unsigned int*, so (on a
   *  platform where the two differ in size) it can't just alias `a`'s
   *  storage: `buf` (n size_t's, owned by the caller) backs `gp` instead,
   *  populated from `a`. Call convert_back afterwards to copy any change
   *  GSL made in `buf` (e.g. after gsl_linalg_LU_decomp) back into `a`. */
  template <dim_t n>
  void convert(const array<dim_t,n>& a, gsl_permutation& gp, size_t* buf) {
		gp.size=n;
		gp.data=buf;
		for(dim_t i=0;i<n;++i) buf[i]=a[i];
	}

	template <dim_t n>
	void convert_back(const size_t* buf, array<dim_t,n>& a) {
		for(dim_t i=0;i<n;++i) a[i]=static_cast<dim_t>(buf[i]);
	}

	//! Converts std::vector to gsl_permutation -- see the array<dim_t,n>
	//! overload above for why a separate size_t buffer is needed.
	inline
  void convert(const vector<dim_t>& a, gsl_permutation& gp, size_t* buf) {
		gp.size=a.size();
		gp.data=buf;
		for(size_t i=0;i<a.size();++i) buf[i]=a[i];
	}

	inline
	void convert_back(const size_t* buf, vector<dim_t>& a) {
		for(size_t i=0;i<a.size();++i) a[i]=static_cast<dim_t>(buf[i]);
	}
	
	//! Performs LU decomposition of a in place, with permutation p
  template <dim_t n>
  int LU_decomp(array<vec<n,double>,n>& a, array<dim_t,n>& p) {
		gsl_permutation gp;
		gsl_matrix gm;
		convert(a,gm);
		size_t pbuf[n];
		convert(p,gp,pbuf);
    int s;
    gsl_linalg_LU_decomp(&gm,&gp,&s);
    convert_back(pbuf,p);
    return s;
	}

	//! Performs LU decomposition of a in place, with permutation p
	inline
  int LU_decomp(dim_t n, vector<double>& a, vector<dim_t>& p) {
		gsl_permutation gp;
		gsl_matrix gm;
		convert(n,n,a,gm);
		vector<size_t> pbuf(p.size());
		convert(p,gp,&pbuf[0]);
    int s;
    gsl_linalg_LU_decomp(&gm,&gp,&s);
    convert_back(&pbuf[0],p);
    return s;
	}

	//! Solves Ax=b of previously LU factored matrix A
  template <dim_t n>
  void LU_solve(array<vec<n,double>,n>& a, array<dim_t,n>& p,
								vec<n,double>& b, vec<n,double>& x) {
		gsl_permutation gp;
		gsl_matrix gm;
		gsl_vector gvb;
		gsl_vector gvx;
		convert(a,gm);
		size_t pbuf[n];
		convert(p,gp,pbuf);
		convert(b,gvb);
		convert(x,gvx);
		gsl_linalg_LU_solve(&gm,&gp,&gvb,&gvx);
	}

	//! Solves Ax=b of previously LU factored matrix A
	inline
  void LU_solve(dim_t n, vector<double>& a, vector<dim_t>& p,
								vector<double>& b, vector<double>& x) {
		gsl_permutation gp;
		gsl_matrix gm;
		gsl_vector gvb;
		gsl_vector gvx;
		convert(n,n,a,gm);
		vector<size_t> pbuf(p.size());
		convert(p,gp,&pbuf[0]);
		convert(b,gvb);
		convert(x,gvx);
		gsl_linalg_LU_solve(&gm,&gp,&gvb,&gvx);
	}

//! Solves Ax=b of previously LU factored matrix A in place
  template <dim_t n>
  void LU_solve_in_place(array<vec<n,double>,n>& a, array<dim_t,n>& p,
								vec<n,double>& bx) {
		gsl_permutation gp;
		gsl_matrix gm;
		gsl_vector gvbx;
		convert(a,gm);
		size_t pbuf[n];
		convert(p,gp,pbuf);
		convert(bx,gvbx);
		gsl_linalg_LU_svx(&gm,&gp,&gvbx);
	}

	//! Computes the determinant of previously LU factored matrix A
  template <dim_t n>
	double LU_det(array<vec<n,double>,n>& a, int sig) { 
		gsl_matrix gm;
		convert(a,gm);
		return gsl_linalg_LU_det(&gm,sig);
	}

	//! Inverts previously LU factored matrix A
  template <dim_t n>
  void LU_invert(array<vec<n,double>,n>& a, array<dim_t,n>& p,
			array<vec<n,double>,n>& inv) { 
		gsl_permutation gp;
		gsl_matrix gm;
		gsl_matrix ginv;
		convert(a,gm);
		size_t pbuf[n];
		convert(p,gp,pbuf);
		convert(inv,ginv);
		gsl_linalg_LU_invert(&gm,&gp,&ginv);
	}

	//! Computes a.v
  template <dim_t n, dim_t m>
	vec<m,double> apply(array<vec<n,double>,m>& a, vec<n,double>& v) {
		gsl_matrix gm;
		gsl_vector gv;
		gsl_vector gvr;
		vec<m,double> vr;
		convert(a,gm);
		convert(v,gv);
		convert(vr,gvr);
		gsl_blas_dgemv(CblasNoTrans,1,&gm,&gv,0,&gvr);
		return vr;
	}


	//! Computes a.v (float)
  template <dim_t n, dim_t m>
	vec<m,float> apply(array<vec<n,float>,m>& a, vec<n,float>& v) {
		vec<m,float> vr;
    for(dim_t i=0;i<m;++i) vr[i]=dot(a[i],v);
		return vr;
	}
	

	//! Computes eigenvectors and eingenvalues of a
  template <dim_t n>
	void eigenvalues(array<vec<n,double>,n>& a, 
			 array<vec<n,double>,n>& eigenvectors, 
			 vec<n,double>& eigenvalues, bool asc=true) { 
		gsl_matrix ga;
		gsl_matrix gev;
		gsl_vector gea;
		convert(a,ga);
		convert(eigenvectors,gev);
		convert(eigenvalues,gea);
 		gsl_eigen_symmv_workspace * ws=gsl_eigen_symmv_alloc(n);
 		gsl_eigen_symmv(&ga, &gea, &gev, ws);
 		if(asc) gsl_eigen_symmv_sort(&gea, &gev, GSL_EIGEN_SORT_VAL_ASC);
		else gsl_eigen_symmv_sort(&gea, &gev, GSL_EIGEN_SORT_VAL_DESC);
		gsl_matrix_transpose(&gev);
 		gsl_eigen_symmv_free(ws);
	}

	/*! @} */

}

#endif // VGTL_LINALG_HPP
