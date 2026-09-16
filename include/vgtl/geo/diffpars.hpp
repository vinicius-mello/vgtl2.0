#ifndef VGTL_GEO_DIFFPARS_HPP
#define VGTL_GEO_DIFFPARS_HPP

#include <vgtl/top/incidence.hpp>
#include <vgtl/comb/triangular.hpp>
#include <vgtl/alg/bb.hpp>

/*! \file
 * \brief 
 */


namespace vgtl {

	
	/*! \addtogroup geo 
	 * @{
	 */

	template <class T, class Scalar, dim_t k, dim_t l>
	void diffpars(const T& t, 
								Simplex(T,k) s, bbform<k+1,Scalar>& bb, Simplex(T,l) ss) {
    const dim_t n=k+1;
    const dim_t m=l+1;
    //triangular_traverser<n,int> tr(m);
    triangular_traverser<n> tr(m);
		bb.clear();
		bb.degree_set(m);
		point<m,Scalar> q;
		Scalar f=1;
		f/=m;
		do {
      vec<n,dim_t> i;
			array<dim_t,m> a;
			i=*tr;
			point<n,Scalar> r;
			bool zerone=true;
			int h=0;
      for(dim_t j=0;j<n;++j) {
				if((i[j]!=0)&&(i[j]!=1)) {
					zerone=false;
					break;
				} else {
					if(i[j]==1) {
						a[h]=j;
						++h;
					}
				}
			}
			if(zerone) {
				q=diffpars(t,face_ind(t,s,a));
        for(dim_t j=0;j<m;++j) r[a[j]]=q[j];
			} else {
        for(dim_t j=0;j<n;++j) r[j]=f*i[j];
			}
			bb.push_back(r);
		} while(++tr);
	}

	template <class T, class Scalar, dim_t k, dim_t l>
	void ex_diffpars(const T& t, 
								Simplex(T,k) s, bbform<k+1,Scalar>& bb, Simplex(T,l) ss) {
    const dim_t n=k+1;
    const dim_t m=l+1;
    //triangular_traverser<n,int> tr(m);
    triangular_traverser<n> tr(m);
		bb.clear();
		bb.degree_set(m);
		point<m,Scalar> q;
		Scalar f=1;
		f/=m;
		do {
      vec<n,dim_t> i;
			array<dim_t,m> a;
			i=*tr;
			point<n,Scalar> r;
			bool zerone=true;
			int h=0;
      for(dim_t j=0;j<n;++j) {
				if((i[j]!=0)&&(i[j]!=1)) {
					zerone=false;
					break;
				} else {
					if(i[j]==1) {
						a[h]=j;
						++h;
					}
				}
			}
			if(zerone) {
				q=ex_diffpars(t,face_ind(t,s,a));
        for(dim_t j=0;j<m;++j) r[a[j]]=q[j];
			} else {
        for(dim_t j=0;j<n;++j) r[j]=f*i[j];
			}
			bb.push_back(r);
		} while(++tr);
	}

	template <class T, class Scalar, dim_t k, dim_t l>
	void diffpars(const T& t, 
								Simplex(T,k) s, int p, bbform<k+1,Scalar>& bb, Simplex(T,l) ss) {
    const dim_t n=k+1;
    const dim_t m=l+1;
    //triangular_traverser<n,int> tr(m);
    triangular_traverser<n> tr(m);
		bb.clear();
		bb.degree_set(m);
		const Scalar * q;
		Scalar f=1;
		f/=m;
		do {
      vec<n,dim_t> i;
			array<dim_t,m> a;
			i=*tr;
			point<n,Scalar> r;
			bool zerone=true;
			int h=0;
      for(dim_t j=0;j<n;++j) {
				if((i[j]!=0)&&(i[j]!=1)) {
					zerone=false;
					break;
				} else {
					if(i[j]==1) {
						a[h]=j;
						++h;
					}
				}
			}
			if(zerone) {
				q=diffpars(t,face_ind(t,s,a));
				q=&q[(k+1)*p];
				r[a[l]]=1;
        for(dim_t j=0;j<l;++j) {
					r[a[j]]=q[j];
					r[a[l]]-=q[j];
				}
			} else {
        for(dim_t j=0;j<n;++j) r[j]=f*i[j];
			}
			bb.push_back(r);
		} while(++tr);
	}


	/*! @} */

}

#endif // VGTL_GEO_DIFFPARS_HPP
