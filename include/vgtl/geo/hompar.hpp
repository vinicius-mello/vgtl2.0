#ifndef VGTL_GEO_HOMPAR_HPP
#define VGTL_GEO_HOMPAR_HPP

#include <vgtl/top/incidence.hpp>
#include <vgtl/utl/pair_tie.hpp>
#include <vgtl/alg/bb.hpp>

/*! \file
 * \brief 
 */


namespace vgtl {

	
	/*! \addtogroup geo 
	 * @{
	 */

	template <class T, class Attr, dim_t k>
	typename disable_if<greater_c<k,0>,void>::type
	hompar(const T& t, 
				Simplex(T,k) s, const array<int,k+1>& cs, Attr ** at) {
		*at=hompar(t,s,0);
	}

	/*! @} */
	template <class T, class Attr, dim_t k>
	typename enable_if<greater_c<k,0>,void>::type
	hompar(const T& t, 
				Simplex(T,k) s, const vec<k+1,int>& cs, Attr ** at) {
		for(int i=0;i<=k;++i) {
			if(cs[i]==0) {
				vec<k,int> csp;
				project(i,cs,csp);
				hompar(t,face_op(t,s,i),csp,at);
				return;
			}
		}
		vec<k+1,int> cs1;
		for(int i=0;i<=k;++i) cs1[i]=cs[i]-1;
		*at=hompar(t,s,hom2ind(cs1));
	}

	template <class T>
	double
	rdiff_c1_error(const T& t, Facet(T) f) {
		int deg=degree(t);
		Cell(T) ce[2];
		pair_tie(ce[0],ce[1])=cells(t,f);
		int i0,i1;
		pair_tie(i0,i1)=common_facet_ind(t,ce[0],ce[1]);
    vec<Dim(T),int> v;
		v[0]=deg-1;
		double err=0;
    do {
			double err1=0;
      double * hp;
      vec<Dim(T)+1,int> uv;
      unproject(i1,v,uv);
      uv[i1]=1;
      double ss=0;
      array<Vertex(T),Dim(T)+1> vs;
      vertices(t,ce[1],vs);
      for(int k=0;k<=Dim(T);++k) ss+=uv[k]*scalar_value(t,vs[k]);
      hompar(t,ce[1],uv,&hp);
			err1-=ss*(*hp);
			point<Dim(T),double> p=euclidean_point(t,vs[i1]);
      vertices(t,ce[0],vs);
			array<point<Dim(T),double>,Dim(T)+1> ps;
      for(int k=0;k<=Dim(T);++k) ps[k]=euclidean_point(t,vs[k]);
      point<Dim(T)+1,double> a=barycentric_coordinates(p,ps);
      uv[i1]=0;
      swap(uv[i1],uv[i0]);
      for(int j=0;j<=Dim(T);++j) {
        double * hpp;
        uv[j]+=1;
        double s=0;
        for(int k=0;k<=Dim(T);++k) s+=uv[k]*scalar_value(t,vs[k]);
        hompar(t,ce[0],uv,&hpp);
        err1+=a[j]*(*hpp)*s;
        uv[j]-=1;
      }
			err+=err1*err1;
    } while(next_triangular(deg-1,v.begin(),v.end()));
		return err;
	}

	/*! @} */

}

#endif // VGTL_GEO_HOMPAR_HPP
