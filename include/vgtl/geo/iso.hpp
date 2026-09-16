#ifndef VGTL_GEO_ISO_HPP
#define VGTL_GEO_ISO_HPP

#include <vgtl/top/iso/iso.hpp>
#include <vgtl/alg/point.hpp>

namespace vgtl {

	namespace iso {

		template <class T, class Scalar>
		void basis(const T& t,
							 Cell(T) ce,
							 const array<Vertex(T),Dim(T)+1>& vs,
							 array<point<Dim(T)+1,Scalar>,Dim(T)>& bps) {
      array<array<dim_t,2>,Dim(T)> ai;
			basis(t,ce,vs,ai);
      for(dim_t j=0; j<Dim(T); ++j) {
				Scalar a=scalar_value(t,vs[ai[j][0]]);
				Scalar b=scalar_value(t,vs[ai[j][1]]);
				bps[j][ai[j][0]]=(-b)/(a-b);
				bps[j][ai[j][1]]=(a)/(a-b);
			}
		}

		template <class T, dim_t k, class Scalar>
		void basis(const T& t,
							 const array<Vertex(T),k+1>& vs,
							 array<point<k+1,Scalar>,k>& bps) {
      array<array<dim_t,2>,k> ai;
			basis(t,vs,ai);
      for(dim_t j=0; j<k; ++j) {
				Scalar a=scalar_value(t,vs[ai[j][0]]);
				Scalar b=scalar_value(t,vs[ai[j][1]]);
				bps[j][ai[j][0]]=(-b)/(a-b);
				bps[j][ai[j][1]]=(a)/(a-b);
			}
		}
		
	}

}

#endif // VGTL_GEO_ISO_HPP
