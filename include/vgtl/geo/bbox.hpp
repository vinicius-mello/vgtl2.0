#ifndef VGTL_GEO_BBOX_HPP
#define VGTL_GEO_BBOX_HPP

#include <vector>
#include <algorithm>
#include <vgtl/alg/point.hpp>

/*! \file
 * \brief Bounding box class and functions 
 * \b Example:
 * \include bbox.cpp 
 */


namespace vgtl {

  using namespace std;
	
	/*! \addtogroup geo 
	 * @{
	 */

  template <dim_t n, class Scalar=double>
	struct bbox {
		point<n,Scalar> b[2];

		bbox() {
		}
		bbox(const point<n,Scalar>& p0, const point<n,Scalar>& p1) {
			b[0]=p0;
			b[1]=p1;
		}
		point<n,Scalar>& operator[](int i) {
			return b[i];
		}
		const point<n,Scalar>& operator[](int i) const {
			return b[i];
		}
    dim_t largest_side() const {
      dim_t ls=0;
			Scalar side=b[1][ls]-b[0][ls];
      for(dim_t i=1;i<n;++i) {
				Scalar t=b[1][i]-b[0][i];
				if(t>side) {
					side=t;
					ls=i;
				}
			}
			return ls;
		}
		void set(const bbox<n,Scalar>& b0, const bbox<n,Scalar>& b1) {
      for(dim_t j=0;j<n;++j) {
				b[0][j]=min(b0[0][j],b1[0][j]);
				b[1][j]=max(b0[1][j],b1[1][j]);
			}
		}
		template <class BBoxIt>
		void boxes(BBoxIt begin, BBoxIt end) {
			if(begin==end) return;
			BBoxIt i=begin;
			b[0]=i->b[0];
			b[1]=i->b[1];
			for(++i; i!=end; ++i) {
        for(dim_t j=0;j<n;++j) {
					b[0][j]=min(b[0][j],i->b[0][j]);
					b[1][j]=max(b[1][j],i->b[1][j]);
				}
			}
		}
		template <class PointsIt>
		void points(PointsIt begin, PointsIt end) {
			if(begin==end) return;
			PointsIt i=begin;
			b[0]=*i;
			b[1]=*i;
			for(++i; i!=end; ++i) {
        for(dim_t j=0;j<n;++j) {
					b[0][j]=min(b[0][j],(*i)[j]);
					b[1][j]=max(b[1][j],(*i)[j]);
				}
			}
		}
		Scalar min_distance(const point<n,Scalar>& p) const {
			vec<n,Scalar> d0, d1;
			Scalar dmin;
			d0=p-b[0];
			d1=b[1]-p;
			dmin=0;
      for(dim_t i=0;i<n;++i) {
				d0[i]*=d0[i];
				d1[i]*=d1[i];
				dmin+=p[i]<b[0][i] ? d0[i] : p[i]>b[1][i] ? d1[i] : 0;
			}
			return dmin;
		}
		void minmax_distance(const point<n,Scalar>& p, Scalar& rmin, Scalar& rmax) const {
			vec<n,Scalar> d0, d1, M, m;
			Scalar dmin, dmax;
			d0=p-b[0];
			d1=b[1]-p;
			dmin=0;
      for(dim_t i=0;i<n;++i) {
				d0[i]*=d0[i];
				d1[i]*=d1[i];
				dmin+=p[i]<b[0][i] ? d0[i] : p[i]>b[1][i] ? d1[i] : 0;
				m[i]=min(d0[i],d1[i]);
				M[i]=max(d0[i],d1[i]);
			}
      for(dim_t i=0;i<n;++i) {
				Scalar t=0;
        for(dim_t j=0;j<n;++j) {
					if(j==i) continue;
					t+=M[j];
				}
				if(i==0) dmax=m[i]+t;
				else dmax=min(dmax,m[i]+t);
			}
			rmin=dmin;
			rmax=dmax;
		}
	};

  template <dim_t n, class Scalar=double>
	struct bbox_tree {
		bbox<n,Scalar> bb;
		bbox_tree<n,Scalar> * left;
		bbox_tree<n,Scalar> * right;
		bbox_tree() : left(0), right(0) {}
		~bbox_tree() {
			if(left) delete left;
			if(right) delete right;
		}
	};

  template <dim_t n, class Scalar>
	struct bbox_point_cut {
		Scalar cut;
		int ls;
		bbox_point_cut(Scalar _cut, int _ls) : cut(_cut), ls(_ls) {}
		bool operator()(const point<n,Scalar>& p) const {
			return p[ls]<cut;
		}
	};
	
  template <dim_t n, class Scalar>
	bbox_tree<n,Scalar> *
	build_bbox_tree(typename vector<point<n,Scalar> >::iterator i,
		typename vector<point<n,Scalar> >::iterator end) {
		if(end==i) return 0;
		bbox_tree<n,Scalar> * bbt
						=new bbox_tree<n,Scalar>();
		bbt->bb.points(i,end);
		if((end-i)==1) return bbt;
		//cout<<bbt->bb[0]<<endl;
		//cout<<bbt->bb[1]<<endl;
		int ls=bbt->bb.largest_side();
		//cout<<ls<<endl;
		Scalar cut=(bbt->bb[0][ls]+bbt->bb[1][ls])/2;
		//cout<<cut<<endl;
		typename vector<point<n,Scalar> >::iterator middle;
		middle=partition(i,end,bbox_point_cut<n,Scalar>(cut,ls));
		bbt->left=build_bbox_tree<n,Scalar>(i,middle);
		bbt->right=build_bbox_tree<n,Scalar>(middle,end);
		return bbt;
	}

  template <dim_t n, class Scalar>
	bbox_tree<n,Scalar> *
	build_bbox_tree(vector<point<n,Scalar> >& v) {
		return build_bbox_tree<n,Scalar>(v.begin(),v.end());
	}
	
  template <dim_t n, class Scalar>
	void points_distance(bbox_tree<n,Scalar> * bbt,
										const point<n,Scalar>& p, Scalar& dist) {
		if((bbt->left==0)&&(bbt->right==0)) {
			vec<n,Scalar> r=bbt->bb[0]-p;
			dist=min(dist,dot(r,r));
			return;
		}
		Scalar rmin, rmax;
		bbt->bb.minmax_distance(p,rmin,rmax);
		if(rmin>dist) return;
		dist=min(dist,rmax);
		if(bbt->left) points_distance(bbt->left,p,dist);
		if(bbt->right) points_distance(bbt->right,p,dist);
	}
	
	/*! @} */

}

#endif // VGTL_GEO_BBOX_HPP
