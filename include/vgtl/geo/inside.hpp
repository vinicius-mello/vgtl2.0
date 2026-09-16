#ifndef VGTL_INSIDE_HPP
#define VGTL_INSIDE_HPP

/*! \file
 * \brief Reports if a point is inside of a cell
 */

#include <utility>
#include <vgtl/alg/point.hpp>
#include <vgtl/top/euclidean.hpp>
#include <vgtl/top/maubach.hpp>

namespace vgtl {
	
	using std::pair;
	
	/*! \addtogroup geo 
	 * @{
	 */

	//! Tests if a point is inside a cell
	template <class T, class Point>
	bool
	inside(const T& t,
				 Cell(T) s,
				 const Point& p) {
		const dim_t dim=Dim(T);
		array<Point, dim+1> ps;
		euclidean_points(t,s,ps);
		return inside(p,ps);
	}

	template <class T, class Scalar>
	Cell(T)
	side(const T& t, Cell(T) ce0, Cell(T) ce1,
			 const point<Dim(T),Scalar>& p) {
		pair<int,int> ind=common_facet_ind(t,ce0,ce1);
		Facet(T) f=face_op(t,ce0,ind.first);
		array<point<Dim(T),Scalar>,Dim(T)> psf;
		euclidean_points(t,f,psf);
		array<point<Dim(T),Scalar>,Dim(T)+1> ps;
    for(int i=0; i<ind.first; ++i) ps[i]=psf[i];
    ps[ind.first]=p;
    for(int i=ind.first; i<Dim(T); ++i) ps[i+1]=psf[i];
		Scalar o=orientation(t,ce0)*((ind.first%2)?-1:1);
		if((o*bracket(ps))>0) return ce0;
		else return ce1;
	}

	
	
	//! Returns the first cell in <code>[begin_initial,end_initial)</code> containing \c p
	template <class T, class Iterator, class Point>
	Cell(T) 
	locate(const T& t,
				 Iterator begin_initial, Iterator end_initial,
				 const Point& p) {
		for(Iterator i=begin_initial;i!=end_initial;++i) {
			Cell(T) cv=*i;
			bool flag=inside(t,cv,p);
			while(flag) {
				if(!is_current(t,cv)) {
					Cell(T) ch[2];
					ch[0]=child(t,cv,0);
					ch[1]=child(t,cv,1);
					cv=side(t,ch[0],ch[1],p);
				} else {
					return cv;
				}
			}
		}
	  //throw "Locate Error";
		Cell(T) r;
		empty_set(r);
		return r;
	}
	
	template <class T, class Point>
	Cell(T)
	locate(const T& t,
				 const Point& p) {
		return locate(t,initial_cells(t).begin(),
			initial_cells(t).end(),p);
	}
	
	//! Performs a sequence of splitings on <code>cv=locate(t,begin_initial,end_initial,p)</code> until <code>test(t,cv,p)==true</code>
	template <class T, class Iterator, class Point, class Apply, class Test>
	Cell(T) 
	maubach_subdivide_until(T& t,
				 Iterator begin_initial, Iterator end_initial,
				 const Point& p, Test& test, Apply& app) {
		Cell(T) cv=locate(t,begin_initial,end_initial,p);
		while(!test(t,cv,p)) {
			maubach_subdivide(t,cv,app);
			Cell(T) ch[2];
			ch[0]=child(t,cv,0);
			ch[1]=child(t,cv,1);
			cv=side(t,ch[0],ch[1],p);
		}
	}
	
	template <class T, class Point, class Apply, class Test>
	Cell(T) 
	maubach_subdivide_until(T& t,
				 const Point& p, Test& test, Apply& app) {
		Cell(T) cv=locate(t,p);
		while(!test(t,cv,p)) {
			maubach_subdivide(t,cv,app);
			Cell(T) ch[2];
			ch[0]=child(t,cv,0);
			ch[1]=child(t,cv,1);
			cv=side(t,ch[0],ch[1],p);
		}
	}
	
	/*! @} */

}

#endif // VGTL_INSIDE_HPP
