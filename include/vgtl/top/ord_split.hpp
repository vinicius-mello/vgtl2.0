#ifndef VGTL_ORD_SPLIT_HPP
#define VGTL_ORD_SPLIT_HPP

/*! \file
 * \brief Generic functions to split a simplex
 */

#include <list>
#include <set>
#include <map>
#include <vgtl/top/mark.hpp>
#include <vgtl/top/decompose.hpp>
#include <vgtl/top/unlink.hpp>

namespace vgtl {

	using std::map;

	/*! \addtogroup top 
	 * @{
	 */
	
	//! Performs a stellar subdivision on \c s
	template <class T, dim_t k, class StarIterator>
	void
	ord_split(T& t,
				Simplex(T,k) s, 
				StarIterator begin, StarIterator end,
				Vertex(T) v, 
				dim_t ind,
				complex_buffer<T>& cb_in,
				complex_buffer<T>& cb, apm_tag) {
		const dim_t dim=Dim(T);
		complex_buffer<T> cb_notin(t);
		decompose(t,s,begin,end,cb_in,cb_notin);
		unlink(cb_in,cb_notin);
		complex_builder<T> cbl(t);
		buffer_to_build(cb_notin,cbl);
	  for(StarIterator i=begin;i!=end; ++i) {
			Cell(T) cv=*i;
			is_current_set(t,cv,false);
			unvisit(t,cv);
			for(int j=0;j<=dim;++j) {
				Facet(T) ce=face_op(t,cv,j);
				if(exists(cb_notin,ce)) {
					array<Vertex(T),dim+1> vs;
					vertices(t,ce,reinterpret_cast<array<Vertex(T),dim>&>(vs));
					for(dim_t r=dim; r>ind; --r) vs[r]=vs[r-1];
					vs[ind]=v;
					Cell(T) cvn=add(cbl,vs);
					fill(cb,cvn,v);
					orientation_set(t,cvn,
													orientation(t,cv)*(j%2 ? -1 : 1)*(ind%2 ? -1 : 1));
					level_set(t,cvn,level(t,cv)+1);
					is_current_set(t,cvn,true);
				}
			}
		}
	}
	
	template <class T, dim_t k, class StarIterator>
	void
	ord_split(T& t,
				Simplex(T,k) s, 
				StarIterator begin, StarIterator end,
				Vertex(T) v, 
				dim_t ind,
				complex_buffer<T>& cb_in,
				complex_buffer<T>& cb, map<Cell(T),Cell(T)>& parent) {
		const dim_t dim=Dim(T);
		complex_buffer<T> cb_notin(t);
		decompose(t,s,begin,end,cb_in,cb_notin);
		unlink(cb_in,cb_notin);
		complex_builder<T> cbl(t);
		buffer_to_build(cb_notin,cbl);
	  for(StarIterator i=begin;i!=end; ++i) {
			Cell(T) cv=*i;
			is_current_set(t,cv,false);
			unvisit(t,cv);
			for(int j=0;j<=dim;++j) {
				Facet(T) ce=face_op(t,cv,j);
				if(exists(cb_notin,ce)) {
					array<Vertex(T),dim+1> vs;
					vertices(t,ce,reinterpret_cast<array<Vertex(T),dim>&>(vs));
					for(dim_t r=dim; r>ind; --r) vs[r]=vs[r-1];
					vs[ind]=v;
					Cell(T) cvn=add(cbl,vs);
					parent[cvn]=cv;
					fill(cb,cvn,v);
					orientation_set(t,cvn,
													orientation(t,cv)*(j%2 ? -1 : 1)*(ind%2 ? -1 : 1));
					level_set(t,cvn,level(t,cv)+1);
					is_current_set(t,cvn,true);
				}
			}
		}
	}
	
	//! Performs a stellar subdivision on \c s
	template <class T, dim_t k, class StarIterator>
	void
	ord_split(T& t,
				Simplex(T,k) s, 
				StarIterator begin, StarIterator end,
				Vertex(T) v, 
				dim_t ind,
				complex_buffer<T>& cb_in,
				complex_buffer<T>& cb, mt_tag) {
		const dim_t dim=Dim(T);
		complex_buffer<T> cb_notin(t);
		decompose(t,s,begin,end,cb_in,cb_notin);
		unlink(cb_in,cb_notin);
		complex_builder<T> cbl(t);
		buffer_to_build(cb_notin,cbl);
	  for(StarIterator i=begin;i!=end; ++i) {
			Cell(T) cv=*i;
			is_current_set(t,cv,false);
			unvisit(t,cv);
			int l=0;
			for(int j=dim;j>=0;--j) {
				Facet(T) ce=face_op(t,cv,j);
				if(exists(cb_notin,ce)) {
					array<Vertex(T),dim+1> vs;
					vertices(t,ce,reinterpret_cast<array<Vertex(T),dim>&>(vs));
					for(dim_t r=dim; r>ind; --r) vs[r]=vs[r-1];
					vs[ind]=v;
					Cell(T) cvn=add(cbl,vs);
					parent_set(t,cvn,cv);
					child_set(t,cv,l,cvn);
					++l;
					orientation_set(t,cvn,
													orientation(t,cv)*(j%2 ? -1 : 1)*(ind%2 ? -1 : 1));
					level_set(t,cvn,level(t,cv)+1);
					is_current_set(t,cvn,true);
					fill(cb,cvn,v);
				}
			}
		}
	}

	template <class T, dim_t k, class StarIterator>
	void
	ord_split(T& t,
				Simplex(T,k) s, 
				StarIterator begin, StarIterator end,
				Vertex(T) v, 
				dim_t ind,
				complex_buffer<T>& cb_in,
				complex_buffer<T>& cb) {
		ord_split(t,s,begin,end,v,ind,cb_in,cb,typename sc_traits<T>::sc_category());
	}
	/*! @} */

}

#endif // VGTL_ORD_SPLIT_HPP
