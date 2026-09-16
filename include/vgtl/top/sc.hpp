#ifndef VGTL_SC_HPP
#define VGTL_SC_HPP

/*! \file
 * \brief Simplicial Complex definitions
 */	

#include <cstddef>
#include <vgtl/utl/constraint.hpp>

//@{
#define Dim(T) vgtl::sc_traits<T>::dim
#define Simplex(T,N) vgtl::simplex_descriptor<T,N>
#define Vertex(T) Simplex(T,0)
#define Edge(T) Simplex(T,1)
#define Facet(T) Simplex(T,Dim(T)-1)
#define Cell(T) Simplex(T,Dim(T))
#define Simplex_it(T,N) vgtl::simplex_iterator<T,N>
#define Vertex_it(T) Simplex_it(T,0)
#define Edge_it(T) Simplex_it(T,1)
#define Facet_it(T) Simplex_it(T,Dim(T)-1)
#define Cell_it(T) Simplex_it(T,Dim(T))
#define SplitDim(T) vgtl::sc_traits<T>::split_dim
#define SplitSimplex(T) Simplex(T,SplitDim(T))
//@}
	
namespace vgtl {
				
	/*! \addtogroup top 
	 * @{
	 */
	
	//! Simplex descriptor type
	template <class T, dim_t d>
	struct simplex_descriptor {
		simplex_descriptor();
		simplex_descriptor& operator=(const simplex_descriptor& s);
		bool operator<(const simplex_descriptor& s) const;
		bool operator==(const simplex_descriptor& s) const;
	};

	//! Simplex iterator type
	template <class T, dim_t d>
	struct simplex_iterator {
		simplex_iterator();
		simplex_iterator& operator++();
		simplex_descriptor<T,d> operator*();
		bool operator==(const simplex_iterator& i);
	};

	//! Simplicial Complex trait class
	template <class T>
	struct sc_traits {
		static const dim_t dim=T::dim;
		static const dim_t split_dim=T::split_dim;
		typedef typename T::sc_category sc_category; 
	};

	//! Abstract Simplicial Complex Concept tag
	struct asc_tag {};
	// Valid expressions:
	// void
	// simplices(const T& t,
	// 					 simplex_iterator<T,k>& begin, simplex_iterator<T,k>& end);
	// simplex_descriptor<T,i-1> 
	// face_op(const T& t, simplex_descriptor<T,i> s, int j)

	//! Abstract Pseudo-Manifold Concept tag
	struct apm_tag : virtual public asc_tag {};
	// Valid expressions:
	// simplex_descriptor<T,i+1> 
	// up_simplex(const T& t, simplex_descriptor<T,i> s)
	// 
	// pair<sc_traits<T>::cell_descriptor,
	//      sc_traits<T>::cell_descriptor> 
	// cells(const T& t, sc_traits<T>::facet_descriptor ce)

	//! Multitriangulations tag
	struct mt_tag : virtual public apm_tag {};
	
	/*! \name Constraints
	 * Dimensional constraints
	 */
	//@{
	//! Vertex constraint
	template <class T, dim_t k>
	struct vertex_c {
		static const bool value=(k==0);
	};

	//! Edge constraint
	template <class T, dim_t k>
	struct edge_c {
		static const bool value=(k==1);
	};

	//! Edge not facet constraint
	template <class T, dim_t k>
	struct edge_not_facet_c {
		static const bool value=(k==1)&&(Dim(T)!=2);
	};

	//! Cell constraint
	template <class T, dim_t k>
	struct cell_c {
		static const bool value=(k==Dim(T));
	};

	//! Not vertex or cell constraint
	template <class T, dim_t k>
	struct not_vertex_or_cell_c {
		static const bool value=(k!=Dim(T))&&(k!=0);
	};

	//! Facet constraint
	template <class T, dim_t k>
	struct facet_c {
		static const bool value=(k==(Dim(T)-1));
	};

	//! Facet or cell constraint
	template <class T, dim_t k>
	struct facet_or_cell_c {
		static const bool value=(k==(Dim(T)-1))||(k==Dim(T));
	};

	//! Otherwise constraint
	template <class T, dim_t k>
	struct otherwise_c {
		static const bool value=((k>1)&&(k<(Dim(T)-1)));
	};

	//! SplitSimplex constraint
	template <class T, dim_t k>
	struct splitsimplex_facet_c {
		static const bool value=((SplitDim(T)==k)&&((Dim(T)-1)==k));
	};

	//! SplitSimplex constraint
	template <class T, dim_t k>
	struct splitsimplex_not_facet_c {
		static const bool value=((SplitDim(T)==k)&&((Dim(T)-1)!=k));
	};
	//@}
	enum mark_type {white_mark, black_mark, gray_mark};
	/*! @} */

}

#endif // VGTL_SC_HPP
