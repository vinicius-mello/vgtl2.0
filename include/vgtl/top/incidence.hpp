#ifndef VGTL_INCIDENCE_HPP
#define VGTL_INCIDENCE_HPP

/*! \file
 * \brief Incidence relations
 */

#include <vgtl/comb/combinations.hpp>
#include <vgtl/top/sc.hpp>
#include <vgtl/utl/array.hpp>

namespace vgtl {


	/*! \addtogroup top 
	 * @{
	 */
	
	//! Composition of \c face_op
	/*! \invariant <code>face_op(t,s,_(a,b,c))==
	 * face_op(t,face_op(t,face_op(t,s,c),b,),a)</code>
	 */
	template <class T, dim_t k, dim_t l>
	Simplex(T,k-l) 
	face_op(const T& t, 
				  		Simplex(T,k) s,
				  		const array<dim_t,l>& a) {
		const array<dim_t,l-1>& b=
			reinterpret_cast<const array<dim_t,l-1>&>(a[1]);
		return face_op(t,face_op(t,s,b),a[0]);
	}

	template <class T, dim_t k>
	Simplex(T,k-1) 
	face_op(const T& t, 
				  		Simplex(T,k) s,
				  		const array<dim_t,1>& a) {
		return face_op(t,s,a[0]);
	}

	//! Traverses the faces of a simplex
	/*! Traverses the \e l dimensional faces of a \e k dimensional simplex.
	 *  It has the traverser semantics.
	 */
	template <class T, dim_t k, dim_t l>
	class face_traverser {
		const T& t;
		combinations_lex_rev_traverser<k-l> clt;
		Simplex(T,k) gs;
		public:
		face_traverser(const T& _t, Simplex(T,k) _gs)
                   : t(_t), clt(k+1), gs(_gs) {}
		bool operator++() {
			return ++clt;
		}
		Simplex(T,l) operator*() {
			return face_op(t,gs,*clt);
		}
		
		Simplex(T,l) operator*() const {
			return face_op(t,gs,*clt);
		}
		
	};

	//! Traverses the faces of a simplex
	/*! Traverses the \e l dimensional faces of a \e k dimensional simplex.
	 *  It has the traverser semantics.
	 */
	template <class T, dim_t k, dim_t l>
	class face_traverser_rev {
		const T& t;
		combinations_lex_traverser<k-l> clt;
		Simplex(T,k) gs;
		public:
		face_traverser_rev(const T& _t, Simplex(T,k) _gs)
						       : t(_t), gs(_gs), clt(k+1) {}
		bool operator++() {
			return ++clt;
		}
		Simplex(T,l) operator*() {
			return face_op(t,gs,*clt);
		}
		
		Simplex(T,l) operator*() const {
			return face_op(t,gs,*clt);
		}
		
	};

	//! Reports the vertices of a simplex
	/*! \c a is updated by reference.
	 */
	template <class T, dim_t k>
	void 
	vertices(const T& t, 
					 Simplex(T,k) s,
					 array<Vertex(T),k+1>& a) {
		face_traverser<T,k,0> ft(t,s);
		dim_t i=0;
		do {
			a[i]=*ft;
			++i;
		} while(++ft);
	}
	
	template <class T>
	void 
	vertices(const T& t, 
					 Vertex(T) s,
					 array<Vertex(T),1>& a) {
		a[0]=s;
	}
	
	//! Tests if \c ss is a face of \c sg
	template <dim_t k, dim_t l, class T>
	typename enable_if<greater_c<k,l>,bool>::type 
	in(const T& t, 
		 Simplex(T,k) sg,
		 Simplex(T,l) ss) {
		face_traverser<T,k,l> ft(t,sg);
		do {
			if(*ft==ss) return true;
		} while(++ft);
		return false;
	}
	
	template <dim_t k, dim_t l, class T>
	typename enable_if<lesser_c<k,l>,bool>::type 
	in(const T& t, 
		 Simplex(T,k) sg,
		 Simplex(T,l) ss) {
		return false;
	}
	
	template <dim_t k, dim_t l, class T>
	typename enable_if<equal_c<k,l>,bool>::type 
	in(const T& t, 
		 Simplex(T,k) sg,
		 Simplex(T,l) ss) {
		return sg==ss;
	}
	
	//! Returns the opposite face of \c ss in \c sg
	template <dim_t k, dim_t l, class T>
	Simplex(T,k-l-1) 
	opposite(const T& t, 
		 			 Simplex(T,k) sg,
					 Simplex(T,l) ss) {
		face_traverser<T,k,l> ft(t,sg);
		face_traverser_rev<T,k,k-l-1> ft_rev(t,sg);
		do {
			if(*ft==ss) return *ft_rev;
			++ft_rev;
		} while(++ft);
		return Simplex(T,k-l-1)();
		
	}

	//! Returns the indices of \c ss in \c sg
	/*! \invariant <code>face_op(t,s,indices(t,s,f))==f</code>
	 */
	template <dim_t k, dim_t l, class T>
	array<dim_t,k-l>
	indices(const T& t,
					Simplex(T,k) sg,
					Simplex(T,l) ss) {
		combinations_lex_traverser<k-l> clt(k+1);
		do {
			Simplex(T,l) st=face_op(t,sg,*clt);
			if(st==ss) {
				return *clt;
			}
		} while(++clt);
	}
	
	template <dim_t k, dim_t l, class T>
	typename enable_if<equal_c<k+1,l>,Simplex(T,l-1)>::type
	face_ind(const T& t,
					Simplex(T,k) s,
					const array<dim_t,l>& a) {
		return s;
	}

	//! Returns the subsimplex of \c s with indices given by \c a
	template <dim_t k, dim_t l, class T>
	typename enable_if<greater_c<k+1,l>,Simplex(T,l-1)>::type
	face_ind(const T& t,
					Simplex(T,k) s,
					const array<dim_t,l>& a) {
		array<dim_t,k+1-l> b;
    dim_t i=0;
    dim_t ii=0;
    for(dim_t j=0;j<=k;++j) {
			if((i<l)&&(j==a[i])) ++i;
			else {
				b[ii]=j;
				++ii;
			}
		}
		return face_op(t,s,b);
	}

	/*! @} */

}

#endif // VGTL_INCIDENCE_HPP
