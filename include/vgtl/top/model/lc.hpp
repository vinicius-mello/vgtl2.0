#ifndef VGTL_LC_HPP
#define VGTL_LC_HPP

/*! \file
 * \brief Models a APM using lists
 */

#include <list>
#include <utility>
#include <vgtl/top/incidence.hpp>
#include <vgtl/top/add.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/top_cell.hpp>
#include <vgtl/utl/array.hpp>
#include <vgtl/top/model/bs.hpp>

namespace vgtl {

	using std::list;
	using std::pair;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <dim_t d,
						template <dim_t> class Data>
	struct lc;
					
	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	struct simplex_descriptor<lc<d,Data>,k> {
		typename list<bs<lc<d,Data>,k,Data> >::iterator desc;
		simplex_descriptor& operator=(const simplex_descriptor& _s) {
			desc=_s.desc;
			return (*this);
		}
		bool operator<(const simplex_descriptor& _s) const {
			return (&*desc)<(&*_s.desc);
		}
		bool operator==(const simplex_descriptor& _s) const {
			return desc==_s.desc;
		}
		bool operator!=(const simplex_descriptor& _s) const {
			return desc!=_s.desc;
		}
		simplex_descriptor() {
			void * * ptr=reinterpret_cast<void * *>(&desc);
			*ptr=0;
		}
		simplex_descriptor(typename list<bs<lc<d,Data>,k,Data> >::iterator _desc)
						: desc(_desc) {}
	};

	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	inline
	bool
	empty(simplex_descriptor<lc<d,Data>,k> s) {
		void * * ptr=reinterpret_cast<void * *>(&s.desc);
		return (*ptr==0);
	}

	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	inline
	void
	empty_set(simplex_descriptor<lc<d,Data>,k>& s) {
		void * * ptr=reinterpret_cast<void * *>(&s.desc);
		*ptr=0;
	}

	template <dim_t d, 
	          template <dim_t> class Data>
	struct sc_traits<lc<d,Data> > {
		static const dim_t dim=d;
		static const dim_t split_dim=0;
		typedef apm_tag sc_category; 
	};

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	struct lc_rep : lc_rep <k-1,d,Data> {
		list<bs<lc<d,Data>,k,Data> > a;
	};

	template <dim_t d,
	          template <dim_t> class Data>
	struct lc_rep<0,d,Data> {
		list<bs<lc<d,Data>,0,Data> > a;
	};

	//! Models a Mutable APM
	template <dim_t d,
						template <dim_t> class Data=extra_data>
	struct lc : lc_rep<d,d,Data> , Data<d+1> {
		static const dim_t dim=d;
		typedef apm_tag sc_category;
	};
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	list<bs<lc<d,Data>,k,Data> >& 
	simplex_container(lc<d,Data>& t) {
		return t.lc_rep<k,d,Data>::a; 
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	const list<bs<lc<d,Data>,k,Data> >& 
	simplex_container(const lc<d,Data>& t) {
		return t.lc_rep<k,d,Data>::a; 
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	simplex_descriptor<lc<d,Data>,k>
	add(lc<d,Data>& t, simplex_descriptor<lc<d,Data>,k>& sd) {
		bs<lc<d,Data>,k,Data> s;
		simplex_container<k,d,Data>(t).push_front(s);
		sd=simplex_descriptor<lc<d,Data>,k>(simplex_container<k,d,Data>(t).begin());
		return sd;
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	bs<lc<d,Data>,k,Data> *
	attr(lc<d,Data>& t, simplex_descriptor<lc<d,Data>,k> s) {
		return &*s.desc;
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	const bs<lc<d,Data>,k,Data> *
	attr(const lc<d,Data>& t, simplex_descriptor<lc<d,Data>,k> s) {
		return &*s.desc;
	}

	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	struct simplex_iterator<lc<d,Data>, k> {
		simplex_descriptor<lc<d,Data>, k> p;
		simplex_iterator(typename list<bs<lc<d,Data>,k,Data> >::iterator i)
						: p(i) {}
		simplex_iterator() {}
		simplex_iterator& operator++() {
			++p.desc;
			return *this;
		}
		simplex_descriptor<lc<d,Data>, k> operator*() {
			return p;
		}
		bool operator!=(const simplex_iterator& i) {
			return p!=i.p;
		}
		bool operator==(const simplex_iterator& i) {
			return p==i.p;
		}
	};

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	simplex_descriptor<lc<d,Data>,k-1>
	face_op(const lc<d,Data>& t,
					simplex_descriptor<lc<d,Data>,k> s,
	        dim_t i) {
		return attr(t,s)->fop[i]; 
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	simplices(const lc<d,Data>& t,
						simplex_iterator<lc<d,Data>,k>& begin,
						simplex_iterator<lc<d,Data>,k>& end) {
		simplex_descriptor<lc<d,Data>,k> b(const_cast<list<bs<lc<d,Data>,k,Data> >&>(simplex_container<k,d,Data>(t)).begin()); 
		simplex_descriptor<lc<d,Data>,k> e(const_cast<list<bs<lc<d,Data>,k,Data> >&>(simplex_container<k,d,Data>(t)).end()); 
		begin.p=b;
		end.p=e;
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	simplices(lc<d,Data>& t,
						simplex_iterator<lc<d,Data>,k>& begin,
						simplex_iterator<lc<d,Data>,k>& end) {
		simplex_descriptor<lc<d,Data>,k> b(simplex_container<k,d,Data>(t).begin()); 
		simplex_descriptor<lc<d,Data>,k> e(simplex_container<k,d,Data>(t).end()); 
		begin.p=b;
		end.p=e;
	}

	template <dim_t d,
	          template <dim_t> class Data>
	inline 
	std::pair<simplex_descriptor<lc<d,Data>,d>,
						simplex_descriptor<lc<d,Data>,d> >
	cells(const lc<d,Data>& t,
				simplex_descriptor<lc<d,Data>,d-1> s) {
		return std::make_pair(attr(t,s)->cov[0], attr(t,s)->cov[1]);
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline 
	typename disable_if<equal_c<k+1,d>,
		simplex_descriptor<lc<d,Data>,k+1> >::type
	up_simplex(const lc<d,Data>& t,
						 simplex_descriptor<lc<d,Data>,k> s) {
		return attr(t,s)->up;
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline 
	typename enable_if<equal_c<k+1,d>,
		simplex_descriptor<lc<d,Data>,k+1> >::type
	up_simplex(const lc<d,Data>& t,
						 simplex_descriptor<lc<d,Data>,k> s) {
		return attr(t,s)->cov[0];
	}
	
  //add
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void face_op_set(lc<d,Data>& t,
						 			 simplex_descriptor<lc<d,Data>,k> s,
									 dim_t i,
						 			 simplex_descriptor<lc<d,Data>,k-1> r) {
		attr(t,s)->fop[i]=r; 
	}
	
	template <dim_t d,
	          template <dim_t> class Data>
	inline 
	void
	cells_set(lc<d,Data>& t,
						simplex_descriptor<lc<d,Data>,d-1> s,
						simplex_descriptor<lc<d,Data>,d> r0,
						simplex_descriptor<lc<d,Data>,d> r1) {
		attr(t,s)->cov[0]=r0;
		attr(t,s)->cov[1]=r1;
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	up_simplex_set(lc<d,Data>& t,
						     simplex_descriptor<lc<d,Data>,k> s,
								 simplex_descriptor<lc<d,Data>,k+1> r) {
		attr(t,s)->up=r;
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	int orientation(const lc<d,Data>& t,
						   	  simplex_descriptor<lc<d,Data>,d> s) {
		return attr(t,s)->ori ? -1 : 1;
	}

	template <dim_t d, 
	          template <dim_t> class Data>
	inline
	void orientation_set(lc<d,Data>& t,
						  	  		 simplex_descriptor<lc<d,Data>,d> s, int o) {
		attr(t,s)->ori= (o<0) ? true : false;
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	mark_type mark(const lc<d,Data>& t,
					 simplex_descriptor<lc<d,Data>,d> s) {
		return attr(t,s)->mk;
	}

	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void mark_set(const lc<d,Data>& t,
						 	  simplex_descriptor<lc<d,Data>,d> s, mark_type i) {
		attr(const_cast<lc<d,Data>&>(t),s)->mk=i;
	}


	//del
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	void
	del(lc<d,Data>& t, simplex_descriptor<lc<d,Data>,k> s) {
		simplex_container<k,d,Data>(t).erase(s.desc);
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	int
	level(const lc<d,Data>& t,
	  	  simplex_descriptor<lc<d,Data>,d> s) {
		return (int)attr(t,s)->lvl;
	}

	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void
	level_set(lc<d,Data>& t,
	  	 			simplex_descriptor<lc<d,Data>,d> s, int i) {
		attr(t,s)->lvl=(short)i;
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	bool is_current(const lc<d,Data>& t,
						  	  simplex_descriptor<lc<d,Data>,d> s) {
		return attr(t,s)->cur;
	}

	template <dim_t k, dim_t d, 
	          template <dim_t> class Data>
	inline
	typename disable_if<cell_c<lc<d,Data>,k>,bool>::type
	is_current(const lc<d,Data>& t,
		       	 simplex_descriptor<lc<d,Data>,k> s) {
		return true;
	}
	
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void is_current_set(lc<d,Data>& t,
						  	  		simplex_descriptor<lc<d,Data>,d> s, bool b) {
		attr(t,s)->cur=b;
	}


	/*! @} */

}

#endif // VGTL_LC_HPP
