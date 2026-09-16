#ifndef VGTL_MODEL_LMT_HPP
#define VGTL_MODEL_LMT_HPP

/*! \file
 * \brief Models a Retangular Triangulation
 */

#include <list>
#include <utility>
#include <vgtl/utl/array.hpp>
#include <vgtl/top/model/ms.hpp>
#include <vgtl/top/incidence.hpp>
#include <vgtl/top/add.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/top_cell.hpp>
#include <vgtl/top/mt.hpp>

namespace vgtl {

	using std::list;
	using std::pair;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <dim_t d, dim_t sd,
						template <dim_t> class Data>
	struct lmt;
	
						
	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	struct simplex_descriptor<lmt<d,sd,Data>,k> {
		typename list<ms<lmt<d,sd,Data>,k,Data> >::iterator desc;
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
		simplex_descriptor(typename list<ms<lmt<d,sd,Data>,k,Data> >::iterator _desc)
						: desc(_desc) {}
	};

	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	inline
	bool
	empty(simplex_descriptor<lmt<d,sd,Data>,k> s) {
		void * * ptr=reinterpret_cast<void * *>(&s.desc);
		return (*ptr==0);
	}

	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	inline
	void
	empty_set(simplex_descriptor<lmt<d,sd,Data>,k>& s) {
		void * * ptr=reinterpret_cast<void * *>(&s.desc);
		*ptr=0;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	struct sc_traits<lmt<d,sd,Data> > {
		static const dim_t dim=d;
		static const dim_t split_dim=sd;
		typedef mt_tag sc_category; 
	};

	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	struct lmt_rep : lmt_rep <k-1,sd,d,Data> {
		list<ms<lmt<d,sd,Data>,k,Data> > a;
	};

	template <dim_t sd, dim_t d,
	          template <dim_t> class Data>
	struct lmt_rep<0,sd,d,Data> {
		list<ms<lmt<d,sd,Data>,0,Data> > a;
	};

	//! Models a Mutable Retangular Triangulation
	template <dim_t d, dim_t sd=1,
						template <dim_t> class Data=extra_data>
	struct lmt : lmt_rep<d,sd,d,Data>, Data<d+1> {
		static const dim_t dim=d;
		static const dim_t split_dim=sd;
		typedef mt_tag sc_category;
	};
	
	template <dim_t k, dim_t sd,
						dim_t d,
	          template <dim_t> class Data>
	list<ms<lmt<d,sd,Data>,k,Data> >& 
	simplex_container(lmt<d,sd,Data>& t) {
		return t.lmt_rep<k,sd,d,Data>::a; 
	}

	template <dim_t k, dim_t sd,
						dim_t d,
	          template <dim_t> class Data>
	const list<ms<lmt<d,sd,Data>,k,Data> >& 
	simplex_container(const lmt<d,sd,Data>& t) {
		return t.lmt_rep<k,sd,d,Data>::a; 
	}

	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	simplex_descriptor<lmt<d,sd,Data>,k>
	add(lmt<d,sd,Data>& t, simplex_descriptor<lmt<d,sd,Data>,k>& si) {
		ms<lmt<d,sd,Data>,k,Data> s;
		simplex_container<k,sd,d,Data>(t).push_front(s);
		si=simplex_descriptor<lmt<d,sd,Data>,k>(simplex_container<k,sd,d,Data>(t).begin());
		return si;
	}

	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	ms<lmt<d,sd,Data>,k,Data> *
	attr(lmt<d,sd,Data>& t, simplex_descriptor<lmt<d,sd,Data>,k> s) {
		return &*s.desc;
	}

	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	const ms<lmt<d,sd,Data>,k,Data> *
	attr(const lmt<d,sd,Data>& t, simplex_descriptor<lmt<d,sd,Data>,k> s) {
		return &*s.desc;
	}


	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	struct simplex_iterator<lmt<d,sd,Data>, k> {
		simplex_descriptor<lmt<d,sd,Data>, k> p;
		simplex_iterator(typename list<ms<lmt<d,sd,Data>,k,Data> >::iterator i)
						: p(i) {}
		simplex_iterator() {}
		simplex_iterator& operator++() {
			++p.desc;
			return *this;
		}
		simplex_descriptor<lmt<d,sd,Data>, k> operator*() {
			return p;
		}
		bool operator!=(const simplex_iterator& i) {
			return p!=i.p;
		}
		bool operator==(const simplex_iterator& i) {
			return p==i.p;
		}
	};

	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline
	simplex_descriptor<lmt<d,sd,Data>,k-1>
	face_op(const lmt<d,sd,Data>& t,
					simplex_descriptor<lmt<d,sd,Data>,k> s,
	        dim_t i) {
		return attr(t,s)->fop[i]; 
	}
	
	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	simplices(const lmt<d,sd,Data>& t,
						simplex_iterator<lmt<d,sd,Data>,k>& begin,
						simplex_iterator<lmt<d,sd,Data>,k>& end) {
		simplex_descriptor<lmt<d,sd,Data>,k> b(const_cast<list<ms<lmt<d,sd,Data>,k,Data> >&>(simplex_container<k,sd,d,Data>(t)).begin()); 
		simplex_descriptor<lmt<d,sd,Data>,k> e(const_cast<list<ms<lmt<d,sd,Data>,k,Data> >&>(simplex_container<k,sd,d,Data>(t)).end()); 
		begin.p=b;
		end.p=e;
	}

	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	simplices(lmt<d,sd,Data>& t,
						simplex_iterator<lmt<d,sd,Data>,k>& begin,
						simplex_iterator<lmt<d,sd,Data>,k>& end) {
		simplex_descriptor<lmt<d,sd,Data>,k> b(simplex_container<k,sd,d,Data>(t).begin()); 
		simplex_descriptor<lmt<d,sd,Data>,k> e(simplex_container<k,sd,d,Data>(t).end()); 
		begin.p=b;
		end.p=e;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline 
	std::pair<simplex_descriptor<lmt<d,sd,Data>,d>,
						simplex_descriptor<lmt<d,sd,Data>,d> >
	cells(const lmt<d,sd,Data>& t,
				simplex_descriptor<lmt<d,sd,Data>,d-1> s) {
		return std::make_pair(attr(t,s)->cov[0], attr(t,s)->cov[1]);
	}
	
	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline 
	typename disable_if<equal_c<k+1,d>,
		simplex_descriptor<lmt<d,sd,Data>,k+1> >::type
	up_simplex(const lmt<d,sd,Data>& t,
						 simplex_descriptor<lmt<d,sd,Data>,k> s) {
		return attr(t,s)->up;
	}
	
	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline 
	typename enable_if<equal_c<k+1,d>,
		simplex_descriptor<lmt<d,sd,Data>,k+1> >::type
	up_simplex(const lmt<d,sd,Data>& t,
						 simplex_descriptor<lmt<d,sd,Data>,k> s) {
		return attr(t,s)->cov[0];
	}
	
  //add
	
	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline
	void face_op_set(lmt<d,sd,Data>& t,
						 			 simplex_descriptor<lmt<d,sd,Data>,k> s,
									 dim_t i,
						 			 simplex_descriptor<lmt<d,sd,Data>,k-1> r) {
		attr(t,s)->fop[i]=r; 
	}
	
	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline 
	void
	cells_set(lmt<d,sd,Data>& t,
						simplex_descriptor<lmt<d,sd,Data>,d-1> s,
						simplex_descriptor<lmt<d,sd,Data>,d> r0,
						simplex_descriptor<lmt<d,sd,Data>,d> r1) {
		attr(t,s)->cov[0]=r0;
		attr(t,s)->cov[1]=r1;
	}
	
	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	up_simplex_set(lmt<d,sd,Data>& t,
						     simplex_descriptor<lmt<d,sd,Data>,k> s,
								 simplex_descriptor<lmt<d,sd,Data>,k+1> r) {
		attr(t,s)->up=r;
	}
						
	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	int orientation(const lmt<d,sd,Data>& t,
						    simplex_descriptor<lmt<d,sd,Data>,d> s) {
		return attr(t,s)->ori ? -1 : 1;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void orientation_set(lmt<d,sd,Data>& t,
						  	  		 simplex_descriptor<lmt<d,sd,Data>,d> s, int o) {
		attr(t,s)->ori= (o<0) ? true : false;
	}
						
	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	mark_type mark(const lmt<d,sd,Data>& t,
					 simplex_descriptor<lmt<d,sd,Data>,d> s) {
		return attr(t,s)->mk;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void mark_set(const lmt<d,sd,Data>& t,
						 	  simplex_descriptor<lmt<d,sd,Data>,d> s, mark_type i) {
		attr(const_cast<lmt<d,sd,Data>&>(t),s)->mk=i;
	}

	
	//del
	
	template <dim_t k, dim_t sd, dim_t d,
	          template <dim_t> class Data>
	inline
	void del(lmt<d,sd,Data>& t,
					 simplex_descriptor<lmt<d,sd,Data>,k> s) {
		simplex_container<k,sd,d,Data>(t).erase(s.desc);
	}


	//lmt


	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	bool is_current(const lmt<d,sd,Data>& t,
						  	  simplex_descriptor<lmt<d,sd,Data>,d> s) {
		return attr(t,s)->cur;
	}

	template <dim_t k, dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	typename disable_if<cell_c<lmt<d,sd,Data>,k>,bool>::type
	is_current(const lmt<d,sd,Data>& t,
		       	 simplex_descriptor<lmt<d,sd,Data>,k> s) {
		return is_current(t,top_cell(t,s));
	}
	
	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void is_current_set(lmt<d,sd,Data>& t,
						  	  		simplex_descriptor<lmt<d,sd,Data>,d> s, bool b) {
		attr(t,s)->cur=b;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	simplex_descriptor<lmt<d,sd,Data>,d>
	child(const lmt<d,sd,Data>& t,
	  	  simplex_descriptor<lmt<d,sd,Data>,d> s,int i) {
		return attr(t,s)->child[i];
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void
	child_set(lmt<d,sd,Data>& t,
	  	  	  simplex_descriptor<lmt<d,sd,Data>,d> s,
	  	 		 	int i,
	  	 		 	simplex_descriptor<lmt<d,sd,Data>,d> s1) {
		attr(t,s)->child[i]=s1;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	simplex_descriptor<lmt<d,sd,Data>,d>
	parent(const lmt<d,sd,Data>& t,
	  	 	 simplex_descriptor<lmt<d,sd,Data>,d> s) {
		return attr(t,s)->parent;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void
	parent_set(lmt<d,sd,Data>& t,
	  	 	 		 simplex_descriptor<lmt<d,sd,Data>,d> s,
	  	 	 		 simplex_descriptor<lmt<d,sd,Data>,d> p) {
		attr(t,s)->parent=p;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	simplex_descriptor<lmt<d,sd,Data>,sd>
	split_simplex(const lmt<d,sd,Data>& t,
	  	 			 simplex_descriptor<lmt<d,sd,Data>,0> s) {
		return attr(t,s)->split;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void
	split_simplex_set(lmt<d,sd,Data>& t,
	  	 			 		 simplex_descriptor<lmt<d,sd,Data>,0> s,
								 simplex_descriptor<lmt<d,sd,Data>,sd> e) {
		attr(t,s)->split=e;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	simplex_descriptor<lmt<d,sd,Data>,0>
	split_vertex(const lmt<d,sd,Data>& t,
	  	 			 	 simplex_descriptor<lmt<d,sd,Data>,sd> s) {
		return attr(t,s)->split;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void
	split_vertex_set(lmt<d,sd,Data>& t,
	  	 			 	 		 simplex_descriptor<lmt<d,sd,Data>,sd> s,
									 simplex_descriptor<lmt<d,sd,Data>,0> v) {
		attr(t,s)->split=v;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	int
	level(const lmt<d,sd,Data>& t,
	  	  simplex_descriptor<lmt<d,sd,Data>,d> s) {
		return (int)attr(t,s)->lvl;
	}

	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	inline
	void
	level_set(lmt<d,sd,Data>& t,
	  	 			simplex_descriptor<lmt<d,sd,Data>,d> s, int i) {
		attr(t,s)->lvl=(short)i;
	}
									
	/*! @} */

}

#endif // VGTL_MODEL_LMT_HPP
