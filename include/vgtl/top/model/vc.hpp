#ifndef VGTL_VC_HPP
#define VGTL_VC_HPP

/*! \file
 * \brief Models a APM using vectors
 */

#include <vector>
#include <vgtl/utl/array.hpp>
#include <vgtl/top/incidence.hpp>
#include <vgtl/top/add.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/top_cell.hpp>
#include <vgtl/top/model/bs.hpp>

namespace vgtl {

	using std::vector;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <dim_t d,
						template <dim_t> class Data>
	struct vc;

	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	struct simplex_descriptor<vc<d,Data>,k> {
		int desc;
		simplex_descriptor& operator=(const simplex_descriptor& _s) {
			desc=_s.desc;
			return (*this);
		}
		bool operator<(const simplex_descriptor& _s) const {
			return desc<_s.desc;
		}
		bool operator==(const simplex_descriptor& _s) const {
			return desc==_s.desc;
		}
		bool operator!=(const simplex_descriptor& _s) const {
			return desc!=_s.desc;
		}
		simplex_descriptor() : desc(-1) {}
		simplex_descriptor(int _desc) : desc(_desc) {}
	};

	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	inline
	bool
	empty(simplex_descriptor<vc<d,Data>,k> s) {
	  return s.desc==-1;
	}

	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	inline
	void
	empty_set(simplex_descriptor<vc<d,Data>,k>& s) {
	  s.desc=-1;
	}

	template <dim_t d, 
	          template <dim_t> class Data>
	struct sc_traits<vc<d,Data> > {
		static const dim_t dim=d;
		static const dim_t split_dim=0;
		typedef apm_tag sc_category; 
	};

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	struct vc_rep : vc_rep <k-1,d,Data> {
		vector<bs<vc<d,Data>,k,Data> > a;
	};

	template <dim_t d,
	          template <dim_t> class Data>
	struct vc_rep<0,d,Data> {
		vector<bs<vc<d,Data>,0,Data> > a;
	};

	//! Models a Mutable APM
	template <dim_t d,
						template <dim_t> class Data=extra_data>
	struct vc : vc_rep<d,d,Data>, Data<d+1> {
		static const dim_t dim=d;
		typedef apm_tag sc_category;
	};
					
	template <dim_t k,
						dim_t d,
	          template <dim_t> class Data>
	inline
	const vector<bs<vc<d,Data>,k,Data> >& 
	simplex_container(const vc<d,Data>& t) {
		return t.vc_rep<k,d,Data>::a; 
	}

	template <dim_t k,
						dim_t d,
	          template <dim_t> class Data>
	inline
	vector<bs<vc<d,Data>,k,Data> >& 
	simplex_container(vc<d,Data>& t) {
		return t.vc_rep<k,d,Data>::a; 
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	simplex_descriptor<vc<d,Data>,k>
	add(vc<d,Data>& t, simplex_descriptor<vc<d,Data>,k>& sd) {
		bs<vc<d,Data>,k,Data> s;
		sd=simplex_descriptor<vc<d,Data>,k>(simplex_container<k,d,Data>(t).size());
		simplex_container<k,d,Data>(t).push_back(s);
		return sd;
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	bs<vc<d,Data>,k,Data> *
	attr(vc<d,Data>& t, simplex_descriptor<vc<d,Data>,k> s) {
		return &(simplex_container<k,d,Data>(t)[s.desc]);
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	const bs<vc<d,Data>,k,Data> *
	attr(const vc<d,Data>& t, simplex_descriptor<vc<d,Data>,k> s) {
		return &(simplex_container<k,d,Data>(t)[s.desc]);
	}


	template <dim_t d, dim_t k,
	          template <dim_t> class Data>
	struct simplex_iterator<vc<d,Data>, k> {
		simplex_descriptor<vc<d,Data>, k> p;
		simplex_iterator(int i) : p(i) {}
		simplex_iterator() {}
		simplex_iterator& operator++() {
			++p.desc;
			return *this;
		}
		simplex_descriptor<vc<d,Data>, k> operator*() {
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
	simplex_descriptor<vc<d,Data>,k-1>
	face_op(const vc<d,Data>& t,
					simplex_descriptor<vc<d,Data>,k> s,
	        dim_t i) {
		return attr(t,s)->fop[i]; 
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	simplices(const vc<d,Data>& t,
						simplex_iterator<vc<d,Data>,k>& begin,
						simplex_iterator<vc<d,Data>,k>& end) {
		simplex_descriptor<vc<d,Data>,k> b(0);
		simplex_descriptor<vc<d,Data>,k> e(simplex_container<k,d,Data>(t).size()); 
		begin.p=b;
		end.p=e;
	}

	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	simplices(vc<d,Data>& t,
						simplex_iterator<vc<d,Data>,k>& begin,
						simplex_iterator<vc<d,Data>,k>& end) {
		
		simplex_descriptor<vc<d,Data>,k> b(0);
		simplex_descriptor<vc<d,Data>,k> e(simplex_container<k,d,Data>(t).size()); 
		begin.p=b;
		end.p=e;
	}

	template <dim_t d,
	          template <dim_t> class Data>
	inline 
	std::pair<simplex_descriptor<vc<d,Data>,d>,
						simplex_descriptor<vc<d,Data>,d> >
	cells(const vc<d,Data>& t,
				simplex_descriptor<vc<d,Data>,d-1> s) {
		return std::make_pair(attr(t,s)->cov[0], attr(t,s)->cov[1]);
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline 
	typename disable_if<equal_c<k+1,d>,
		simplex_descriptor<vc<d,Data>,k+1> >::type
	up_simplex(const vc<d,Data>& t,
						 simplex_descriptor<vc<d,Data>,k> s) {
		return attr(t,s)->up;
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline 
	typename enable_if<equal_c<k+1,d>,
		simplex_descriptor<vc<d,Data>,k+1> >::type
	up_simplex(const vc<d,Data>& t,
						 simplex_descriptor<vc<d,Data>,k> s) {
		return attr(t,s)->cov[0];
	}
	
  //add
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void face_op_set(vc<d,Data>& t,
						 			 simplex_descriptor<vc<d,Data>,k> s,
									 dim_t i,
						 			 simplex_descriptor<vc<d,Data>,k-1> r) {
		attr(t,s)->fop[i]=r; 
	}
	
	template <dim_t d,
	          template <dim_t> class Data>
	inline 
	void
	cells_set(vc<d,Data>& t,
						simplex_descriptor<vc<d,Data>,d-1> s,
						simplex_descriptor<vc<d,Data>,d> r0,
						simplex_descriptor<vc<d,Data>,d> r1) {
		attr(t,s)->cov[0]=r0;
		attr(t,s)->cov[1]=r1;
	}
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void
	up_simplex_set(vc<d,Data>& t,
						     simplex_descriptor<vc<d,Data>,k> s,
								 simplex_descriptor<vc<d,Data>,k+1> r) {
		attr(t,s)->up=r;
	}
							
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	int orientation(const vc<d,Data>& t,
						   	  simplex_descriptor<vc<d,Data>,d> s) {
		return attr(t,s)->ori ? -1 : 1;
	}

	template <dim_t d, 
	          template <dim_t> class Data>
	inline
	void orientation_set(vc<d,Data>& t,
						  	  		 simplex_descriptor<vc<d,Data>,d> s, int o) {
		attr(t,s)->ori= (o<0) ? true : false;
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	mark_type mark(const vc<d,Data>& t,
					 simplex_descriptor<vc<d,Data>,d> s) {
		return attr(t,s)->mk;
	}

	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void mark_set(const vc<d,Data>& t,
						 	  simplex_descriptor<vc<d,Data>,d> s, mark_type i) {
		attr(const_cast<vc<d,Data>&>(t),s)->mk=i;
	}


	template <dim_t d,
	          template <dim_t> class Data>
	int
	level(const vc<d,Data>& t,
	  	  simplex_descriptor<vc<d,Data>,d> s) {
		return (int)attr(t,s)->lvl;
	}

	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void
	level_set(vc<d,Data>& t,
	  	 			simplex_descriptor<vc<d,Data>,d> s, int i) {
		attr(t,s)->lvl=(short)i;
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	bool is_current(const vc<d,Data>& t,
						  	  simplex_descriptor<vc<d,Data>,d> s) {
		return attr(t,s)->cur;
	}

	template <dim_t k, dim_t d, 
	          template <dim_t> class Data>
	inline
	typename disable_if<cell_c<vc<d,Data>,k>,bool>::type
	is_current(const vc<d,Data>& t,
		       	 simplex_descriptor<vc<d,Data>,k> s) {
		return true;
	}
	
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void is_current_set(vc<d,Data>& t,
						  	  		simplex_descriptor<vc<d,Data>,d> s, bool b) {
		attr(t,s)->cur=b;
	}

	//del
	
	template <dim_t k, dim_t d,
	          template <dim_t> class Data>
	inline
	void del(vc<d,Data>& t,
					 simplex_descriptor<vc<d,Data>,k> s) {
	}
						
	template <dim_t d,
	          template <dim_t> class Data>
	inline
	void del(vc<d,Data>& t,
					 simplex_descriptor<vc<d,Data>,d> s) {
		is_current_set(t,s,false);
	}
						
	/*! @} */

}

#endif // VGTL_VC_HPP
