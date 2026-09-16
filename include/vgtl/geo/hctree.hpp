#ifndef VGTL_GEO_HCTREE_HPP
#define VGTL_GEO_HCTREE_HPP

#include <vector>
#include <algorithm>
#include <vgtl/alg/point.hpp>

/*! \file
 * \brief 
 */


namespace vgtl {

	using namespace std;
	
	/*! \addtogroup geo 
	 * @{
	 */

	struct hctree_nothing {};

  template <dim_t n, class Scalar=double, class Data=hctree_nothing>
	struct hctree_node : public Data {
		short lvl;
		point<n,Scalar> center;
		hctree_node<n,Scalar,Data> * children[1<<n];
		hctree_node<n,Scalar,Data> * parent;
		hctree_node() {
			lvl=0;
			parent=0;
      for(dim_t i=0;i<(1<<n);++i) children[i]=0;
      for(dim_t i=0;i<n;++i) {
				center[i]=1.0;center[i]=center[i]/2;
			}
		}
		~hctree_node() {
      for(dim_t i=0;i<(1<<n);++i) if(children[i]) delete children[i];
		}
	};

  template <dim_t n, class Scalar, class Data>
	int level(hctree_node<n,Scalar,Data> * node) {
		return static_cast<int>(node->lvl);
	}

  template <dim_t n, class Scalar, class Data>
	Scalar diameter(hctree_node<n,Scalar,Data> * node) {
		int l=level(node);
		Scalar d=1.0;
		d=d/(1<<l);
		return d;
	}

  template <dim_t n, class Scalar, class Data>
	void level_set(hctree_node<n,Scalar,Data> * node, int lvl) {
		node->lvl=(short)lvl;
	}

  template <dim_t n, class Scalar, class Data>
	point<n,Scalar> center(hctree_node<n,Scalar,Data> * node) {
		return node->center;
	}

  template <dim_t n, class Scalar, class Data>
	void center_set(hctree_node<n,Scalar,Data> * node,
		const point<n,Scalar>& p) {
		node->center=p;
	}

  template <dim_t n, class Scalar, class Data>
	bool has_parent(hctree_node<n,Scalar,Data> * node) {
		return node->lvl!=0;
	}

  template <dim_t n, class Scalar, class Data>
	hctree_node<n,Scalar,Data> *
	parent(hctree_node<n,Scalar,Data> * node) {
		return node->parent;
	}

  template <dim_t n, class Scalar, class Data>
	void
	parent_set(hctree_node<n,Scalar,Data> * node,
		hctree_node<n,Scalar,Data> * nodep) {
		node->parent=nodep;
	}

  template <dim_t n, class Scalar, class Data>
	bool has_children(hctree_node<n,Scalar,Data> * node) {
		return node->children[0]!=0;
	}

  template <dim_t n, class Scalar, class Data>
	hctree_node<n,Scalar,Data> *
	child(hctree_node<n,Scalar,Data> * node, int i) {
		return node->children[i];
	}

  template <dim_t n, class Scalar, class Data>
	void
	child_set(hctree_node<n,Scalar,Data> * node, int i,
		hctree_node<n,Scalar,Data> * nodec) {
		node->children[i]=nodec;
	}

  template <dim_t n, class Scalar, class Data>
	void
	split(hctree_node<n,Scalar,Data> * node) {
		//array<int,n> v;
		point<n,Scalar> p;
		int l=level(node);
		Scalar d=diameter(node);
    for(dim_t c=0;c<(1<<n);++c) {
			p=center(node);
			hctree_node<n,Scalar,Data> * nodec=
				new hctree_node<n,Scalar,Data>();
			level_set(nodec,l+1);
			parent_set(nodec,node);
			child_set(node,c,nodec);
      for(dim_t i=0;i<n;++i) {
				p[i]+= ((1<<i)&c) ? (-d/4) : (d/4);
			}
			center_set(nodec,p);
		}
	}

  template <dim_t n, class Scalar, class Data>
	hctree_node<n,Scalar,Data> * 
	locate_node(hctree_node<n,Scalar,Data> * node, const point<n,Scalar>& p) {
		if(!has_children(node)) return node;
		point<n,Scalar> q;
		q=p-center(node);
		int c=0;
    for(dim_t i=0;i<n;++i) {
			if(q[i]<0) c=((1<<i)|c);
		}
		return locate_node(child(node,c),p);
	}

	/*! @} */

}

#endif // VGTL_GEO_HCTREE_HPP
