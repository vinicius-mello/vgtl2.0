#ifndef VGTL_RIVARA_HPP
#define VGTL_RIVARA_HPP

/*! \file
 * \brief Rivara Scheme
 */

#include <algorithm>
#include <vgtl/top/maubach.hpp>

namespace vgtl {

	using std::list;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	//! 
	template <class T>
	void 
	build_lepp(const T& t, list<pair<Cell(T),Edge(T)> >& lp) {
		pair<Cell(T),Edge(T)> pce=lp.front();
		if(empty(pce.second)) {//find initial edge
			lp.pop_front();
			Edge(T) ea[3];
			double la[3];
			for(int i=0;i<3;++i) {
				ea[i]=face_op(t,pce.first,i);
				la[i]=length(t,ea[i]);
			}
			//sort edges
ilikegoto: //sorry...
			if(la[0]<=la[1]) {
				if(la[1]>la[2]) {
					double temp=la[2];la[2]=la[1];la[1]=temp;
					swap(ea[1],ea[2]);
					goto ilikegoto;
				}
			} else {
				double temp=la[1];la[1]=la[0];la[0]=temp;
				swap(ea[0],ea[1]);
				goto ilikegoto;
			}
			if(la[0]==la[2]) {//equilateral
				list<pair<Cell(T),Edge(T)> > lpa[3];
				for(int i=0;i<3;++i) {
					lpa[i].push_front(make_pair(pce.first,ea[i]));
					build_lepp(t,lpa[i]);
				}
				if(lpa[0].size()<=lpa[1].size()) {
					if(lpa[0].size()<=lpa[2].size()) 
						lp.splice(lp.begin(),lpa[0]);
					else
						lp.splice(lp.begin(),lpa[2]);
				} else {
					if(lpa[1].size()<=lpa[2].size()) 
						lp.splice(lp.begin(),lpa[1]);
					else
						lp.splice(lp.begin(),lpa[2]);
				}
				return;
			} else if(la[1]==la[2]) {//isosceles
				list<pair<Cell(T),Edge(T)> > lpa[2];
				for(int i=0;i<2;++i) {
					lpa[i].push_front(make_pair(pce.first,ea[i+1]));
					build_lepp(t,lpa[i]);
				}
				if(lpa[0].size()<=lpa[1].size()) 
					lp.splice(lp.begin(),lpa[0]);
				else
					lp.splice(lp.begin(),lpa[1]);
				return;
			}
			pce.second=ea[2];
			lp.push_front(pce);
		}
		while(true) {
			Cell(T) ceo=adjacent(t,pce.first,pce.second);
			if(ceo==pce.first) {//boundary
				Edge(T) e;
				empty_set(e);
				lp.push_front(make_pair(ceo,e));
				return;
			}
			Edge(T) ea[3];
			double la[3];
			int j;
			for(int i=0;i<3;++i) {
				ea[i]=face_op(t,ceo,i);
				la[i]=length(t,ea[i]);
				if(ea[i]==pce.second) j=i;
			}
			if(j!=0) {
				swap(ea[0],ea[j]);
				double temp=la[j];la[j]=la[0];la[0]=temp;
			}
			if((la[0]>=la[1])&&(la[0]>=la[2])) {//common
				Edge(T) e;
				empty_set(e);
				lp.push_front(make_pair(ceo,e));
				return;
			}
			if(la[1]==la[2]) {//isosceles
				list<pair<Cell(T),Edge(T)> > lpa[2];
				for(int i=0;i<2;++i) {
					lpa[i].push_front(make_pair(ceo,ea[i+1]));
					build_lepp(t,lpa[i]);
				}
				if(lpa[0].size()<=lpa[1].size()) 
					lp.splice(lp.begin(),lpa[0]);
				else
					lp.splice(lp.begin(),lpa[1]);
				return;
			} 
			// follow path
			if(la[1]>la[2]) 
				pce.second=ea[1];
			else 
				pce.second=ea[2];
			pce.first=ceo;
			lp.push_front(pce);
		}
	}

	//! 
	template <class T, class ApplyNew, class ApplyOld>
	void
	lepp_subdivide(T& t, list<pair<Cell(T),Edge(T)> >& lp, ApplyNew& app, ApplyOld& appold)  {
		while(true) {
			Cell(T) ca[2];
			pair<Cell(T),Edge(T)> pce=lp.front();
			lp.pop_front();
			ca[0]=pce.first;
			pce=lp.front();
			lp.pop_front();
			Edge(T) e=pce.second;
			ca[1]=pce.first;
			maubach_edge_split(t,e,0,ca,ca+((ca[0]==ca[1])?1:2),app,appold);
			if(lp.empty()) break;
			else build_lepp(t,lp);
		}
	}

	template <class T, class ApplyNew, class ApplyOld>
	void
	rivara_subdivide(T& t, Cell(T) ce, ApplyNew& app, ApplyOld& appold)  {
		Edge(T) e;
		empty_set(e);
		list<pair<Cell(T),Edge(T)> > lp;
		lp.push_front(make_pair(ce,e));
		build_lepp(t,lp);
		lepp_subdivide(t,lp,app,appold);
	}

	template <class T, class ApplyNew>
	void
	rivara_subdivide(T& t, Cell(T) ce, ApplyNew& app)  {
		do_nothing dn;
		rivara_subdivide(t,ce,app,dn);
	}

	template <class T, class ApplyNew, class ApplyOld>
	void
	fourtle_subdivide(T& t, Cell(T) ce, ApplyNew& app, ApplyOld& appold)  {
		Edge(T) e;
		empty_set(e);
		list<pair<Cell(T),Edge(T)> > lp;
		lp.push_front(make_pair(ce,e));
		build_lepp(t,lp);
		lepp_subdivide(t,lp,app,appold);
		Cell(T) ch[2];
		for(int i=0;i<2;++i) {
			ch[i]=child(t,ce,i);
			lp.push_front(make_pair(ch[i],face_op(t,ch[i],0)));
			build_lepp(t,lp);
			lepp_subdivide(t,lp,app,appold);
		}
	}

	template <class T, class ApplyNew>
	void
	fourtle_subdivide(T& t, Cell(T) ce, ApplyNew& app)  {
		do_nothing dn;
		fourtle_subdivide(t,ce,app,dn);
	}

	template <class T>
	Vertex(T)
	rivara_weld_vertex(const T& t, Cell(T) ce) {
		dim_t weld_ind=0;
		return face_ind(t,ce,_(weld_ind));
	}
	
	template <class T, class ApplyOld>
	void
	rivara_weld(T& t, Cell(T) ce, ApplyOld& app)  {
		Vertex(T) v=rivara_weld_vertex(t,ce);
		recursive_weld(t,v,rivara_weld_vertex<T>,app);
	}

	/*! @} */

}

#endif // VGTL_RIVARA_HPP
