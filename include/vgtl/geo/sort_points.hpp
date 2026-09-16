#ifndef VGTL_SORT_POINTS_HPP
#define VGTL_SORT_POINTS_HPP

#include <vgtl/alg/barychange.hpp>
#include <vgtl/geo/inside.hpp>
#include <vgtl/utl/pair_tie.hpp>
#include <vector>

namespace vgtl {

	using std::vector;

	template <class T>
	struct points_cmp {
		const T& t;
		static const dim_t d=Dim(T);
		points_cmp(const T& _t) : t(_t) {}
		bool operator()(const point<d,double>& p0, const point<d,double>& p1) const {
			return locate(t,p0)<locate(t,p1);
		}
	};
				
	template <class T>
	struct point_cmp1 {
		const T& t;
		static const dim_t d=Dim(T);
		point_cmp1(const T& _t) : t(_t) {}
		bool operator()(const point<d,double>& p,
										Cell(T) cv) const {
			return locate(t,p)<cv;
		}
	};
				
	template <class T>
	struct point_cmp2 {
		const T& t;
		static const dim_t d=Dim(T);
		point_cmp2(const T& _t) : t(_t) {}
		bool operator()(Cell(T) cv,
										const point<d,double>& p) const {
			return cv<locate(t,p);
		}
	};
	
	template <class T>
	void classify_points(T& t) {
		sort(points(t).begin(),points(t).end(),points_cmp<T>(t));
		for(typename vector<Cell(T)>::iterator  i=initial_cells(t).begin();
			i!=initial_cells(t).end(); ++i) {
			Cell(T) cv=*i;
			points_begin_set(t,cv,
				lower_bound(points(t).begin(),points(t).end(),cv,point_cmp1<T>(t)));
			points_end_set(t,cv,
				upper_bound(points(t).begin(),points(t).end(),cv,point_cmp2<T>(t)));
		}
	}
				
	template <dim_t d>
	struct point_in_cmp {
		barychange<d,double> bch;
		point_in_cmp(const array<point<d,double>,d+1>& ps) : bch(ps) {}
		bool operator()(const point<d,double>& p) {
			point<d+1,double> bp=bch(p);
			bool t=true;
			for(int i=0;i<=d;++i) if(bp[i]<0) {t=false;break;}
			return t;
		}
	};
	
	template <class T>
	void test_points(T& t, Cell(T) cv) {
		array<point<Dim(T),double>,Dim(T)+1> ps; // rever
		euclidean_points(t,cv,ps);
		barychange<Dim(T),double> bch(ps);
		for(typename vector<point<Dim(T),double> >::iterator i=points_begin(t,cv);
				i!=points_end(t,cv);++i) {
			point<Dim(T)+1,double> bp=bch(*i);
			bool t=true;
			for(int j=0;j<=Dim(T);++j) if(bp[j]<0) {t=false;break;}
			using namespace std;
			if(!t) cout<<bp<<endl;
		}
	}

	template <class T>
	void classify_points(T& t, Cell(T) cv) {
		Cell(T) cv0,cv1;
		cv0=child(t,cv,0);
		cv1=child(t,cv,1);
		array<point<Dim(T),double>,Dim(T)+1> ps; // rever
		euclidean_points(t,cv0,ps);
		points_begin_set(t,cv0,points_begin(t,cv));
		points_end_set(t,cv0,
			partition(points_begin(t,cv),points_end(t,cv),
							point_in_cmp<Dim(T)>(ps)));
		points_begin_set(t,cv1,points_end(t,cv0));
		points_end_set(t,cv1,points_end(t,cv));
	}
	
	template <class T>
	void reclassify_points(T& t, Cell(T) cv) {
		if(!has_children(t,cv)) return;
		Cell(T) cv0,cv1;
		cv0=child(t,cv,0);
		cv1=child(t,cv,1);
		array<point<Dim(T),double>,Dim(T)+1> ps; // rever
		euclidean_points(t,cv0,ps);
		points_begin_set(t,cv0,points_begin(t,cv));
		points_end_set(t,cv0,
			partition(points_begin(t,cv),points_end(t,cv),
							point_in_cmp<Dim(T)>(ps)));
		points_begin_set(t,cv1,points_end(t,cv0));
		points_end_set(t,cv1,points_end(t,cv));
		reclassify_points(t,cv0);
		reclassify_points(t,cv1);
	}

	template <class T>
	void reclassify_points(T& t) {
		for(typename vector<Cell(T)>::iterator  i=initial_cells(t).begin();
			i!=initial_cells(t).end(); ++i) {
			Cell(T) cv=*i;
			reclassify_points(t,cv);
		}
	}

}

#endif // VGTL_SORT_POINTS_HPP
