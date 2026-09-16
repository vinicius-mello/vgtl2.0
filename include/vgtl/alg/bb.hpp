#ifndef VGTL_BB_HPP
#define VGTL_BB_HPP

/*! \file
 * \brief BB definitions
 */

#include <vgtl/alg/point.hpp>
#include <vgtl/comb/binomial.hpp>
#include <vgtl/comb/triangular.hpp>
#include <vgtl/utl/constraint.hpp>
#include <vector>

namespace vgtl {
		
	using std::vector;
	
	/*! \addtogroup alg 
	 * @{
	 */

  template <dim_t n>
	struct de_Casteljau {
		int operator()(int i, int j) const;
	};

  template <>
	struct de_Casteljau<1> {
		int operator()(int i, int j) const {
			return 0;
		}
	};
	
  template <>
	struct de_Casteljau<2> {
		int operator()(int i, int j) const {
			return i+j;
		}
	};
	
  template <>
	struct de_Casteljau<3> {
		int operator()(int i, int j) const {
			static const unsigned char lut[15][3]={
				{0 ,1 ,2},
				{1 ,3 ,4},
				{2 ,4 ,5},
				{3 ,6 ,7},
				{4 ,7 ,8},
				{5 ,8 ,9},
				{6 ,10 ,11},
				{7 ,11 ,12},
				{8 ,12 ,13},
				{9 ,13 ,14},
				{10 ,15 ,16},
				{11 ,16 ,17},
				{12 ,17 ,18},
				{13 ,18 ,19},
				{14 ,19 ,20}
			};
			if(i<15) {
				return lut[i][j];
			} else {
				int k=0;
				int ii=i+1;
				do {
					k++;
					for(int l=0;(l<k)&&(ii>0);++l) --ii;
				} while(ii>0);
				if(j==0) return i;
				if(j==1) return i+k;
				if(j==2) return i+k+1;
			}
		}
	};
	
  template <>
	struct de_Casteljau<4> {
		int operator()(int i, int j) const {
			static const unsigned char lut[56][4]={
				{0 ,1 ,2 ,3},
				{1 ,4 ,5 ,6},
				{2 ,5 ,7 ,8},
				{3 ,6 ,8 ,9},
				{4 ,10 ,11 ,12},
				{5 ,11 ,13 ,14},
				{6 ,12 ,14 ,15},
				{7 ,13 ,16 ,17},
				{8 ,14 ,17 ,18},
				{9 ,15 ,18 ,19},
				{10 ,20 ,21 ,22},
				{11 ,21 ,23 ,24},
				{12 ,22 ,24 ,25},
				{13 ,23 ,26 ,27},
				{14 ,24 ,27 ,28},
				{15 ,25 ,28 ,29},
				{16 ,26 ,30 ,31},
				{17 ,27 ,31 ,32},
				{18 ,28 ,32 ,33},
				{19 ,29 ,33 ,34},
				{20 ,35 ,36 ,37},
				{21 ,36 ,38 ,39},
				{22 ,37 ,39 ,40},
				{23 ,38 ,41 ,42},
				{24 ,39 ,42 ,43},
				{25 ,40 ,43 ,44},
				{26 ,41 ,45 ,46},
				{27 ,42 ,46 ,47},
				{28 ,43 ,47 ,48},
				{29 ,44 ,48 ,49},
				{30 ,45 ,50 ,51},
				{31 ,46 ,51 ,52},
				{32 ,47 ,52 ,53},
				{33 ,48 ,53 ,54},
				{34 ,49 ,54 ,55}, // Up to degree 5
				{35 ,56 ,57 ,58},
				{36 ,57 ,59 ,60},
				{37 ,58 ,60 ,61},
				{38 ,59 ,62 ,63},
				{39 ,60 ,63 ,64},
				{40 ,61 ,64 ,65},
				{41 ,62 ,66 ,67},
				{42 ,63 ,67 ,68},
				{43 ,64 ,68 ,69},
				{44 ,65 ,69 ,70},
				{45 ,66 ,71 ,72},
				{46 ,67 ,72 ,73},
				{47 ,68 ,73 ,74},
				{48 ,69 ,74 ,75},
				{49 ,70 ,75 ,76},
				{50 ,71 ,77 ,78},
				{51 ,72 ,78 ,79},
				{52 ,73 ,79 ,80},
				{53 ,74 ,80 ,81},
				{54 ,75 ,81 ,82},
				{55 ,76 ,82 ,83} // Up to degree 6
			};
			return lut[i][j];
		}
	};
	
  template <dim_t n>
	int hom2ind1(array<int,n>& c) {
		de_Casteljau<n> cj;
    dim_t i=0;
		for(;i<n;++i) if(c[i]!=0) break;
		if(i==n) return 0;
		c[i]-=1;
		int t=hom2ind1(c);
		return cj(t,i);
	}

  template <dim_t n>
  typename enable_if<greater_c<n,0>,int>::type
	hom2ind(const array<int,n>& c) {
		array<int,n> ct=c;
		return hom2ind1(ct);
	}

  template <dim_t n>
  typename disable_if<greater_c<n,0>,int>::type
	hom2ind(const array<int,n>& c) {
		return 0;
	}

  template <dim_t n>
	int degree(const vec<n,int>& m) {
		int r=0;
    for(dim_t i=0; i<n; ++i) r+=m[i];
		return r;
	}
	
  template <dim_t n, class Scalar>
	class bbform : public vector<point<n,Scalar> > {
		int m;
		public:
		bbform() : m(0) {}
		bbform(int _m) : m(_m) {
			reset();
		}	
		point<n,Scalar>& operator[](const vec<n,int>& v) {
			triangular_traverser<n,int> t(m);
			int j=0;
			do {
				vec<n,int> i;
				i=*t;
				if(i==v) return vector<point<n,Scalar> >::operator[](j);
				++j;
			} while(++t);
		}
		const point<n,Scalar>& operator[](const vec<n,int>& v) const {
			triangular_traverser<n,int> t(m);
			int j=0;
			do {
				vec<n,int> i;
				i=*t;
				if(i==v) return vector<point<n,Scalar> >::operator[](j);
				++j;
			} while(++t);
		}
		int degree() const {
			return m;
		}
		int degree_set(int _m) {
			return m=_m;
		}
		void reset() {	
			triangular_traverser<n,int> t(m);
			vector<point<n,Scalar> >::clear();
			Scalar f=1;
			f/=m;
			do {
				vec<n,int> i;
				i=*t;
				point<n,Scalar> r;
        for(dim_t j=0;j<n;++j) r[j]=f*i[j];
				this->push_back(r);
			} while(++t);
		}
		void resetq() {	
			triangular_traverser<n-1,int> t(m);
			vector<point<n,Scalar> >::clear();
			Scalar f=1;
			f/=m;
			do {
				vec<n-1,int> i;
				i=*t;
				point<n,Scalar> r;
        for(dim_t j=0;j<(n-1);++j) r[j]=f*i[j];
				r[n-1]=1;
				this->push_back(r);
			} while(++t);
		}
	};

  template <dim_t n, dim_t l, class Scalar>
	void eval_dc(const vector<point<l,Scalar> >& bb,
		int degree, const point<n,Scalar>& p,
		vector<point<l,Scalar> >& temp, bool homogeneous=false) {
		typename vector<point<l,Scalar> >::const_iterator cbb=bb.begin();
		array<point<l,Scalar>,n> pb;
		de_Casteljau<n> cj;
		int B=0;
		temp.clear();
		for(int m=degree-1;m>=0;--m) {
			int N=binomial(m+n-1,m);
			for(int j=0; j<N; ++j) {
        for(dim_t i=0;i<n;++i) {
					pb[i]=*(cbb+cj(j,i));
					if((B==0)&&(homogeneous)) {
            for(dim_t k=0;k<(l-1);++k) pb[i][k]*=pb[i][l-1];
					}
				}
				temp.push_back(barycentric_combination(p,pb));
			}
			cbb=temp.begin()+B;
			B=temp.size();
		}
	}

  template <dim_t n, dim_t l, class Scalar>
	void  diff_dc(const vector<point<l,Scalar> >& bb,
		int degree, const point<n,Scalar>& p,
		array<vec<n,Scalar>,l>& J, bool homogeneous=false) {
		const vector<point<l,Scalar> >& cbb=bb;
		point<l,Scalar> q;
		de_Casteljau<n> cj;
		int m=degree;
		int N=binomial(m-1+n-1,n-1);
		vector<point<l,Scalar> > temp, temp2; 
		temp.reserve(N);
		temp2.reserve(N);
    for(dim_t i=0;i<n;++i) {
			for(int j=0; j<N; ++j) {
				temp.push_back(cbb[cj(j,i)]);
			}
			eval_dc(temp,degree-1,p,temp2,homogeneous);
			q=temp2.back();
			temp.clear();
      for(dim_t j=0;j<l;++j) {
				J[j][i]=m*q[j];
			}
		}
	}

  template <dim_t n, dim_t l, class Scalar>
	point<l,Scalar> eval(const bbform<l,Scalar>& bb, 
											 const point<n,Scalar>& p) {
		const vector<point<l,Scalar> >& vbb=bb;
		vector<point<l,Scalar> > temp;
		temp.reserve(vbb.size());
		eval_dc(vbb,bb.degree(),p,temp,false);
		return temp.back();
	}

  template <dim_t n, dim_t l, class Scalar>
	point<l-1,Scalar> evalq(const bbform<l,Scalar>& bb, 
											 const point<n,Scalar>& p) {
		const vector<point<l,Scalar> >& vbb=bb;
		vector<point<l,Scalar> > temp;
		temp.reserve(vbb.size());
		eval_dc(vbb,bb.degree(),p,temp,true);
		point<l,Scalar> r=temp.back();
		point<l-1,Scalar> rr;
    for(dim_t i=0;i<(l-1);++i) rr[i]=r[i]/r[l-1];
		return rr;
	}

  template <dim_t n, dim_t l, class Scalar>
	point<l,Scalar> jacobian(const bbform<l,Scalar>& bb, 
													 const point<n,Scalar>& p,
				 									 array<vec<n,Scalar>,l>& J) {
		const vector<point<l,Scalar> >& cbb=bb;
		diff_dc(cbb,bb.degree(),p,J,false);
		return eval(bb,p);
	}

  template <dim_t n, dim_t l, class Scalar>
	point<l-1,Scalar> jacobianq(const bbform<l,Scalar>& bb, 
													 const point<n,Scalar>& p,
				 									 array<vec<n,Scalar>,l-1>& J) {
		const vector<point<l,Scalar> >& cbb=bb;
		array<vec<n,Scalar>,l> JJ;
		diff_dc(cbb,bb.degree(),p,JJ,true);
		vector<point<l,Scalar> > temp;
		temp.reserve(cbb.size());
		eval_dc(cbb,bb.degree(),p,temp,true);
		point<l,Scalar> r=temp.back();
    for(dim_t i=0;i<n;++i) {
      for(dim_t j=0;j<(l-1);++j) {
				J[j][i]=(JJ[j][i]-r[j]*JJ[l-1][i]/r[l-1])/r[l-1];
			}
		}
		point<l-1,Scalar> rr;
    for(dim_t i=0;i<(l-1);++i) rr[i]=r[i]/r[l-1];
		return rr;
	}

  template <dim_t n, dim_t k, class Scalar>
	point<n,Scalar> jacobian2(const bbform<n,Scalar>& bb, 
													 const point<n,Scalar>& p,
				 									 array<vec<n,Scalar>,k>& vs) {
		array<vec<n,Scalar>,k> vst;
		const vector<point<n,Scalar> >& cbb=bb;
		bbform<n,Scalar> temp; 
		point<n,Scalar> q;
		de_Casteljau<n> cj;
		int m=bb.degree();
		temp.degree_set(m-1);
		int N=binomial(m-1+n-1,n-1);
    for(dim_t i=0;i<n;++i) {
			for(int j=0; j<N; ++j) {
				temp.push_back(cbb[cj(j,i)]);
			}
			q=eval(temp,p);
			temp.clear();
      for(dim_t l=0;l<k;++l) {
        for(dim_t j=0;j<n;++j) {
					vst[l][j]+=m*q[j]*vs[l][i];
				}
			}
		}
    for(dim_t l=0;l<k;++l) {
			vs[l]=vst[l];
		}
		return eval(bb,p);
	}

  template <dim_t n, class Scalar>
	point<n,Scalar> snap(const vec<n,int>& m, const point<n,Scalar>& p) {
		vec<n,Scalar> t;
		Scalar s=0;
    for(dim_t i=0;i<n;++i) {
			if((p[i]<0)||(m[i]==0)) t[i]=0;
			else { t[i]=p[i]; s+=t[i];}
		}
		t*=1/s;
		return t;
	}

	/*! @} */

}

#endif //VGTL_BB_HPP

