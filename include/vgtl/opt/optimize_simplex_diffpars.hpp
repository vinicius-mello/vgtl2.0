#ifndef VGTL_OPT_OPTIMIZE_SIMPLEX_DIFFPARS_HPP
#define VGTL_OPT_OPTIMIZE_SIMPLEX_DIFFPARS_HPP

#include <iostream>
#include <vgtl/opt/lbfgsb01.hpp>
#include <vgtl/top/sc.hpp>

namespace vgtl {

	using namespace std;

	template <class T, dim_t k,
		class CompSimplexErr=double (*)(const T&, Simplex(T,k))>
	class optimize_simplex_error {
		double * base_x;
		double * base_grad;
		int n_;
		int m_;
		T& t;
		Simplex(T,k) s;
		int p;
		CompSimplexErr cserr;
		double eps;
		public:
		optimize_simplex_error(T& _t, Simplex(T,k) _s, int _p,
			CompSimplexErr _cserr, double * grad, double _eps=0.000001) : 
			eps(_eps), s(_s), t(_t), p(_p), cserr(_cserr), base_grad(grad) {
			n_=(k+1)*p;
			m_=(k+1);	
			base_x=diffpars(t,s);
		}
		double f() {
			return cserr(t,s);
		}
		void update_grad() {
			double err=f();
			double epsk=eps/(k+1);
			double epsk1=k*eps/(k+1);
			for(int l=0;l<p;++l) {
				double * grad=&base_grad[(k+1)*l];
				double * xp=&base_x[(k+1)*l];
				grad[k]=0;
				for(int i=0;i<k;++i) {
					for(int j=0;j<=k;++j) {
						if(j!=i) xp[j]-=epsk;
						else xp[j]+=epsk1;
					}		
					double err_eps=f();
					grad[i]=(err_eps-err)/eps;
					grad[k]-=grad[i];
					for(int j=0;j<=k;++j) {
						if(j!=i) xp[j]+=epsk;
						else xp[j]-=epsk1;
					}
				}
			}
		}
		int n() {
			return n_;
		}
		int m() {
			return m_;
		}
		double * x() {
			return base_x;
		}
		double * grad() {
			return base_grad;
		}
		void init() {
		}
		void reset() {
		}
		void eps_set(double d) {
			eps=d;
		}
	};

}

#endif // VGTL_OPT_OPTIMIZE_SIMPLEX_DIFFPARS_HPP
