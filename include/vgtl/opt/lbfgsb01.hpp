#ifndef VGTL_LBFGSB01_HPP
#define VGTL_LBFGSB01_HPP

#include <cstring>
#include "lbfgsb01.h"

namespace vgtl {

	class lbfgsb01 {
		int max_n;
		int max_m;
		double * wa;
		int * iwa;
		char task[60];
		int iprint;
		char csave[60];
		int lsave[4];
		int isave[44];
		double dsave[29];
		double lower,upper;
		double factr;
		double pgtol;
		int max_it;
		public:
		lbfgsb01(int _n, int _m) : max_n(_n), max_m(_m) {
			wa=new double[2*max_m*max_n+4*max_n+11*max_m*max_m+8*max_m];
			iwa=new int[3*max_n];
			lower=0;
			upper=1;
			factr=1.0e+4;
			pgtol=1.0e-5;
			iprint=-1;
			max_it=300;
		}
		~lbfgsb01() {
			delete [] wa;
			delete [] iwa;
		}

		void max_iterations_set(int i) {
			max_it=i;
		}
		void factr_set(double d) {
			factr=d;
		}
		void pgtol_set(double d) {
			pgtol=d;
		}
		void print_set(bool p) {
			if(p) iprint=1;
			else iprint=-1;
		}
		template <class Driver>
		double optimize(Driver& drv) {
			drv.init();
			int n=drv.n();
			int m=drv.m();
			if((n>max_n)||(m>max_m)) {
				max_n=n;
				delete [] iwa;
				iwa=new int[3*max_n];
				max_m=m;
				delete [] wa;
				wa=new double[2*max_m*max_n+4*max_n+11*max_m*max_m+8*max_m];
			}
			strcpy(task,"START");
			memset(task+5,' ',55);
			double f;
			for(int i=0;i<max_it;++i) {
				setb01_(&n,&m,drv.x(),&lower,&upper,&f,drv.grad(),
					&factr,&pgtol,wa,
					iwa, task, &iprint, csave,
					lsave,isave,dsave);
				if(!strncmp("FG",task,2)) {
					f=drv.f();
					drv.update_grad();
				} else if(!strncmp("NEW_X",task,5)) {
				} else {
					break;
				}
			}
			drv.reset();
			return f;
		}
	};
}

#endif // VGTL_LBFGSB01_HPP
