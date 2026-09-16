#ifndef VGTL_OPT_OPTIMIZE_SCALAR_CB_HPP
#define VGTL_OPT_OPTIMIZE_SCALAR_CB_HPP

#include <iostream>
#include <list>
#include <vgtl/top/star.hpp>
#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/opt/lbfgsb01.hpp>


namespace vgtl {

	using std::list;
	using std::set;

	template <class T, class CompCellErr>
	double compute_scalar_error(complex_buffer<T>& cb, CompCellErr ccerr) {
		double r=0;
		int c=0;
		for(typename set<Cell(T)>::iterator cei=cb.template begin<Dim(T)>();
				cei!=cb.template end<Dim(T)>(); ++cei) {
			Cell(T) ce=*cei;
			double err=ccerr(cb.t,ce);
			scalar_error_set(cb.t,ce,err);
			r+=err;
			c++;
		}
		return r/c;
	}
	
	template <class T, class CompCellErr> 
	void compute_scalar_grad(complex_buffer<T>& cb, Vertex(T) v,
		double * grad, double eps, CompCellErr ccerr) {
		list<Cell(T)> st;
		star(cb.t,v,back_inserter(st));
		double err=0;
		int c=0;
		for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
			Cell(T) cv=*j;
			if(exists(cb,cv)) {
				err+=scalar_error(cb.t,cv);
				c++;
			}
		}
		if(c==0) {
			grad[0]=0;
			return;
		}
		double * xp=scalar_value_x(cb.t,v);
		xp[0]+=eps;
		double err_eps=0;
		for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
			Cell(T) cv=*j;
			if(exists(cb,cv))
				err_eps+=ccerr(cb.t,cv);
		}
		grad[0]=(err_eps-err)/c/eps;
		xp[0]-=eps;
	}
					
	
	template <class T, class CompCellErr>
	void
	compute_scalar_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		for(typename set<Simplex(T,0)>::iterator vi=cb.template begin<0>();
				vi!=cb.template end<0>();++vi) {
			Vertex(T) v=*vi;
			compute_scalar_grad(cb,v,base_grad+(scalar_value_x(cb.t,v)-base_x),eps,ccerr);
		}
	}
	
	template <class T>
	double *
	init_scalar(complex_buffer<T>& cb, double * base_x) {
		for(typename set<Simplex(T,0)>::iterator vi=cb.template begin<0>();
				vi!=cb.template end<0>();++vi) {
			Simplex(T,0) v=*vi;
			scalar_value_x_set(cb.t,v,base_x);
			base_x+=1;
		}
		return base_x;
	}
	
	template <class T>
	void
	reset_scalar(complex_buffer<T>& cb) {
		for(typename set<Simplex(T,0)>::iterator vi=cb.template begin<0>();
				vi!=cb.template end<0>();++vi) {
			Vertex(T) v=*vi;
			scalar_value_x_set(cb.t,v,0);
		}
	}

	template <class T, class CompCellErr=double (*)(const T&, Cell(T))>
	class optimize_scalar_error_cb {
		double * base_x;
		double * base_grad;
		int n_;
		int m_;
		complex_buffer<T>& cb;
		CompCellErr ccerr;
		double eps_;
		public:
		optimize_scalar_error_cb(complex_buffer<T>& _cb,
								CompCellErr _ccerr, int _m, double * base_x_, 
			double * base_grad_, double _eps=0.000001) : m_(_m), 
			eps_(_eps), base_x(base_x_), base_grad(base_grad_), cb(_cb),
			ccerr(_ccerr) {
		}
		double f() {
		//	cout<<"f()"<<endl;
			return compute_scalar_error(cb,ccerr);
		}
		void update_grad() {
		//	cout<<"update_grad()"<<endl;
		  compute_scalar_grad(cb,base_x,base_grad,eps_,ccerr);
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
			n_=init_scalar(cb,base_x)-base_x;
		}
		void reset() {
			reset_scalar(cb);
		}
		void eps_set(double d) {
			eps_=d;
		}
	};

}

#endif // VGTL_OPT_OPTIMIZE_SCALAR_CB_HPP
