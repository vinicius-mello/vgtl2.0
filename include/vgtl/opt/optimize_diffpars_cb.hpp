#ifndef VGTL_OPT_OPTIMIZE_DIFFPARS_CB_HPP
#define VGTL_OPT_OPTIMIZE_DIFFPARS_CB_HPP

#include <iostream>
#include <list>
#include <vgtl/top/star.hpp>
#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/opt/lbfgsb01.hpp>


namespace vgtl {

	using std::list;
	using namespace std;

	template <class T, class CompCellErr>
	double compute_cell_error(complex_buffer<T>& cb, CompCellErr ccerr) {
		double r=0;
		int c=0;
		for(typename set<Cell(T)>::iterator cei=cb.template begin<Dim(T)>();
				cei!=cb.template end<Dim(T)>(); ++cei) {
			Cell(T) ce=*cei;
			double err=ccerr(cb.t,ce);
			cell_error_set(cb.t,ce,err);
			r+=err;
			c++;
		}
		return r;
		//return r/c;
	}
	
	template <class T, class CompFacetErr>
	double compute_facet_error(complex_buffer<T>& cb, CompFacetErr cferr) {
		double r=0;
		int c=0;
		for(typename set<Facet(T)>::iterator fi=cb.template begin<Dim(T)-1>();
				fi!=cb.template end<Dim(T)-1>(); ++fi) {
			Facet(T) f=*fi;
			double err=cferr(cb.t,f);
			facet_error_set(cb.t,f,err);
			r+=err;
			c++;
		}
		return r;
		//return r/c;
	}

	template <class T, dim_t k, class CompCellErr> 
	void compute_cell_grad(complex_buffer<T>& cb, Simplex(T,k) s,
		double * grad, double eps, CompCellErr ccerr) {
		list<Cell(T)> st;
		star(cb.t,s,back_inserter(st));
		double err=0;
		int c=0;
		for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
			Cell(T) cv=*j;
			if(exists(cb,cv)) {
				err+=cell_error(cb.t,cv);
				c++;
			}
		}
		if(c==0) {
			for(int i=0;i<=k;++i) grad[i]=0;
			return;
		}
		double * xp=diffpars_x(cb.t,s);
		grad[k]=0;
		double epsk=eps/(k+1);
		double epsk1=k*eps/(k+1);
		for(int i=0;i<k;++i) {
			for(int j=0;j<=k;++j) {
				if(j!=i) xp[j]-=epsk;
				else xp[j]+=epsk1;
			}
			double err_eps=0;
			for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
				Cell(T) cv=*j;
				if(exists(cb,cv))
					err_eps+=ccerr(cb.t,cv);
			}
			grad[i]=(err_eps-err)/eps;
			//grad[i]=(err_eps-err)/c/eps;
			grad[k]-=grad[i];
			for(int j=0;j<=k;++j) {
				if(j!=i) xp[j]+=epsk;
				else xp[j]-=epsk1;
			}
		}	
	}
					
	template <class T, dim_t k, class CompFacetErr> 
	void compute_facet_grad(complex_buffer<T>& cb, Simplex(T,k) s,
		double * grad, double eps, CompFacetErr cferr) {
		const dim_t d=Dim(T);
		list<Cell(T)> st;
		star(cb.t,s,back_inserter(st));
		set<Facet(T)> sf;
		for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
			Cell(T) ce=*j;
		//	if(!exists(cb,ce)) continue;
      for(dim_t i=0;i<=d;++i) {
				Facet(T) f=face_op(cb.t,ce,i);
				if(exists(cb,f))
					sf.insert(f);
			}
		}
		double err=0;
		int c=0;
	  for(typename set<Facet(T)>::iterator j=sf.begin();
				j!=sf.end(); ++j) {
			Facet(T) f=*j;
			err+=facet_error(cb.t,f);
			c++;
		}
		if(c==0) {
      for(dim_t i=0;i<=k;++i) grad[i]=0;
			return;
		}
		double * xp=diffpars_x(cb.t,s);
		grad[k]=0;
		double epsk=eps/(k+1);
		double epsk1=k*eps/(k+1);
    for(dim_t i=0;i<k;++i) {
      for(dim_t j=0;j<=k;++j) {
				if(j!=i) xp[j]-=epsk;
				else xp[j]+=epsk1;
			}
			double err_eps=0;
	  	for(typename set<Facet(T)>::iterator j=sf.begin();
					j!=sf.end(); ++j) {
				Facet(T) f=*j;
				err_eps+=cferr(cb.t,f);
			}
			grad[i]=(err_eps-err)/eps;
			//grad[i]=(err_eps-err)/c/eps;
			grad[k]-=grad[i];
      for(dim_t j=0;j<=k;++j) {
				if(j!=i) xp[j]+=epsk;
				else xp[j]-=epsk1;
			}
		}	
	}
	
	template <class T, dim_t k, class CompCellErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_cell_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr);

	template <class T, dim_t k, class CompCellErr>
	typename disable_if<cell_c<T,k>,void>::type
	compute_cell_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			compute_cell_grad(cb,s,base_grad+(diffpars_x(cb.t,s)-base_x),eps,ccerr);
		}
		compute_cell_grad<T,k+1>(cb,base_x,base_grad,eps,ccerr);
	}
	
	template <class T, dim_t k, class CompCellErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_cell_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			compute_cell_grad(cb,s,base_grad+(diffpars_x(cb.t,s)-base_x),eps,ccerr);
		}
	}
	
	template <class T, class CompCellErr>
	void compute_cell_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		compute_cell_grad<T,1>(cb,base_x,base_grad,eps,ccerr);
	}
	
	template <class T, dim_t k, class CompFacetErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_facet_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr);

	template <class T, dim_t k, class CompFacetErr>
	typename disable_if<cell_c<T,k>,void>::type
	compute_facet_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			compute_facet_grad(cb,s,base_grad+(diffpars_x(cb.t,s)-base_x),eps,cferr);
		}
		compute_facet_grad<T,k+1>(cb,base_x,base_grad,eps,cferr);
	}
	
	template <class T, dim_t k, class CompFacetErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_facet_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			compute_facet_grad(cb,s,base_grad+(diffpars_x(cb.t,s)-base_x),eps,cferr);
		}
	}
	
	template <class T, class CompFacetErr>
	void compute_facet_grad(complex_buffer<T>& cb, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr) {
		compute_facet_grad<T,1>(cb,base_x,base_grad,eps,cferr);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,double *>::type
	init_diffpars(complex_buffer<T>& cb, double * base_x);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,double *>::type
	init_diffpars(complex_buffer<T>& cb, double * base_x) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			diffpars_x_set(cb.t,s,base_x);
			base_x+=k+1;
		}
		return init_diffpars<T,k+1>(cb,base_x);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,double *>::type
	init_diffpars(complex_buffer<T>& cb, double * base_x) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			diffpars_x_set(cb.t,s,base_x);
			base_x+=k+1;
		}
		return base_x;
	}

	template <class T>
	double * init_diffpars(complex_buffer<T>& cb, double * base_x) {
		return init_diffpars<T,1>(cb,base_x);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	reset_diffpars(complex_buffer<T>& cb);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	reset_diffpars(complex_buffer<T>& cb) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			diffpars_x_set(cb.t,s,0);
		}
		reset_diffpars<T,k+1>(cb);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	reset_diffpars(complex_buffer<T>& cb) {
		for(typename set<Simplex(T,k)>::iterator si=cb.template begin<k>();
				si!=cb.template end<k>();++si) {
			Simplex(T,k) s=*si;
			diffpars_x_set(cb.t,s,0);
		}
	}

	template <class T>
	void reset_diffpars(complex_buffer<T>& cb) {
		reset_diffpars<T,1>(cb);
	}

	template <class T, class CompCellErr=double (*)(const T&, Cell(T))>
	class optimize_cell_error_cb {
		double * base_x;
		double * base_grad;
		int n_;
		int m_;
		complex_buffer<T>& cb;
		CompCellErr ccerr;
		double eps_;
		public:
		optimize_cell_error_cb(complex_buffer<T>& _cb,
								CompCellErr _ccerr, int _m, double * base_x_, 
			double * base_grad_, double _eps=0.000001) : m_(_m), 
			eps_(_eps), base_x(base_x_), base_grad(base_grad_), cb(_cb),
			ccerr(_ccerr) {
		}
		double f() {
		//	cout<<"f()"<<endl;
			return compute_cell_error(cb,ccerr);
		}
		void update_grad() {
		//	cout<<"update_grad()"<<endl;
		  compute_cell_grad(cb,base_x,base_grad,eps_,ccerr);
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
			n_=init_diffpars(cb,base_x)-base_x;
		}
		void reset() {
			reset_diffpars(cb);
		}
		void eps_set(double d) {
			eps_=d;
		}
	};

	template <class T, class CompFacetErr=double (*)(const T&, Facet(T))>
	class optimize_facet_error_cb {
		int n_;
		int m_;
    double eps_;
    double * base_x;
    double * base_grad;
		complex_buffer<T>& cb;
		CompFacetErr cferr;
		public:
		optimize_facet_error_cb(complex_buffer<T>& _cb,
										CompFacetErr _cferr, int _m, double * base_x_, 
			double * base_grad_, double _eps=0.000001) : m_(_m), 
			eps_(_eps), base_x(base_x_), base_grad(base_grad_), cb(_cb),
			cferr(_cferr) {
		}
		double f() {
		//	cout<<"f()"<<endl;
			return compute_facet_error(cb,cferr);
		}
		void update_grad() {
		//	cout<<"update_grad()"<<endl;
		  compute_facet_grad(cb,base_x,base_grad,eps_,cferr);
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
			n_=init_diffpars(cb,base_x)-base_x;
		}
		void reset() {
			reset_diffpars(cb);
		}
		void eps_set(double d) {
			eps_=d;
		}
	};

}

#endif // VGTL_OPT_OPTIMIZE_DIFFPARS_CB_HPP
