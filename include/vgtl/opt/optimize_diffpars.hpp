#ifndef VGTL_OPT_OPTIMIZE_DIFFPARS_HPP
#define VGTL_OPT_OPTIMIZE_DIFFPARS_HPP

#include <iostream>
#include <list>
#include <vgtl/top/star.hpp>
#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/opt/lbfgsb01.hpp>


namespace vgtl {

	using std::list;
	using namespace std;

	template <class T, class CompCellErr>
	double compute_cell_error(T& t, CompCellErr ccerr) {
		double r=0;
		Cell_it(T) cei, ceend;
		for(simplices(t,cei,ceend); cei!=ceend; ++cei) {
			Cell(T) ce=*cei;
			if(!is_current(t,ce)) continue;
			double err=ccerr(t,ce);
			cell_error_set(t,ce,err);
			r+=err;
		}
		return r;
	}
	
	template <class T, class CompFacetErr>
	double compute_facet_error(T& t, CompFacetErr cferr) {
		double r=0;
		Facet_it(T) fi, fend;
		for(simplices(t,fi,fend); fi!=fend; ++fi) {
			Facet(T) f=*fi;
			if(!is_current(t,f)) continue;
			double err=cferr(t,f);
			facet_error_set(t,f,err);
			r+=err;
		}
		return r;
	}

	template <class T, dim_t k, class CompCellErr> 
	void compute_cell_grad(T& t, Simplex(T,k) s,
		double * grad, double eps, CompCellErr ccerr) {
		list<Cell(T)> st;
		star(t,s,back_inserter(st));
		double err=0;
		for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
			Cell(T) cv=*j;
			err+=cell_error(t,cv);
		}
		double * xp=diffpars_x(t,s);
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
				err_eps+=ccerr(t,cv);
			}
			grad[i]=(err_eps-err)/eps;
			grad[k]-=grad[i];
			for(int j=0;j<=k;++j) {
				if(j!=i) xp[j]+=epsk;
				else xp[j]-=epsk1;
			}
		}	
	}
					
	template <class T, dim_t k, class CompFacetErr> 
	void compute_facet_grad(T& t, Simplex(T,k) s,
		double * grad, double eps, CompFacetErr cferr) {
		const dim_t d=Dim(T);
		list<Cell(T)> st;
		star(t,s,back_inserter(st));
		set<Facet(T)> sf;
		for(typename list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j) {
			Cell(T) ce=*j;
			for(int i=0;i<=d;++i) {
				Facet(T) f=face_op(t,ce,i);
				sf.insert(f);
			}
		}
		double err=0;
	  for(typename set<Facet(T)>::iterator j=sf.begin();
				j!=sf.end(); ++j) {
			Facet(T) f=*j;
			err+=facet_error(t,f);
		}
		double * xp=diffpars_x(t,s);
		grad[k]=0;
		double epsk=eps/(k+1);
		double epsk1=k*eps/(k+1);
		for(int i=0;i<k;++i) {
			for(int j=0;j<=k;++j) {
				if(j!=i) xp[j]-=epsk;
				else xp[j]+=epsk1;
			}
			double err_eps=0;
	  	for(typename set<Facet(T)>::iterator j=sf.begin();
					j!=sf.end(); ++j) {
				Facet(T) f=*j;
				err_eps+=cferr(t,f);
			}
			grad[i]=(err_eps-err)/eps;
			grad[k]-=grad[i];
			for(int j=0;j<=k;++j) {
				if(j!=i) xp[j]+=epsk;
				else xp[j]-=epsk1;
			}
		}	
	}
	
	template <class T, dim_t k, class CompCellErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_cell_grad(T& t, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr);

	template <class T, dim_t k, class CompCellErr>
	typename disable_if<cell_c<T,k>,void>::type
	compute_cell_grad(T& t, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			compute_cell_grad(t,s,base_grad+(diffpars_x(t,s)-base_x),eps,ccerr);
		}
		compute_cell_grad<T,k+1>(t,base_x,base_grad,eps,ccerr);
	}
	
	template <class T, dim_t k, class CompCellErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_cell_grad(T& t, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			compute_cell_grad(t,s,base_grad+(diffpars_x(t,s)-base_x),eps,ccerr);
		}
	}
	
	template <class T, class CompCellErr>
	void compute_cell_grad(T& t, double * base_x,
		double * base_grad, double eps, CompCellErr ccerr) {
		compute_cell_grad<T,1>(t,base_x,base_grad,eps,ccerr);
	}


	
	template <class T, dim_t k, class CompFacetErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_facet_grad(T& t, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr);

	template <class T, dim_t k, class CompFacetErr>
	typename disable_if<cell_c<T,k>,void>::type
	compute_facet_grad(T& t, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			compute_facet_grad(t,s,base_grad+(diffpars_x(t,s)-base_x),eps,cferr);
		}
		compute_facet_grad<T,k+1>(t,base_x,base_grad,eps,cferr);
	}
	
	template <class T, dim_t k, class CompFacetErr>
	typename enable_if<cell_c<T,k>,void>::type
	compute_facet_grad(T& t, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			compute_facet_grad(t,s,base_grad+(diffpars_x(t,s)-base_x),eps,cferr);
		}
	}
	
	template <class T, class CompFacetErr>
	void compute_facet_grad(T& t, double * base_x,
		double * base_grad, double eps, CompFacetErr cferr) {
		compute_facet_grad<T,1>(t,base_x,base_grad,eps,cferr);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,double *>::type
	init_diffpars(T& t, double * base_x);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,double *>::type
	init_diffpars(T& t, double * base_x) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			diffpars_x_set(t,s,base_x);
			base_x+=k+1;
		}
		return init_diffpars<T,k+1>(t,base_x);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,double *>::type
	init_diffpars(T& t, double * base_x) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			diffpars_x_set(t,s,base_x);
			base_x+=k+1;
		}
		return base_x;
	}

	template <class T>
	double * init_diffpars(T& t, double * base_x) {
		return init_diffpars<T,1>(t,base_x);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	reset_diffpars(T& t);

	template <class T, dim_t k>
	typename disable_if<cell_c<T,k>,void>::type
	reset_diffpars(T& t) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			diffpars_x_set(t,s,0);
		}
		reset_diffpars<T,k+1>(t);
	}
	
	template <class T, dim_t k>
	typename enable_if<cell_c<T,k>,void>::type
	reset_diffpars(T& t) {
		Simplex_it(T,k) si, send;
		for(simplices(t,si,send);si!=send;++si) {
			Simplex(T,k) s=*si;
			if(!is_current(t,s)) continue;
			diffpars_x_set(t,s,0);
		}
	}

	template <class T>
	void reset_diffpars(T& t) {
		reset_diffpars<T,1>(t);
	}

	template <class T, class CompCellErr=double (*)(const T&, Cell(T))>
	class optimize_cell_error {
		double * base_x;
		double * base_grad;
		int n_;
		int m_;
		T& t;
		CompCellErr ccerr;
		double eps_;
		public:
		optimize_cell_error(T& _t, CompCellErr _ccerr, int _m, double * base_x_, 
			double * base_grad_, double _eps=0.000001) : m_(_m), 
			eps_(_eps), base_x(base_x_), base_grad(base_grad_), t(_t), ccerr(_ccerr) {
		}
		double f() {
		//	cout<<"f()"<<endl;
			return compute_cell_error(t,ccerr);
		}
		void update_grad() {
		//	cout<<"update_grad()"<<endl;
		  compute_cell_grad(t,base_x,base_grad,eps_,ccerr);
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
			n_=init_diffpars(t,base_x)-base_x;
		}
		void reset() {
			reset_diffpars(t);
		}
		void eps_set(double d) {
			eps_=d;
		}
	};

	template <class T, class CompFacetErr=double (*)(const T&, Facet(T))>
	class optimize_facet_error {
		double * base_x;
		double * base_grad;
		int n_;
		int m_;
		T& t;
		CompFacetErr cferr;
		double eps_;
		public:
		optimize_facet_error(T& _t, CompFacetErr _cferr, int _m, double * base_x_, 
			double * base_grad_, double _eps=0.000001) : m_(_m), 
			eps_(_eps), base_x(base_x_), base_grad(base_grad_), t(_t), cferr(_cferr) {
		}
		double f() {
		//	cout<<"f()"<<endl;
			return compute_facet_error(t,cferr);
		}
		void update_grad() {
		//	cout<<"update_grad()"<<endl;
		  compute_facet_grad(t,base_x,base_grad,eps_,cferr);
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
			n_=init_diffpars(t,base_x)-base_x;
		}
		void reset() {
			reset_diffpars(t);
		}
		void eps_set(double d) {
			eps_=d;
		}
	};

}

#endif // VGTL_OPT_OPTIMIZE_DIFFPARS_HPP
