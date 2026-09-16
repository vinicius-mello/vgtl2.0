#include <iostream>
#include <vector>
#include <vgtl/opt/simulated_annealing.hpp>
#include <vgtl/alg/point.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <gsl/gsl_rng.h>

using namespace std;
using namespace vgtl;

struct uniform_01 {
  gsl_rng * r;
  uniform_01() {
    r=gsl_rng_alloc(gsl_rng_mt19937);
  }
  // boltzmann_criteria/tsallis_criteria store their RandomGenerator by
  // value, so this type must be safely copyable: a plain shallow copy of
  // `r` would leave two owners of the same gsl_rng*, and whichever is
  // destroyed second double-frees it. gsl_rng_clone gives each copy its
  // own independent generator (same state, further draws diverge, which
  // is fine here since only one distribution ever consumes it going
  // forward).
  uniform_01(const uniform_01& other) : r(gsl_rng_clone(other.r)) {}
  uniform_01& operator=(const uniform_01& other) {
    if(this!=&other) {
      gsl_rng_free(r);
      r=gsl_rng_clone(other.r);
    }
    return *this;
  }
  ~uniform_01() {
    gsl_rng_free(r);
  }

  double operator()() {
    return gsl_rng_uniform(r);
  }
};

struct uniform_smallint {
  unsigned long n;
  gsl_rng * r;
  uniform_smallint(unsigned long _n) : n(_n) {
    r=gsl_rng_alloc(gsl_rng_mt19937);
  }
  // see uniform_01's copy constructor for why this is needed
  uniform_smallint(const uniform_smallint& other) : n(other.n), r(gsl_rng_clone(other.r)) {}
  uniform_smallint& operator=(const uniform_smallint& other) {
    if(this!=&other) {
      n=other.n;
      gsl_rng_free(r);
      r=gsl_rng_clone(other.r);
    }
    return *this;
  }
  ~uniform_smallint() {
    gsl_rng_free(r);
  }

  int operator()() {
    return gsl_rng_uniform_int(r,n);
  }
};

int figure=1;
vector<point<2,double> > sites(20);

void fill_sites() {
  uniform_01 r;
	for(int i=0;i<sites.size();++i) {
		point<2,double> p;
		p[0]=10*r();
		p[1]=10*r();
		sites[i]=p;
	}
}

void print_sites_initial() {
	cout<<"beginfig("<<figure<<")"<<endl;
	cout<<"\tnumeric u;"<<endl;
	cout<<"\tu=1cm;"<<endl;
	cout<<"\tpickup pencircle scaled 4bp;"<<endl;
	for(int i=0;i<sites.size();++i) {
		point<2,double> p=sites[i];
		cout<<"\tdraw ("<<p[0]<<"u,"<<p[1]<<"u);"<<endl;
	}
	cout<<"endfig;"<<endl;
	figure++;
}

void print_sites(const vector<int>& route) {
	cout<<"beginfig("<<figure<<")"<<endl;
	cout<<"\tnumeric u;"<<endl;
	cout<<"\tu=1cm;"<<endl;
	cout<<"\tdraw ";
	for(int i=0;i<sites.size();++i) {
		point<2,double> p=sites[route[i]];
		cout<<"("<<p[0]<<"u,"<<p[1]<<"u)..";
	}
	cout<<"cycle;"<<endl;
	cout<<"\tpickup pencircle scaled 4bp;"<<endl;
	for(int i=0;i<sites.size();++i) {
		point<2,double> p=sites[i];
		cout<<"\tdraw ("<<p[0]<<"u,"<<p[1]<<"u);"<<endl;
	}
	cout<<"endfig;"<<endl;
	figure++;
}

class tsp {
	int nsites;
  uniform_smallint r;
	vector<int> route;
	vector<int> best_route;
	vgtl::array<int,2> trans;
	double * distance_matrix;
	public:
	tsp(const vector<point<2,double> >& cities)
          : nsites(cities.size()), r(nsites),
	route(nsites), best_route(nsites) {
		for(int i=0;i<nsites;++i) route[i]=i;
		distance_matrix=new double[nsites*nsites];
		for(int i=0;i<nsites;++i) 
			for(int j=0;j<nsites;++j) {
				vec<2,double> d=cities[i]-cities[j];
				distance_matrix[i*nsites+j]=dot(d,d);
			}

	}
	~tsp() {
		delete[] distance_matrix;
	}
	double energy();
	double change_configuration();
	void save_configuration();
	void restore_configuration();
	void accept_configuration();
};

double tsp::energy()
{
	double E=0;
  for (int i = 0; i < nsites; ++i) {
    E += distance_matrix[route[i]*nsites+route[(i + 1) % nsites]];
  }
  return E;
}

double tsp::change_configuration()
{
  trans[0]=r();
  do {
    trans[1]=r();
  } while (trans[0]==trans[1]);
	double Ec=energy();
  swap(route[trans[0]],route[trans[1]]);
	Ec=energy()-Ec;
  swap(route[trans[0]],route[trans[1]]);
	return Ec;
}

void tsp::accept_configuration()
{
  swap(route[trans[0]],route[trans[1]]);
}

void tsp::save_configuration()
{
  for(int i=0;i<nsites;++i) {
		best_route[i]=route[i];
		cerr<<best_route[i]<<" ";
	}
	cerr<<endl;
}

void tsp::restore_configuration()
{
  for(int i=0;i<nsites;++i) route[i]=best_route[i];
	print_sites(best_route);
}

struct never_criteria {
	bool operator()(double ech, double T) { return false; }
};

struct allways_criteria {
	bool operator()(double ech, double T) { return true; }
};

int main(int argc, char * argv[])
{
  gsl_rng_env_setup();
	fill_sites();
	cout<<"prologues:=2;"<<endl;
	print_sites_initial();
	double K=1.0;
  uniform_01 r;
  boltzmann_criteria<uniform_01> bc(K,r);
  tsallis_criteria<uniform_01> tc(2,K,r);
	never_criteria nc;
	allways_criteria ac;
	geometric_cooling_schedule cs(5000,0.0001,0.95);
	simple_tries_handler tr(3000);
	tsp t(sites);
//	double min=simulated_annealing(t,cs,tr,tc);
  cerr<<"BP"<<endl;
	double min=simulated_annealing(t,cs,tr,bc);
//	double min=simulated_annealing(t,cs,tr,ac);
	cout<<"bye;"<<endl;
}


