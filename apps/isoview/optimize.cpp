#include <iostream>
#include <cmath>
#include "app.hpp"
#include <vgtl/utl/pair_tie.hpp>
#include <vgtl/geo/iso.hpp>
#include <vgtl/geo/diffpars.hpp>
#include <vgtl/opt/optimize_diffpars_cb.hpp>
#include <vgtl/alg/bb.hpp>
#include <vgtl/alg/linalg.hpp>
#include <vgtl/alg/barygrad.hpp>
#include <vgtl/alg/bb_solve.hpp>

#define MAX_VAR 5000


using namespace std;
using namespace vgtl;

void eval(const bbform<4,double>& bb2,const bbform<4,double>& bb3,
          const bbform<4,double>& bb4, point<4,double>& bp,
          vgtl::array<vec<4,double>,2>& bv) {
  for(int k=0;k<app::kc3;++k)
    bp=jacobian2(bb4,bp,bv);
  for(int k=0;k<app::kc2;++k)
    bp=jacobian2(bb3,bp,bv);
  for(int k=0;k<app::kc1;++k)
    bp=jacobian2(bb2,bp,bv);
}

double comp_facet_err(const T& t, Facet(T) f) {
  vgtl::array<Vertex(T),3> vsf;
  vertices(t,f,vsf);
  int m,p;
  pair_tie(p,m)=iso::dim(t,vsf);
  if(p==-1 || m==-1) return 0;
  Cell(T) ce[2];
  pair_tie(ce[0],ce[1])=cells(t,f);
  if(ce[0]==ce[1]) return 0;
  //static const double sample[]={0,.333,.666,1};
  static const double sample[]={0,.25,.5,.75,1};
  //static const double sample[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  const index_t sz=sizeof(sample)/sizeof(double);
  bbform<4,double> bb2;
  bbform<4,double> bb3;
  bbform<4,double> bb4;
  vec<3,double> n[2][sz];
  vgtl::array<point<3,double>,2> bpsf;
  vgtl::array<point<4,double>,2> bpce;
  vgtl::array<point<3,double>,4> ps;
  vgtl::array<Vertex(T),4> vs;

  iso::basis(t,vsf,bpsf);
  double err=0;
  for(int i=0;i<2;++i) {
    diffpars(t,ce[i],bb2,Edge(T)());
    diffpars(t,ce[i],bb3,Facet(T)());
    diffpars(t,ce[i],bb4,Cell(T)());
    vertices(t,ce[i],vs);
    for(int j=0;j<=3;++j) ps[j]=euclidean_point(t,vs[j]);
    vgtl::array<point<4,double>,3> bps;
    iso::basis(t,ce[i],vs,bps);
    int ind=0;
    while(face_op(t,ce[i],ind)!=f) ++ind;
    unproject(ind,bpsf[0],bpce[0]);
    unproject(ind,bpsf[1],bpce[1]);
    for(index_t j=0;j<sz;++j) {
      point<4,double> bp=bpce[0]+sample[j]*(bpce[1]-bpce[0]);
      vgtl::array<vec<4,double>,2> bv;
      bv[0]=bps[1]-bps[0];
      bv[1]=bps[2]-bps[0];
      eval(bb2,bb3,bb4,bp,bv);
      vec<3,double> tg[2];
      tg[0]=barycentric_combination(bv[0],ps);
      tg[1]=barycentric_combination(bv[1],ps);
      n[i][j]=cross(tg[0],tg[1]);
      normalize(n[i][j]);
    }
  }
  for(index_t j=0;j<sz;++j) {
    double ne=dot(n[0][j],n[1][j]);
    //ne=1-ne*ne;
    //err+=ne;
    ne=1-ne;
    err+=ne*ne;
  }
  return err/sz;
}

vector<double> base_x(MAX_VAR);
vector<double> base_grad(MAX_VAR);
lbfgsb01 opt(MAX_VAR,20);

void smooth(complex_buffer<T>& cb, int it) {
  optimize_facet_error_cb<T> drv(cb,comp_facet_err,
                                 20,&base_x[0],&base_grad[0],1.0e-6);
  opt.pgtol_set(1.0e-7);
  opt.factr_set(1.0e4);
  opt.print_set(false);
  opt.max_iterations_set(it);
  opt.optimize(drv);
}

