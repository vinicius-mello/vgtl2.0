#include <iostream>
#include <vgtl/alg/linalg.hpp>
#include <vgtl/alg/barygrad.hpp>
#include <vgtl/alg/barychange.hpp>

using namespace std;
using namespace vgtl;

int main() {
	vgtl::array<vec<3,double>,3> a;
  vgtl::array<dim_t,3> p;
	a[0]=_<double>(2,0,0);
	a[1]=_<double>(0,2,0);
	a[2]=_<double>(0,0,2);
	LU_decomp(a,p);
	vec<3,double> x;
	vec<3,double> b=_<double>(1,1,4);
	LU_solve(a,p,b,x);
	cout<<x<<endl;
	vgtl::array<vec<3,double>,3> ai;
	LU_invert(a,p,ai);
	cout<<ai<<endl;
	vgtl::array<point<2,double>,3> ps;
	ps[0]=_<double>(0,0);
	ps[1]=_<double>(2,0);
	ps[2]=_<double>(0,1);
	vec<3,double> f=_<double>(+1,-1,+1);

	vec<2,double> bg;
	bg=barygrad(ps,f);
	cout<<dot(bg,bg)<<endl;
	cout<<bg<<endl;
	
	//cout<<barycentric_combination(bg,ps)<<endl;
	a[0]=_<double>(2,1,-1);
	a[1]=_<double>(1,2,5);
	a[2]=_<double>(-1,5,2);
	vgtl::array<vec<3,double>,3> eigvec;
	vec<3,double> eigval;
	eigenvalues(a,eigvec,eigval);
	cout<<eigval<<endl;
}

