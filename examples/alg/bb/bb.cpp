#include <iostream>
#include <vgtl/alg/bb_solve.hpp>
#include <vgtl/alg/bb.hpp>
#include <vgtl/comb/triangular.hpp>
#include <vgtl/utl/array_cons.hpp>

using namespace std;
using namespace vgtl;

int main() {
	bbform<3,double> bb(3);
	bb[_(3,0,0)]=_(1.0,0.0,0.0);
	bb[_(0,3,0)]=_(0.0,1.0,0.0);
	bb[_(0,0,3)]=_(0.0,0.0,1.0);
	bb[_(2,1,0)]=_(0.7,0.3,0.0);
	bb[_(2,0,1)]=_(0.8,0.0,0.2);
	bb[_(0,2,1)]=_(0.0,0.3,0.7);
	bb[_(0,1,2)]=_(0.0,0.1,0.9);
	bb[_(1,0,2)]=_(0.3,0.0,0.7);
	bb[_(1,2,0)]=_(0.4,0.6,0.0);
	bb[_(1,1,1)]=_(0.1,0.1,0.8);
	vec<3,double> ones= _(1.0,1.0,1.0);
	{
		int i=0;
		triangular_traverser<3,double> tt(90);
		do {
			vec<3,double> v=*tt;
			vec<3,double> p=1.0/93.0*(ones+v);
			eval(bb,reinterpret_cast<const point<3,double>& >(p));
			++i;
		} while(++tt);
		cout<<"Iterations: "<<i<<endl;
	}
	
	{
		int i=0;
		vgtl::array<vec<3,double>,3> J;
		triangular_traverser<3,double> tt(90);
		do {
			vec<3,double> v=*tt;
			vec<3,double> p=1.0/93.0*(ones+v);
			solveJ(bb,reinterpret_cast<const point<3,double>& >(p),J);
			++i;
		} while(++tt);
		cout<<"Iterations: "<<i<<endl;
	}
}

