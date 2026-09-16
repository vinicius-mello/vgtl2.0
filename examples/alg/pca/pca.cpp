#include <iostream>
#include <list>
#include <vgtl/alg/pca.hpp>

using namespace std;
using namespace vgtl;

int main() {
	vgtl::array<vec<2,double>,2> axes;
	vec<2,double> sizes;
	list<point<2,double> > lp;
	point<2,double> bc;
	lp.push_back(_<double>(-1,0));
	lp.push_back(_<double>(1,0));
	lp.push_back(_<double>(0,-0.5));
	lp.push_back(_<double>(0,0.5));
	pca(lp.begin(),lp.end(),bc,axes,sizes);
	cout<<bc<<endl;
	cout<<axes[0]<<endl;
	cout<<axes[1]<<endl;
	cout<<sizes[0]<<endl;
	cout<<sizes[1]<<endl;
}

