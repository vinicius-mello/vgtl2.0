#include <iostream>
#include <vgtl/alg/point.hpp>

using namespace std;
using namespace vgtl;

int main() {
	point<2,float> p;
	point<3,float> q;
	point<2,float> a(_(1,1));
	point<2,float> b(_(2,1));
	point<2,float> c(_(2,2));
	cout<<"a="<<a<<endl;
	cout<<"b="<<b<<endl;
	cout<<"c="<<c<<endl;
	cout<<"a+3*(b-c)="<<a+3.0f*(b-c)<<endl;
	cout<<"bracket(a,b,c)="<<bracket(_(a,b,c))<<endl;
	cout<<"p=barycenter(a,b,c)="<<(p=barycenter(_(a,b,c)))<<endl;
	cout<<"q=barycentric_coordinates(p|a,b,c)="<<
					(q=barycentric_coordinates(p,_(a,b,c)))<<endl;
	cout<<"barycentric_combination(q|a,b,c)="<<
					barycentric_combination(q,_(a,b,c))<<endl;
	point<2,double> r;
	r=c;
	cout<<r<<endl;
}
