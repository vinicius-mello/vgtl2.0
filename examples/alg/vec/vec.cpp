#include <iostream>
#include <vgtl/alg/vec.hpp>

using namespace std;
using namespace vgtl;

int main() {
	vec<3,float> a(_(1,1,0));
	vec<3,float> b(_(2,0,1));
	vec<3,float> c(_(0,0,1));
	cout<<"a="<<a<<endl;
	cout<<"b="<<b<<endl;
	cout<<"c="<<c<<endl;
	cout<<"3*a+b="<<a*3+b<<endl;
	cout<<"dot(a,b)="<<dot(a,b)<<endl;
	cout<<"cross(a,b)="<<cross(a,b)<<endl;
	cout<<"cross(a,b)="<<cross(_(a,b))<<endl;
	cout<<"det(a,b,c)="<<det(_(a,b,c))<<endl;
	vgtl::array<vec<3,float>,3> cf;
	cofat(_(a,b,c),cf);
	cout<<"cofat(a,b,c)="<<cf<<endl;
	vec<3> d;
	d=c;
	cout<<d<<endl;
	vec<4> e(_(2,1,1,1));
	vec<4> f(_(2,3,1,1));
	vec<4> g(_(2,1,6,1));
	vec<4> h=cross(_(e,f,g));
	cout<<"cross(e,f,g)="<<h<<endl;
	cout<<"dot(e,h)="<<dot(e,h)<<endl;
	cout<<"dot(f,h)="<<dot(f,h)<<endl;
	cout<<"dot(g,h)="<<dot(g,h)<<endl;

}
