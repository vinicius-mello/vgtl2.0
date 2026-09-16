#include <iostream>
#include <vgtl/alg/quat.hpp>

using namespace std;
using namespace vgtl;

int main() {
	quat<float> a(_(1,1,0,1));
	quat<float> b(_(2,0,1,0));
	quat<float> c(_(3,3,0,2));
	cout<<"a="<<a<<endl;
	cout<<"b="<<b<<endl;
	cout<<"c="<<c<<endl;
	cout<<"3*a+b*c="<<3.0f*a+b*c<<endl;
	cout<<"a/b="<<a/b<<endl;
}
