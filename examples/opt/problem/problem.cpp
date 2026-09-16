#include <iostream>
#include <vgtl/opt/lbfgsb01.hpp>

using namespace std;

class driver {
	double x_[3];
	double grad_[3];
	public:
	double f() {
		double s=x_[0]+x_[1]+x_[2];
		double u=x_[0];//s;
		double v=x_[1];//s;
		return (u-.1)*(u-.1)+(v-.2)*(v-.2);
	}
	double update_grad() {
		double eps=0.00000001;
		double t=f();
		x_[0]+=2*eps/3;x_[1]-=eps/3;x_[2]-=eps/3;
		grad_[0]=(f()-t)/eps;
		x_[0]-=2*eps/3;x_[1]+=eps/3;x_[2]+=eps/3;
		x_[0]-=eps/3;x_[1]+=2*eps/3;x_[2]-=eps/3;
		grad_[1]=(f()-t)/eps;
		x_[0]+=eps/3;x_[1]-=2*eps/3;x_[2]+=eps/3;
		//x_[0]-=eps/2;x_[1]-=eps/2;x_[2]+=eps;
		//grad_[2]=(f()-t)/eps;
		//x_[0]+=eps/2;x_[1]+=eps/2;x_[2]-=eps;
		grad_[2]=-grad_[0]-grad_[1];
    //exit(1);
		return t;
	}
	int n() {
		return 3;
	}
	int m() {
		return 2;
	}
	double * x() {
		return (double *)x_;
	}
	double * grad() {
		return (double *)grad_;
	}
	void init() {
		x_[0]=x_[1]=x_[2]=1.0/3;
	}
	void reset() {
	}
};

int main() {
	vgtl::lbfgsb01 op(10,3);
	op.factr_set(1.0e+4);
	op.pgtol_set(1.0e-9);
	driver drv;
	op.optimize(drv);
	double * x=drv.x();
	//double s=x[0]+x[1]+x[2];
	//cout<<(x[0]/s)<<" "<<(x[1]/s)<<" "<<(x[2]/s)<<endl;
	cout<<(x[0])<<" "<<(x[1])<<" "<<(x[2])<<endl;
}
