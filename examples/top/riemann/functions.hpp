#ifndef VGTL_RIEMANN_FUNCTIONS_HPP
#define VGTL_RIEMANN_FUNCTIONS_HPP

// Catalog of test curves F(w,z) = w^n + f_{n-1}(z) w^{n-1} + ... + f_0(z)
// (monic in w) for the Riemann-surface triangulation example (riemann.cpp),
// selected from the command line by index (--function N; --list-functions
// prints this list). Add a new curve by appending an entry to
// function_catalog() below -- nothing in riemann.cpp is specific to a
// particular degree or curve.

#include <complex>
#include <vector>
#include <string>
#include <iostream>

typedef std::complex<double> cx;

// Each f_i(z) is given as a list of complex coefficients, lowest degree
// first -- the same input shape the method was specified with.
struct poly_F {
	int n;
	std::vector<std::vector<cx> > c; // c[i] = coefficients of f_i(z), i=0..n-1

	cx f(int i, cx z) const {
		const std::vector<cx>& p=c[i];
		cx r(0,0);
		for(int j=(int)p.size()-1;j>=0;--j) r=r*z+p[j];
		return r;
	}
	cx fprime(int i, cx z) const {
		const std::vector<cx>& p=c[i];
		cx r(0,0);
		for(int j=(int)p.size()-1;j>=1;--j) r=r*z+double(j)*p[j];
		return r;
	}
};

struct catalog_entry {
	std::string name;
	std::string description;
	poly_F F;
};

// Index 0 is the default (used when --function is not given).
inline std::vector<catalog_entry>& function_catalog() {
	static std::vector<catalog_entry> cat;
	if(!cat.empty()) return cat;

	{
		catalog_entry e;
		e.name="elliptic";
		e.description="w^2 - z^3 + z : smooth elliptic curve (genus 1); "
			"branch points at z=0,+1,-1.";
		e.F.n=2;
		e.F.c.resize(2);
		e.F.c[1].push_back(cx(0,0));                  // f_1(z) = 0
		cx f0[4]={cx(0,0),cx(1,0),cx(0,0),cx(-1,0)};   // f_0(z) = z - z^3
		e.F.c[0].assign(f0,f0+4);
		cat.push_back(e);
	}
	{
		catalog_entry e;
		e.name="parabola";
		e.description="w^2 - z : smooth parabola; one branch point at "
			"(0,0), which lands exactly on the octahedron seed vertex z=0.";
		e.F.n=2;
		e.F.c.resize(2);
		e.F.c[1].push_back(cx(0,0));                   // f_1(z) = 0
		cx f0[2]={cx(0,0),cx(-1,0)};                    // f_0(z) = -z
		e.F.c[0].assign(f0,f0+2);
		cat.push_back(e);
	}

	return cat;
}

inline void print_function_catalog(std::ostream& out) {
	std::vector<catalog_entry>& cat=function_catalog();
	for(size_t i=0;i<cat.size();++i)
		out<<"  "<<i<<": "<<cat[i].name<<" -- "<<cat[i].description<<"\n";
}

#endif // VGTL_RIEMANN_FUNCTIONS_HPP
