#ifndef RIEMANN_PC2_FUNCTIONS_HPP
#define RIEMANN_PC2_FUNCTIONS_HPP

// Catalog of test curves for the CP^2 approach: each is a single
// homogeneous polynomial F(X,Y,Z) of degree d, given as a sparse list
// of monomials (i,j,k,coeff) with i+j+k=d. Unlike examples/top/riemann/
// functions.hpp (which needed per-curve f(i,z)/fprime(i,z) callbacks
// because F was only monic-in-w, not fully homogeneous), everything
// here -- F itself and its gradient -- is derived generically from the
// monomial list by symbolic differentiation (diff_poly3), so a new
// curve is just a new monomial list, nothing else to write.

#include <vector>
#include <complex>
#include <string>
#include <iostream>

typedef std::complex<double> cx;

struct term3 { int e[3]; cx c; }; // X^e0 Y^e1 Z^e2, coefficient c

struct poly_F3 {
	int d; // degree (0 for the zero polynomial)
	std::vector<term3> t;
};

inline cx ipow(cx base, int e) {
	cx r(1,0);
	for(int k=0;k<e;++k) r*=base;
	return r;
}

inline cx eval_poly3(const poly_F3& p, cx X, cx Y, cx Z) {
	cx r(0,0);
	for(size_t k=0;k<p.t.size();++k)
		r += p.t[k].c * ipow(X,p.t[k].e[0]) * ipow(Y,p.t[k].e[1]) * ipow(Z,p.t[k].e[2]);
	return r;
}

// d/dvar (var=0,1,2 for X,Y,Z), by the power rule term-by-term.
inline poly_F3 diff_poly3(const poly_F3& p, int var) {
	poly_F3 q;
	q.d = (p.d>0) ? p.d-1 : 0;
	for(size_t k=0;k<p.t.size();++k) {
		int e=p.t[k].e[var];
		if(e==0) continue;
		term3 nt=p.t[k];
		nt.c *= double(e);
		nt.e[var]=e-1;
		q.t.push_back(nt);
	}
	return q;
}

struct catalog_entry_pc2 {
	std::string name;
	std::string description;
	poly_F3 F;
};

inline poly_F3 mono(int i,int j,int k,cx c) {
	poly_F3 p; p.d=i+j+k; term3 t; t.e[0]=i; t.e[1]=j; t.e[2]=k; t.c=c; p.t.push_back(t);
	return p;
}
inline poly_F3 operator+(const poly_F3& a, const poly_F3& b) {
	poly_F3 r=a;
	r.t.insert(r.t.end(), b.t.begin(), b.t.end());
	r.d = a.d>b.d ? a.d : b.d;
	return r;
}

inline std::vector<catalog_entry_pc2>& function_catalog_pc2() {
	static std::vector<catalog_entry_pc2> cat;
	if(cat.empty()) {
		{
			// Smooth conic XY - Z^2 = 0 -- genus 0, constant curvature.
			// Doesn't pass through any of the 9 Hesse seed vertices
			// (checked directly), so no seed-degeneracy concerns.
			catalog_entry_pc2 e;
			e.name="conic";
			e.description="XY-Z^2=0 : smooth conic (genus 0), constant curvature";
			e.F = mono(1,1,0,cx(1,0)) + mono(0,0,2,cx(-1,0));
			cat.push_back(e);
		}
		{
			// Fermat cubic X^3+Y^3+Z^3=0 -- genus 1. Every member of the
			// Hesse pencil (this one included) passes through the same 9
			// base points used as the seed mesh's own vertices (see
			// hesse_points() in riemann_pc2.cpp) -- a real extraction
			// stress test: those 9 vertices sit exactly ON the curve
			// forever (Maubach never touches seed vertices), the direct
			// CP^2 analogue of examples/top/riemann/'s "branch point
			// pinned to a seed vertex" issue.
			catalog_entry_pc2 e;
			e.name="fermat_cubic";
			e.description="X^3+Y^3+Z^3=0 : smooth elliptic curve (genus 1); passes exactly through "
				"the seed mesh's 9 Hesse-point vertices (stress test, not a clean baseline)";
			e.F = mono(3,0,0,cx(1,0)) + mono(0,3,0,cx(1,0)) + mono(0,0,3,cx(1,0));
			cat.push_back(e);
		}
		{
			// Another Hesse-pencil member, less symmetric than Fermat but
			// with the SAME seed-vertex degeneracy (every pencil member
			// shares the 9 base points) -- useful for checking a result
			// isn't an artifact of the Fermat cubic's extra symmetry, NOT
			// for avoiding the seed-degeneracy issue itself.
			catalog_entry_pc2 e;
			e.name="cubic_generic";
			e.description="X^3+Y^3+Z^3-3XYZ=0 : Hesse pencil member away from Fermat, still passes "
				"through all 9 seed vertices, no extra symmetry with the ambient mesh beyond that";
			e.F = mono(3,0,0,cx(1,0)) + mono(0,3,0,cx(1,0)) + mono(0,0,3,cx(1,0)) + mono(1,1,1,cx(-3,0));
			cat.push_back(e);
		}
		{
			// Smooth quartic, genus 3 -- a more complex test curve, and
			// (checked directly) does NOT pass through any of the 9 seed
			// vertices, so a clean baseline unlike the two cubics above.
			catalog_entry_pc2 e;
			e.name="quartic";
			e.description="X^4+Y^4+Z^4=0 : smooth quartic (genus 3), no special relationship to "
				"the seed mesh's 9 vertices";
			e.F = mono(4,0,0,cx(1,0)) + mono(0,4,0,cx(1,0)) + mono(0,0,4,cx(1,0));
			cat.push_back(e);
		}
		{
			// The same "parabola" curve as examples/top/riemann/'s catalog,
			// F(w,z)=w^2-z, homogenized into a single CP^2 via w=X/Z, z=Y/Z
			// (Z the denominator/homogenizing coordinate -- matches
			// triangle_intersection's own chart-2 convention, though which
			// chart a given triangle actually uses at runtime is picked
			// per-triangle regardless): (X/Z)^2-(Y/Z)=0, times Z^2.
			catalog_entry_pc2 e;
			e.name="parabola";
			e.description="X^2-YZ=0 : homogenization of examples/top/riemann's w^2-z (w=X/Z,z=Y/Z); "
				"one branch point, at (w,z)=(0,0). Doesn't pass through any of the 9 seed vertices "
				"(checked directly), clean baseline";
			e.F = mono(2,0,0,cx(1,0)) + mono(0,1,1,cx(-1,0));
			cat.push_back(e);
		}
		{
			// The same "elliptic" curve as examples/top/riemann/'s catalog,
			// F(w,z)=w^2-z^3+z, homogenized the same way: (X/Z)^2-(Y/Z)^3+
			// (Y/Z)=0, times Z^3.
			catalog_entry_pc2 e;
			e.name="elliptic";
			e.description="X^2*Z-Y^3+Y*Z^2=0 : homogenization of examples/top/riemann's w^2-z^3+z "
				"(w=X/Z,z=Y/Z); smooth elliptic curve (genus 1), branch points at z=0,+1,-1. Passes "
				"exactly through ONE of the 9 seed vertices ([0:1:-1], checked directly) -- a much "
				"milder version of fermat_cubic's 9-point degeneracy, but not entirely clean either";
			e.F = mono(2,0,1,cx(1,0)) + mono(0,3,0,cx(-1,0)) + mono(0,1,2,cx(1,0));
			cat.push_back(e);
		}
	}
	return cat;
}

inline void print_function_catalog_pc2(std::ostream& out) {
	std::vector<catalog_entry_pc2>& cat=function_catalog_pc2();
	for(size_t i=0;i<cat.size();++i)
		out<<"  "<<i<<": "<<cat[i].name<<" -- "<<cat[i].description<<std::endl;
}

#endif // RIEMANN_PC2_FUNCTIONS_HPP
