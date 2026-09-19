#ifndef IMPLICIT3D_FUNCTIONS_HPP
#define IMPLICIT3D_FUNCTIONS_HPP

// Catalog of test implicit surfaces f(x,y,z)=0 for the plain R^3 case:
// unlike examples/top/riemann(_PC2)'s codimension-2 curves (a complex
// F=0 in C^2/CP^2, real ambient dim 4), this is the classical
// codimension-1 case -- one REAL equation, one REAL ambient dimension
// fewer, and no distinguished chart/pole handling needed at all. Small
// catalog, so gradients are just hand-written analytically (no need
// for functions_pc2.hpp's symbolic-differentiation machinery).

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

typedef double (*scalar_fun3)(double,double,double);

struct catalog_entry_implicit3d {
	std::string name, description;
	scalar_fun3 f, fx, fy, fz;
	double Lx, Ly, Lz; // half-extents of the axis-aligned box seed domain
};

// --- sphere: constant curvature, the trivial baseline. ---
inline double sphere_f (double x,double y,double z){ return x*x+y*y+z*z-1.0; }
inline double sphere_fx(double x,double y,double z){ return 2.0*x; }
inline double sphere_fy(double x,double y,double z){ return 2.0*y; }
inline double sphere_fz(double x,double y,double z){ return 2.0*z; }

// --- ellipsoid: anisotropic curvature (much sharper near the short
// z-axis than around the long x-equator) -- a direct test of whether
// the curvature-adaptive criterion actually concentrates refinement
// where the surface bends more, rather than refining uniformly. ---
inline double ellip_f (double x,double y,double z){
	double a=2.0,b=1.0,c=0.6;
	return (x/a)*(x/a)+(y/b)*(y/b)+(z/c)*(z/c)-1.0;
}
inline double ellip_fx(double x,double y,double z){ double a=2.0; return 2.0*x/(a*a); }
inline double ellip_fy(double x,double y,double z){ double b=1.0; return 2.0*y/(b*b); }
inline double ellip_fz(double x,double y,double z){ double c=0.6; return 2.0*z/(c*c); }

// --- torus: genus 1 (a real hole through the domain), tube curvature
// varies smoothly with position around the ring (tighter on the inner
// equator, gentler on the outer one) -- also a nontrivial topology
// test, unlike the two genus-0 surfaces above. Singular only on the
// z-axis itself (s=0), which the torus (R=1) never comes near. ---
inline double torus_f (double x,double y,double z){
	double R=1.0, r=0.35;
	double s=std::sqrt(x*x+y*y);
	double d=s-R;
	return d*d+z*z-r*r;
}
inline double torus_fx(double x,double y,double z){
	double R=1.0;
	double s=std::sqrt(x*x+y*y);
	if(s<1e-12) return 0.0;
	return 2.0*(s-R)*x/s;
}
inline double torus_fy(double x,double y,double z){
	double R=1.0;
	double s=std::sqrt(x*x+y*y);
	if(s<1e-12) return 0.0;
	return 2.0*(s-R)*y/s;
}
inline double torus_fz(double x,double y,double z){ return 2.0*z; }

inline std::vector<catalog_entry_implicit3d>& function_catalog_implicit3d() {
	static std::vector<catalog_entry_implicit3d> cat;
	if(cat.empty()) {
		{
			catalog_entry_implicit3d e;
			e.name="sphere";
			e.description="x^2+y^2+z^2-1=0 : unit sphere, constant curvature (baseline)";
			e.f=sphere_f; e.fx=sphere_fx; e.fy=sphere_fy; e.fz=sphere_fz;
			e.Lx=1.5; e.Ly=1.5; e.Lz=1.5;
			cat.push_back(e);
		}
		{
			catalog_entry_implicit3d e;
			e.name="ellipsoid";
			e.description="(x/2)^2+y^2+(z/0.6)^2-1=0 : anisotropic curvature -- "
				"refinement should concentrate near the sharper z-poles, not the flatter x-equator";
			e.f=ellip_f; e.fx=ellip_fx; e.fy=ellip_fy; e.fz=ellip_fz;
			e.Lx=2.5; e.Ly=1.5; e.Lz=1.0;
			cat.push_back(e);
		}
		{
			catalog_entry_implicit3d e;
			e.name="torus";
			e.description="(sqrt(x^2+y^2)-1)^2+z^2-0.35^2=0 : genus 1, curvature varies "
				"around the tube (inner vs outer equator)";
			e.f=torus_f; e.fx=torus_fx; e.fy=torus_fy; e.fz=torus_fz;
			e.Lx=1.5; e.Ly=1.5; e.Lz=0.5;
			cat.push_back(e);
		}
	}
	return cat;
}

inline void print_function_catalog_implicit3d(std::ostream& out) {
	std::vector<catalog_entry_implicit3d>& cat=function_catalog_implicit3d();
	for(size_t i=0;i<cat.size();++i)
		out<<"  "<<i<<": "<<cat[i].name<<" -- "<<cat[i].description<<std::endl;
}

#endif // IMPLICIT3D_FUNCTIONS_HPP
