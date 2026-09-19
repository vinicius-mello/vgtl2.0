#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <queue>
#include <set>
#include <map>
#include <complex>
#include <cmath>
#include <random>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/euler.hpp>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/alg/vec.hpp>
#include "functions.hpp"

// Seed mesh for C_infty^2, built as the simplicial product of two
// copies of the standard octahedral triangulation of the Riemann
// sphere C_infty (via stereographic projection), following the
// Eilenberg-Zilber "shuffle" construction: for each pair of triangles
// (one per factor) and each monotone path across the 3x3 grid of
// their vertices, we get one 4-simplex.

#define DIM 4

using namespace std;
using namespace vgtl;

typedef vgtl::nmt<DIM> T;

T t;

namespace vgtl {

	template <>
	struct extra_data<0> {
		vec<3,double> w_sphere;
		vec<3,double> z_sphere;
		int w_label; // -1 for vertices created by refinement
		int z_label;
		std::complex<double> fw, fz; // dF/dw, dF/dz at this vertex
		double tval;                 // 1/(|fw|^2+|fz|^2): refinement priority
	};

	vec<3,double> w_sphere(const T& t, Vertex(T) v) { return attr(t,v)->w_sphere; }
	vec<3,double> z_sphere(const T& t, Vertex(T) v) { return attr(t,v)->z_sphere; }

	void w_sphere_set(T& t, Vertex(T) v, const vec<3,double>& p) { attr(t,v)->w_sphere=p; }
	void z_sphere_set(T& t, Vertex(T) v, const vec<3,double>& p) { attr(t,v)->z_sphere=p; }

	int w_label(const T& t, Vertex(T) v) { return attr(t,v)->w_label; }
	int z_label(const T& t, Vertex(T) v) { return attr(t,v)->z_label; }

	void w_label_set(T& t, Vertex(T) v, int i) { attr(t,v)->w_label=i; }
	void z_label_set(T& t, Vertex(T) v, int i) { attr(t,v)->z_label=i; }

	double tval(const T& t, Vertex(T) v) { return attr(t,v)->tval; }

	// Cache of the F=0 intersection point(s) of a 2-simplex, computed
	// lazily by triangle_intersection(). A triangle can genuinely be
	// crossed by more than one sheet (F has degree n in w), so up to 2
	// roots are kept, each identified by its barycentric coordinates
	// (l1,l2), with l0=1-l1-l2.
	template <>
	struct extra_data<2> {
		bool computed;
		int nroots;
		std::complex<double> w_pt[2], z_pt[2];
		double l1[2], l2[2];
		extra_data() : computed(false), nroots(0) {}
	};

}

// Labels for the 6 octahedron vertices: 0, 1, -1, i, -i, infinity.
// 0 <-> south pole, infinity <-> north pole (stereographic projection).
enum { LZERO=0, LONE=1, LMONE=2, LI=3, LMI=4, LINF=5, NLABELS=6 };

const char* label_name(int l) {
	static const char* names[NLABELS] = { "0", "1", "-1", "i", "-i", "inf" };
	return names[l];
}

vec<3,double> label_sphere(int l) {
	vec<3,double> p;
	switch(l) {
		case LZERO: p[0]=0;  p[1]=0;  p[2]=-1; break;
		case LONE:  p[0]=1;  p[1]=0;  p[2]=0;  break;
		case LMONE: p[0]=-1; p[1]=0;  p[2]=0;  break;
		case LI:    p[0]=0;  p[1]=1;  p[2]=0;  break;
		case LMI:   p[0]=0;  p[1]=-1; p[2]=0;  break;
		case LINF:  p[0]=0;  p[1]=0;  p[2]=1;  break;
	}
	return p;
}

// The 8 octahedron triangles, as ordered vertex-label triples plus
// orientation sign, checked for face-compatibility (d_i s1 = d_j s2 => i=j)
// and read off directly from the seed the method was specified with:
//   +<0,1,i>  -<0,-1,i>  +<0,-1,-i>  -<0,1,-i>
//   +<inf,1,i> -<inf,-1,i> +<inf,-1,-i> -<inf,1,-i>
struct octa_tri { int v[3]; int sign; };

// North hemisphere (infinity) signs flipped relative to the original
// seed: south and north triangles sharing an equatorial edge must be
// oppositely oriented for the octahedron itself to be coherent.
const octa_tri octahedron[8] = {
	{ {LZERO, LONE,  LI },  +1 },
	{ {LZERO, LMONE, LI },  -1 },
	{ {LZERO, LMONE, LMI }, +1 },
	{ {LZERO, LONE,  LMI }, -1 },
	{ {LINF,  LONE,  LI },  -1 },
	{ {LINF,  LMONE, LI },  +1 },
	{ {LINF,  LMONE, LMI }, -1 },
	{ {LINF,  LONE,  LMI }, +1 },
};

// The 6 monotone paths across a 3x3 grid (p=q=2), as sequences of
// grid points (i,j) from (0,0) to (2,2), with the Eilenberg-Zilber
// shuffle sign already worked out for each.
struct path6 { int i[5]; int j[5]; int sign; };

const path6 paths[6] = {
	// WWZZ
	{ {0,1,2,2,2}, {0,0,0,1,2}, +1 },
	// WZWZ
	{ {0,1,1,2,2}, {0,0,1,1,2}, -1 },
	// WZZW
	{ {0,1,1,1,2}, {0,0,1,2,2}, +1 },
	// ZWWZ
	{ {0,0,1,2,2}, {0,1,1,1,2}, +1 },
	// ZWZW
	{ {0,0,1,1,2}, {0,1,1,2,2}, -1 },
	// ZZWW
	{ {0,0,0,1,2}, {0,1,2,2,2}, +1 },
};

// Counts internal 3-simplex facets whose two incident 4-simplices do
// not cancel (i.e. are not coherently oriented).
int count_incoherent(const T& t) {
	int nincoherent=0;
	Simplex_it(T,3) i,end;
	for(simplices(t,i,end); i!=end; ++i) {
		Simplex(T,3) f=*i;
		if(!is_current(t,f) || boundary(t,f)) continue;
		pair<Cell(T),Cell(T)> cs=cells(t,f);
		Cell(T) c0=cs.first, c1=cs.second;
		int sign0=1,sign1=1;
		for(int local=0; local<=DIM; ++local) {
			if(face_op(t,c0,local)==f) sign0=(local%2==0)?1:-1;
			if(face_op(t,c1,local)==f) sign1=(local%2==0)?1:-1;
		}
		if(orientation(t,c0)*sign0 + orientation(t,c1)*sign1 != 0) ++nincoherent;
	}
	return nincoherent;
}

// poly_F, cx and the selectable curve catalog live in functions.hpp so
// new test curves can be added there without touching this file. g_F
// is set from function_catalog()[index] in main(), based on --function.
poly_F g_F;
double g_imz_min=1e300, g_imz_max=-1e300;

cx F_raw(cx w, cx z) {
	cx r(1,0);
	for(int i=g_F.n-1;i>=0;--i) r=r*w+g_F.f(i,z);
	return r;
}
cx Fw_raw(cx w, cx z) {
	cx r=double(g_F.n);
	for(int i=g_F.n-1;i>=1;--i) r=r*w+double(i)*g_F.f(i,z);
	return r;
}
cx Fz_raw(cx w, cx z) {
	cx r(0,0);
	for(int i=g_F.n-1;i>=0;--i) r=r*w+g_F.fprime(i,z);
	return r;
}

// Optional complex-affine change of coordinates, F~(w,z)=F(a w+b, c z+d)
// with |a|=|c|=1 (a rotation+translation of the mesh relative to the
// curve, not of the curve's shape). Off (a=c=1, b=d=0) by default.
// Investigating the "duplicate crossing on a shared mesh edge" failure
// mode showed it concentrates where the curve's own symmetry lines up
// with the octahedral seed mesh's (real coefficients, seed vertices on
// the real/imaginary axes); turning this on breaks that alignment
// without changing the curve being triangulated, to compare the two.
//
// This only ever moves the *finite* landmark vertices (0,1,-1,i,-i),
// though: an affine map always fixes infinity, so the seed's w=infty
// and z=infty vertices can never be displaced this way. For a curve
// whose z=infty is a genuine branch point (e.g. any degree-2 F -- by
// Riemann-Hurwitz a degree-n cover of P^1 needs branch points summing
// to a fixed count, so a curve like w^2=z with only one *finite*
// branch point is forced to have a second one at infinity), that seed
// vertex sits exactly on it forever, at every refinement depth (Maubach
// bisection never touches an original seed vertex) -- confirmed
// directly: branch_gap_corner(zr)->0 as zr->0, and the affected cells'
// failure rate stays flat (~38-40%) at every depth tried, instead of
// shrinking like the ordinary alignment-sensitivity case does. See
// g_generic_rot below for the fix.
bool g_generic_coords=false;
cx GA(1,0), GB(0,0), GC(1,0), GD(0,0);

// Sphere rotation applied to the *seed octahedron's* vertex placement
// (not to F) under --generic, independently for the w- and z-factor --
// a rotation is itself a (unitary) Mobius transform of C_infty, and
// unlike the affine one above it does NOT fix infinity: the vertex
// that used to sit at the literal pole now sits at some ordinary
// finite point, generic relative to the curve, so it no longer pins a
// permanently-coarse corner there. Everything downstream (from_sphere,
// the corner chart, sphere_is_far, spherical-midpoint bisection) reads
// off actual sphere coordinates and needs no awareness of this --
// only the initial vertex construction does. F itself is deliberately
// left untouched here (composing F_raw with a true Mobius map, rather
// than just moving mesh vertices, would introduce a *new* pole of F at
// some other finite mesh point -- see the discussion that led here).
vec<3,double> rotate3(const vec<3,double>& p, const vec<3,double>& axis, double c, double s) {
	double dot=p[0]*axis[0]+p[1]*axis[1]+p[2]*axis[2];
	vec<3,double> cr;
	cr[0]=axis[1]*p[2]-axis[2]*p[1];
	cr[1]=axis[2]*p[0]-axis[0]*p[2];
	cr[2]=axis[0]*p[1]-axis[1]*p[0];
	vec<3,double> r;
	for(int i=0;i<3;++i) r[i]=p[i]*c+cr[i]*s+axis[i]*dot*(1-c);
	return r;
}
struct sphere_rot {
	vec<3,double> axis; double c,s;
	sphere_rot() { axis[0]=0; axis[1]=0; axis[2]=1; c=1; s=0; }
	void set(double ax,double ay,double az,double angle) {
		double n=sqrt(ax*ax+ay*ay+az*az);
		axis[0]=ax/n; axis[1]=ay/n; axis[2]=az/n;
		c=cos(angle); s=sin(angle);
	}
	vec<3,double> apply(const vec<3,double>& p) const { return rotate3(p,axis,c,s); }
};
sphere_rot g_rot_w, g_rot_z; // set from main() when --generic is passed

// Random SO(3) rotation (uniform axis via 3 Gaussians, uniform angle),
// replacing the earlier fixed/hardcoded g_rot_w.set(0.3,0.5,0.8,0.9)
// etc. -- same purpose (some rotation that isn't axis-aligned and
// differs between w and z), but now seed-controlled via --generic-seed
// so the "does this actually fix it, or was it luck of one particular
// rotation" question (raised for the analogous CP^2 --generic in
// examples/top/riemann_PC2/) can be asked here too.
unsigned g_generic_seed=12345;
sphere_rot random_sphere_rot(std::mt19937& rng) {
	std::normal_distribution<double> nd(0.0,1.0);
	std::uniform_real_distribution<double> ud(0.0,2.0*std::acos(-1.0));
	sphere_rot r;
	r.set(nd(rng),nd(rng),nd(rng),ud(rng));
	return r;
}

cx F(cx w, cx z)  { return F_raw(GA*w+GB, GC*z+GD); }
cx Fw(cx w, cx z) { return GA*Fw_raw(GA*w+GB, GC*z+GD); }
cx Fz(cx w, cx z) { return GC*Fz_raw(GA*w+GB, GC*z+GD); }

// Determinant of a complex square matrix by Gaussian elimination with
// partial pivoting (only used for small Sylvester matrices below).
cx complex_det(vector<vector<cx> > m) {
	int n=(int)m.size();
	cx det(1,0);
	for(int col=0;col<n;++col) {
		int piv=col; double best=std::abs(m[col][col]);
		for(int r=col+1;r<n;++r) { double v=std::abs(m[r][col]); if(v>best) { best=v; piv=r; } }
		if(best<1e-300) return cx(0,0);
		if(piv!=col) { std::swap(m[piv],m[col]); det=-det; }
		det*=m[col][col];
		for(int r=col+1;r<n;++r) {
			cx factor=m[r][col]/m[col][col];
			for(int cc=col;cc<n;++cc) m[r][cc]-=factor*m[col][cc];
		}
	}
	return det;
}

// Resultant of p (degree dp, highest-degree coefficient first) and q
// (degree dq, likewise), via the (dp+dq)x(dp+dq) Sylvester matrix.
cx resultant(const vector<cx>& p, int dp, const vector<cx>& q, int dq) {
	int n=dp+dq;
	vector<vector<cx> > m(n,vector<cx>(n,cx(0,0)));
	for(int r=0;r<dq;++r) for(int cc=0;cc<=dp;++cc) m[r][r+cc]=p[cc];
	for(int r=0;r<dp;++r) for(int cc=0;cc<=dq;++cc) m[dq+r][r+cc]=q[cc];
	return complex_det(m);
}

// Estimated smallest pairwise gap between the n roots of F(.,z) in w,
// from Res(F(.,z),Fw(.,z)) = leading_coeff^{2n-2} * prod_{i<j}(r_i-r_j)^2
// (Res = -Disc for monic F; the sign doesn't matter here). Taking the
// (n(n-1))-th root of |Res| is exact for n=2 (a single pair) and an
// order-of-magnitude proxy for n>2 (product of several gaps, not just
// the smallest) -- good enough to drive refinement without solving for
// the roots themselves. Generalizes to any degree, unlike a hand-picked
// vertex-local formula.
double branch_gap(cx z) {
	if(g_F.n<2) return 1e18;
	z=GC*z+GD; // gaps between w-roots are invariant under the |a|=1 part
	vector<cx> p(g_F.n+1), q(g_F.n);
	p[0]=cx(1,0);
	for(int i=0;i<g_F.n;++i) p[g_F.n-i]=g_F.f(i,z);
	for(int i=1;i<=g_F.n;++i) {
		cx coeff = (i==g_F.n) ? cx(g_F.n,0) : double(i)*g_F.f(i,z);
		q[g_F.n-i]=coeff;
	}
	cx res=resultant(p,g_F.n,q,g_F.n-1);
	double mag=std::abs(res);
	return std::pow(mag, 1.0/(g_F.n*(g_F.n-1)));
}

// Inverse stereographic projection, sphere -> C_infty. Blows up at (is a
// placeholder for) the north pole itself; this chart is only used once a
// triangle is known not to need the corner treatment below (see
// triangle_intersection / sphere_is_far).
cx from_sphere(const vec<3,double>& p) {
	double denom=1.0-p[2];
	if(fabs(denom)<1e-12) return cx(1e8,0);
	return cx(p[0]/denom,p[1]/denom);
}

// --- Projective treatment of the (w,z) -> (infty,infty) corner --------
// F is monic in w, so w=infty is never a root of F for finite z: the
// w-homogenized w'^n F(1/w',z) equals 1 at w'=0. So a vertex near just
// *one* of the two poles (w far but z finite, or vice versa) is a region
// this curve family's surface never actually visits -- such triangles
// stay excluded from extraction, as before. The only place a second
// chart is actually needed is the *joint* corner where w and z diverge
// together, which is exactly where every catalog curve closes up (e.g.
// the elliptic curve's point(s) at infinity).

// 1/w, directly from sphere coordinates (finite everywhere except at the
// south pole w=0 -- the complementary chart to from_sphere()).
cx recip_from_sphere(const vec<3,double>& p) {
	double denom=1.0+p[2];
	if(fabs(denom)<1e-12) return cx(1e8,0); // w=0 => 1/w=infty placeholder
	return cx(p[0]/denom,-p[1]/denom);
}

// The --generic affine transform doesn't commute with a plain 1/w (the
// translation GB blows the homogenization up), so the "far" variable
// reciprocates the *already-transformed* coordinate instead:
// wr = 1/(GA*w+GB) = w'/(GA+GB*w') with w'=1/w -- reduces to w'=1/w when
// GB=0 (the non-generic default), and stays finite even exactly at
// w'=0. from_wr/from_zr are the inverse, used once a root is found in
// this chart, to report it back as an ordinary (possibly huge but
// finite) (w,z) point for downstream code.
cx to_wr(cx wprime) { return wprime/(GA+GB*wprime); }
cx to_zr(cx zprime) { return zprime/(GC+GD*zprime); }
// Clamped to the same practical-infinity scale as the placeholder above:
// once the single-far-only band is allowed through (see
// triangle_intersection), a root can land at a wr/zr that's small but
// not below the placeholder cutoff, recovering a w/z far past 1e8 that
// downstream code (in particular the plain-F residual check) never
// expected. Rescaling (not truncating) keeps the direction, only bounds
// the magnitude to what the rest of the code already treats as "at the
// pole" in practice.
cx from_wr(cx wr) {
	if(std::norm(wr)<1e-24) return cx(1e8,0);
	cx v=(cx(1,0)/wr-GB)/GA;
	double m=std::abs(v);
	if(m>1e8) v*=1e8/m;
	return v;
}
cx from_zr(cx zr) {
	if(std::norm(zr)<1e-24) return cx(1e8,0);
	cx v=(cx(1,0)/zr-GD)/GC;
	double m=std::abs(v);
	if(m>1e8) v*=1e8/m;
	return v;
}

// A vertex is "far" (near the north-pole proxy) when the ordinary chart
// value there is uncomfortably large -- same cutoff and intent as the
// old near_pole()/POLE_CUTOFF, just checked to decide chart selection
// rather than to blanket-exclude: a triangle with one modestly-large-w
// vertex (say |w|~50) and otherwise ordinary vertices is still handled
// perfectly well by the plain chart (see triangle_intersection), so this
// must stay a high cutoff -- lowering it would needlessly route (or
// worse, exclude, if only one of w/z crosses it) triangles the ordinary
// Newton already solves correctly.
const double FAR_W_CUTOFF=1e6;
bool sphere_is_far(const vec<3,double>& p) { return std::abs(from_sphere(p))>FAR_W_CUTOFF; }

cx ipow(cx base, int e) {
	cx r(1,0);
	for(int k=0;k<e;++k) r*=base;
	return r;
}

// Set once in main() from g_F: the highest z-degree among all f_i(z),
// needed to homogenize F in z at the corner.
int g_dz=0;

// f_i(z) reindexed to the corner chart: z^dz * f_i(1/z) with z=1/zr,
// i.e. f_i "reversed" and padded to degree g_dz -- finite at zr=0 (z=infty).
cx fi_corner(int i, cx zr) {
	const vector<cx>& p=g_F.c[i];
	cx r(0,0);
	for(int k=(int)p.size()-1;k>=0;--k) r+=p[k]*ipow(zr,g_dz-k);
	return r;
}
cx fi_corner_dz(int i, cx zr) { // d/dzr of fi_corner
	const vector<cx>& p=g_F.c[i];
	cx r(0,0);
	for(int k=0;k<(int)p.size();++k) {
		int e=g_dz-k;
		if(e>=1) r+=p[k]*double(e)*ipow(zr,e-1);
	}
	return r;
}

// K(wr,zr) = wr^n * zr^dz * F_raw(1/wr,1/zr): the doubly-homogenized
// polynomial for the (w,z)->(infty,infty) corner, finite and
// well-defined at (wr,zr)=(0,0). Its coefficients are exactly g_F.c[][],
// just reindexed -- no need to re-expand the --generic affine transform
// symbolically, since wr/zr already absorb it (see to_wr/to_zr above,
// which reciprocate GA*w+GB and GC*z+GD, not w and z directly).
cx F_corner(cx wr, cx zr) {
	cx r=ipow(zr,g_dz);
	for(int i=0;i<g_F.n;++i) r+=fi_corner(i,zr)*ipow(wr,g_F.n-i);
	return r;
}
cx Fwr_corner(cx wr, cx zr) { // dK/dwr
	cx r(0,0);
	for(int i=0;i<g_F.n;++i) {
		int e=g_F.n-i;
		if(e>=1) r+=fi_corner(i,zr)*double(e)*ipow(wr,e-1);
	}
	return r;
}
cx Fzr_corner(cx wr, cx zr) { // dK/dzr
	cx r(0,0);
	if(g_dz>=1) r=double(g_dz)*ipow(zr,g_dz-1);
	for(int i=0;i<g_F.n;++i) r+=fi_corner_dz(i,zr)*ipow(wr,g_F.n-i);
	return r;
}

// Corner-chart analogue of branch_gap(): smallest pairwise gap between
// the n roots of K(.,zr) in wr, via the same resultant technique (K's
// coefficients as a polynomial in wr, highest degree first, are exactly
// fi_corner(0..n-1,zr) then the zr^dz constant term; its wr-derivative's
// likewise from differentiating each wr^(n-i) term). Needed because
// branch_gap(z) itself is a *raw w-space* gap, which genuinely diverges
// as z->infty for a curve like w^2=z^3-z (the two roots are ~+-z^1.5
// apart) -- using it unchanged this close to the corner would keep
// cell_priority() artificially tiny there (huge mingap) and starve
// refinement right where the corner chart above needs it most, even
// though the wr-sphere diameter is already bounded.
double branch_gap_corner(cx zr) {
	if(g_F.n<2) return 1e18;
	vector<cx> p(g_F.n+1), q(g_F.n);
	for(int i=0;i<g_F.n;++i) p[i]=fi_corner(i,zr);
	p[g_F.n]=ipow(zr,g_dz);
	for(int i=0;i<g_F.n;++i) q[i]=fi_corner(i,zr)*double(g_F.n-i);
	cx res=resultant(p,g_F.n,q,g_F.n-1);
	double mag=std::abs(res);
	return std::pow(mag,1.0/(g_F.n*(g_F.n-1)));
}

void compute_vertex_data(T& t, Vertex(T) v) {
	cx w=from_sphere(w_sphere(t,v));
	cx z=from_sphere(z_sphere(t,v));
	attr(t,v)->fw=Fw(w,z);
	attr(t,v)->fz=Fz(w,z);
	// Cached local branch-gap estimate (small = close to a branch
	// point); cell_priority() below combines this with the cell's own
	// w-extent, since "close to branch point" alone doesn't say the
	// cell's triangles are actually at risk of holding >1 sheet.
	// Near the corner (z far), use the corner-chart gap instead of the
	// raw one -- see branch_gap_corner's comment.
	if(sphere_is_far(z_sphere(t,v))) {
		cx zr=to_zr(recip_from_sphere(z_sphere(t,v)));
		attr(t,v)->tval=branch_gap_corner(zr);
		return;
	}
	attr(t,v)->tval=branch_gap(z);
}

// Spherical midpoint (exact SLERP at t=1/2): sum the two unit vectors
// and renormalize. Degenerates when a,b are antipodal.
vec<3,double> slerp_midpoint(const vec<3,double>& a, const vec<3,double>& b) {
	vec<3,double> s;
	for(int i=0;i<3;++i) s[i]=a[i]+b[i];
	double n=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
	if(n<1e-9) {
		cerr<<"warning: near-antipodal edge in spherical bisection, "
				<<"falling back to an arbitrary perpendicular point"<<endl;
		vec<3,double> perp;
		perp[0]=-a[1]; perp[1]=a[0]; perp[2]=0;
		double pn=sqrt(perp[0]*perp[0]+perp[1]*perp[1]+perp[2]*perp[2]);
		if(pn<1e-9) { perp[0]=1; perp[1]=0; perp[2]=0; pn=1; }
		for(int i=0;i<3;++i) perp[i]/=pn;
		return perp;
	}
	vec<3,double> r;
	for(int i=0;i<3;++i) r[i]=s[i]/n;
	return r;
}

// Maubach ApplyNew: computes the spherical-midpoint rule for new
// vertices and collects the new top cells created by one subdivide
// call (both new vertices and new cells pass through here, since
// maubach_edge_split applies this object to every simplex in the
// buffer of newly created entities -- add_edge_vertex is the
// vertex-specific hook, apply(t,Cell(T)) catches the new 4-simplices).
struct refine_app : do_nothing {
	using do_nothing::apply;
	vector<Cell(T)> new_cells;

	Vertex(T) add_edge_vertex(T& t, Edge(T) e) {
		vgtl::array<Vertex(T),2> vs;
		vertices(t,e,vs);
		Vertex(T) v=add(t);
		w_sphere_set(t,v,slerp_midpoint(w_sphere(t,vs[0]),w_sphere(t,vs[1])));
		z_sphere_set(t,v,slerp_midpoint(z_sphere(t,vs[0]),z_sphere(t,vs[1])));
		w_label_set(t,v,-1);
		z_label_set(t,v,-1);
		compute_vertex_data(t,v);
		return v;
	}

	void apply(T& t, Cell(T) s) {
		new_cells.push_back(s);
	}
};

double sphere_dist(const vec<3,double>& a, const vec<3,double>& b) {
	double s=0;
	for(int i=0;i<3;++i) { double d=a[i]-b[i]; s+=d*d; }
	return sqrt(s);
}

// t_sigma = (cell's own diameter in w) / (smallest cached branch-gap
// among its vertices). Large means "this cell is wide relative to how
// close the sheets get here" -- the actual risk factor for a triangle
// inside it seeing more than one root -- rather than just "some vertex
// is near a branch point" (which says nothing about the cell's size).
// Diameter is measured on the w-sphere (chordal distance, bounded in
// [0,2]), not in the from_sphere() chart -- unlike a chart diameter, this
// stays meaningful right up to a pole, so cells near one no longer need
// to be excluded here (see sphere_is_far/the corner chart above for how
// extraction itself handles them).
// A cell whose 5 vertices straddle the sphere_is_far cutoff (some past
// it, some not, in either w or z) gets top priority regardless of
// branch_gap: triangle_intersection picks its chart per-triangle from
// that same cutoff, so a straddling cell has some facets in the
// ordinary chart and others in the corner chart K(wr,zr) -- which,
// being the doubly-homogenized polynomial, can carry spurious
// components/roots along wr=0 or zr=0 that aren't part of the actual
// affine curve. Mixing those with genuine ordinary-chart crossings
// inside one cell breaks the per-tetrahedron crossing count even far
// from any branch point (empirically: this is where nearly all
// "bad" cells landed once the single-far band was let through, see
// the corner-chart fix above, and the count grew rather than shrank
// with depth -- not a sampling limit, so branch_gap-based priority
// alone never drives refinement here). Forcing it to subdivide until
// every vertex lands on one side removes the mixed-chart cell itself;
// only a stopgap until chart selection is made cell-consistent instead
// of per-triangle (see cell_priority's caller for the TODO).
bool straddles_far(const T& t, const vgtl::array<Vertex(T),DIM+1>& vs) {
	bool wf=false, wn=false, zf=false, zn=false;
	for(int i=0;i<=DIM;++i) {
		if(sphere_is_far(w_sphere(t,vs[i]))) wf=true; else wn=true;
		if(sphere_is_far(z_sphere(t,vs[i]))) zf=true; else zn=true;
	}
	return (wf&&wn)||(zf&&zn);
}

double cell_priority(const T& t, Cell(T) cv) {
	vgtl::array<Vertex(T),DIM+1> vs;
	vertices(t,cv,vs);
	if(straddles_far(t,vs)) return 1e18;
	double mingap=tval(t,vs[0]);
	for(int i=0;i<=DIM;++i) if(tval(t,vs[i])<mingap) mingap=tval(t,vs[i]);
	double wdiam=0;
	for(int i=0;i<=DIM;++i)
		for(int j=i+1;j<=DIM;++j) {
			double d=sphere_dist(w_sphere(t,vs[i]),w_sphere(t,vs[j]));
			if(d>wdiam) wdiam=d;
		}
	return wdiam/(mingap+1e-12);
}

// Real inner product of two complex numbers viewed as vectors in R^2
// (Re a*Re b + Im a*Im b) -- used below to build a genuinely orthogonal
// real frame for a triangle living in C^2 = R^4 (dot product of a
// (w,z)-pair is this applied to each component and summed).
double real_dot(cx a, cx b) { return a.real()*b.real()+a.imag()*b.imag(); }

// Per-triangle local frame: x runs along the triangle's longest side, y
// along the perpendicular dropped from the opposite vertex. Because the
// side is chosen to be the *longest*, the foot of that perpendicular
// (Po) lands strictly inside it (the angles at its two endpoints are the
// triangle's two smallest, hence acute) -- so this is a well-defined,
// genuinely orthogonal decomposition of the triangle's affine plane,
// with both basis directions scaled to the triangle's own extent. This
// is used only to condition the Newton iteration below: raw barycentric
// edge vectors (w1-w0,z1-z0),(w2-w0,z2-z0) can be near-parallel/skewed
// for a needle-shaped triangle and ill-condition its Jacobian through
// the parametrization alone, independent of F. (It is not a complex/
// holomorphic reparametrization -- see the discussion that ruled that
// out -- just a better-conditioned real one.)
struct tri_frame {
	int io, ip, iq;     // indices (0,1,2) of the opposite / longest-side vertices
	cx Pow, Poz;        // foot of the perpendicular, Po
	cx Xw, Xz, Yw, Yz;  // frame directions: Xdir=Pq-Po, Ydir=Popp-Po
	double xp;          // x-coordinate of vertex p (the other long-side endpoint)
};

void
build_tri_frame(const cx w[3], const cx z[3], tri_frame& fr) {
	double d01=std::norm(w[1]-w[0])+std::norm(z[1]-z[0]);
	double d12=std::norm(w[2]-w[1])+std::norm(z[2]-z[1]);
	double d02=std::norm(w[2]-w[0])+std::norm(z[2]-z[0]);
	if(d01>=d12 && d01>=d02)      { fr.io=2; fr.ip=0; fr.iq=1; }
	else if(d12>=d01 && d12>=d02) { fr.io=0; fr.ip=1; fr.iq=2; }
	else                          { fr.io=1; fr.ip=0; fr.iq=2; }
	cx dw=w[fr.iq]-w[fr.ip], dz=z[fr.iq]-z[fr.ip];
	cx ow=w[fr.io]-w[fr.ip], oz=z[fr.io]-z[fr.ip];
	double denom=real_dot(dw,dw)+real_dot(dz,dz);
	double t=(denom>1e-300) ? (real_dot(ow,dw)+real_dot(oz,dz))/denom : 0.5;
	fr.Pow=w[fr.ip]+t*dw; fr.Poz=z[fr.ip]+t*dz;
	fr.Xw=w[fr.iq]-fr.Pow; fr.Xz=z[fr.iq]-fr.Poz;
	fr.Yw=w[fr.io]-fr.Pow; fr.Yz=z[fr.io]-fr.Poz;
	double denom2=1.0-t;
	fr.xp=(fabs(denom2)>1e-9) ? -t/denom2 : -1e9;
}

// Barycentric (w.r.t. the triangle's original vertex order 0,1,2) ->
// frame (x,y): both are affine parametrizations of the same plane, so
// x,y are themselves affine (in fact linear) in l0,l1,l2, fixed by their
// values at the 3 vertices (x,y)=(0,1) at "o", (xp,0) at "p", (1,0) at "q".
void
bary_to_xy(const tri_frame& fr, double l1, double l2, double& x, double& y) {
	double l[3]={1.0-l1-l2,l1,l2};
	y=l[fr.io];
	x=l[fr.ip]*fr.xp+l[fr.iq];
}

typedef cx (*cxfun2)(cx,cx);

// Newton's method for Ffun(w(x,y),z(x,y))=0 in the triangle's local
// orthogonal frame (see tri_frame above), from a given seed. Ffun/Fwfun/
// Fzfun default to the ordinary F/Fw/Fz, but triangle_intersection passes
// F_corner/Fwr_corner/Fzr_corner (and correspondingly (w,z)-shaped but
// actually-(wr,zr)-valued frame data) for a triangle at the (w,z) ->
// (infty,infty) corner -- Newton itself doesn't need to know which chart
// it's iterating in. Converts back to barycentric coordinates (w.r.t.
// the original vertex order) at the end, so callers and the domain/
// boundary check are unaffected by the frame or chart used internally.
bool
newton_on_triangle(const tri_frame& fr, double x, double y, double& out_l1, double& out_l2,
										cxfun2 Ffun=F, cxfun2 Fwfun=Fw, cxfun2 Fzfun=Fz) {
	for(int iter=0; iter<20; ++iter) {
		cx w=fr.Pow+x*fr.Xw+y*fr.Yw;
		cx z=fr.Poz+x*fr.Xz+y*fr.Yz;
		cx Fv=Ffun(w,z);
		if(std::norm(Fv)<1e-24) break;
		cx fw=Fwfun(w,z), fz=Fzfun(w,z);
		cx cX=fw*fr.Xw+fz*fr.Xz;
		cx cY=fw*fr.Yw+fz*fr.Yz;
		double j00=cX.real(), j01=cY.real(), j10=cX.imag(), j11=cY.imag();
		double jdet=j00*j11-j01*j10;
		if(fabs(jdet)<1e-300) return false;
		double g0=Fv.real(), g1=Fv.imag();
		x-=(g0*j11-j01*g1)/jdet;
		y-=(j00*g1-g0*j10)/jdet;
		if(fabs(x)>1e6||fabs(y)>1e6) return false; // diverged
	}
	cx w=fr.Pow+x*fr.Xw+y*fr.Yw;
	cx z=fr.Poz+x*fr.Xz+y*fr.Yz;
	if(std::norm(Ffun(w,z))>1e-20) return false;

	double l[3];
	l[fr.io]=y;
	double denom=fr.xp-1.0;
	double lp=(fabs(denom)>1e-12) ? (x-(1.0-y))/denom : 0.0;
	l[fr.ip]=lp;
	l[fr.iq]=1.0-y-lp;

	const double dom_tol=1e-6;
	if(l[0]<-dom_tol||l[0]>1+dom_tol||l[1]<-dom_tol||l[1]>1+dom_tol||l[2]<-dom_tol||l[2]>1+dom_tol) return false;
	out_l1=l[1]; out_l2=l[2];
	return true;
}

// Finds where F=0 crosses a 2-simplex: up to two points, since F has
// degree n in w and a triangle's w-span can in principle still contain
// more than one root even after refinement (the discriminant-based
// criterion in branch_gap() bounds how close two roots of F(.,z) get,
// not whether a triangle happens to straddle two of them). Each root
// is sought by Newton, seeded first from the closed-form solution of
// the *linear* interpolation of F at the 3 vertices (exact for a
// linear field, the affine analogue of edge-crossing linear
// interpolation in marching tetrahedra), then from a fixed set of
// seeds spread over the triangle -- both to catch a second sheet the
// linear seed doesn't lead to, and as a fallback when the linear seed
// itself is degenerate or diverges. Distinct converged roots (up to 2)
// are cached per-triangle, since a triangle is examined from both
// tetrahedra that contain it.
void
triangle_intersection(T& t, Simplex(T,2) tri) {
	extra_data<2>* d=attr(t,tri);
	if(d->computed) return;
	d->computed=true;

	vgtl::array<Vertex(T),3> vs;
	vertices(t,tri,vs);

	// Which chart this triangle needs. F is monic in w, so at the exact
	// point z=infty every root w is also infinite -- but w and z don't
	// reach the *same* fixed FAR_W_CUTOFF at the same point along the
	// curve unless w grows linearly in z: e.g. on the parabola
	// w^2=z, |w|~sqrt(|z|), so z crosses the cutoff while w is still
	// only ~sqrt(FAR_W_CUTOFF) -- a real, curve-traversed band near the
	// corner where exactly one of the two is "far". The old code treated
	// w_far!=z_far as the (literally-at-infinity-only) region "w far,
	// z finite" the surface never visits, and excluded it outright --
	// silently dropping that entire band (confirmed: the excluded-triangle
	// count *grows* with refinement depth instead of shrinking, since
	// finer triangles sample more of the band). Fix: use the corner chart
	// whenever *either* is far, not only when both are. This costs
	// nothing extra for a truly-empty triangle (Newton just finds no
	// root there, same outcome as the old exclusion) and is exact where
	// the curve actually passes through this band, since K(wr,zr) is a
	// valid chart (1/w, 1/z are finite and well-conditioned) whether or
	// not the *other* coordinate happens to be large too.
	bool w_far=false, z_far=false;
	for(int k=0;k<3;++k) {
		if(sphere_is_far(w_sphere(t,vs[k]))) w_far=true;
		if(sphere_is_far(z_sphere(t,vs[k]))) z_far=true;
	}

	cx w0,w1,w2,z0,z1,z2; // either the ordinary (w,z) chart, or (wr,zr)
	cxfun2 Ffun,Fwfun,Fzfun;
	if(w_far||z_far) {
		w0=to_wr(recip_from_sphere(w_sphere(t,vs[0])));
		w1=to_wr(recip_from_sphere(w_sphere(t,vs[1])));
		w2=to_wr(recip_from_sphere(w_sphere(t,vs[2])));
		z0=to_zr(recip_from_sphere(z_sphere(t,vs[0])));
		z1=to_zr(recip_from_sphere(z_sphere(t,vs[1])));
		z2=to_zr(recip_from_sphere(z_sphere(t,vs[2])));
		Ffun=F_corner; Fwfun=Fwr_corner; Fzfun=Fzr_corner;
	} else {
		w0=from_sphere(w_sphere(t,vs[0]));
		w1=from_sphere(w_sphere(t,vs[1]));
		w2=from_sphere(w_sphere(t,vs[2]));
		z0=from_sphere(z_sphere(t,vs[0]));
		z1=from_sphere(z_sphere(t,vs[1]));
		z2=from_sphere(z_sphere(t,vs[2]));
		Ffun=F; Fwfun=Fw; Fzfun=Fz;
	}

	double seeds[8][2]; int nseeds=0;
	{
		cx F0=Ffun(w0,z0), F1=Ffun(w1,z1), F2=Ffun(w2,z2);
		cx a1=F1-F0, a2=F2-F0;
		double m00=a1.real(), m01=a2.real(), m10=a1.imag(), m11=a2.imag();
		double det=m00*m11-m01*m10;
		if(fabs(det)>1e-14) {
			seeds[nseeds][0]=(-F0.real()*m11+m01*F0.imag())/det;
			seeds[nseeds][1]=(-m00*F0.imag()+F0.real()*m10)/det;
			++nseeds;
		}
	}
	static const double extra_seeds[7][2]={
		{1/3.,1/3.},{0.1,0.1},{0.8,0.1},{0.1,0.8},{0.45,0.1},{0.1,0.45},{0.45,0.45}
	};
	for(int s=0;s<7;++s) { seeds[nseeds][0]=extra_seeds[s][0]; seeds[nseeds][1]=extra_seeds[s][1]; ++nseeds; }

	tri_frame fr;
	{ cx w3[3]={w0,w1,w2}, z3[3]={z0,z1,z2}; build_tri_frame(w3,z3,fr); }

	for(int s=0; s<nseeds && d->nroots<2; ++s) {
		double sx,sy; bary_to_xy(fr,seeds[s][0],seeds[s][1],sx,sy);
		double ol1,ol2;
		if(!newton_on_triangle(fr,sx,sy,ol1,ol2,Ffun,Fwfun,Fzfun)) continue;
		bool dup=false;
		for(int r=0;r<d->nroots;++r)
			if(fabs(d->l1[r]-ol1)<1e-7 && fabs(d->l2[r]-ol2)<1e-7) dup=true;
		if(dup) continue;
		double ol0=1.0-ol1-ol2;
		int r=d->nroots;
		d->l1[r]=ol1; d->l2[r]=ol2;
		if(w_far||z_far) {
			// The root was found in (wr,zr); convert back to an ordinary
			// (possibly huge but finite) (w,z) point for downstream code
			// (crossing-node identity itself only ever uses l1,l2, which
			// are already chart-independent).
			cx wr_pt=ol0*w0+ol1*w1+ol2*w2;
			cx zr_pt=ol0*z0+ol1*z1+ol2*z2;
			d->w_pt[r]=from_wr(wr_pt);
			d->z_pt[r]=from_zr(zr_pt);
		} else {
			d->w_pt[r]=ol0*w0+ol1*w1+ol2*w2;
			d->z_pt[r]=ol0*z0+ol1*z1+ol2*z2;
		}
		++d->nroots;
	}
}

int triangle_nroots(const T& t, Simplex(T,2) tri) { return attr(t,tri)->nroots; }

void triangle_point(const T& t, Simplex(T,2) tri, int r, cx& w, cx& z) {
	const extra_data<2>* d=attr(t,tri);
	w=d->w_pt[r]; z=d->z_pt[r];
}

// Barycentric coordinates of root r of an already-computed triangle
// (l0=1-l1-l2).
void triangle_bary(const T& t, Simplex(T,2) tri, int r, double& l0, double& l1, double& l2) {
	const extra_data<2>* d=attr(t,tri);
	l1=d->l1[r]; l2=d->l2[r]; l0=1.0-l1-l2;
}

// A crossing point's identity: which mesh element it sits at, purely
// in terms of barycentric coordinates (dimensionless, so scale
// invariant -- unlike a distance in (w,z), which shrinks together with
// the cell size under refinement and so cannot reliably tell "the same
// point, found independently by two adjacent triangles" apart from
// "two genuinely distinct, nearby points"). A vertex (dim=0) or edge
// (dim=1, plus a parameter along it) node arises when the point found
// inside a triangle actually lies on that triangle's boundary; an
// interior point (dim=2) is identified by the triangle and which root
// (0 or 1) it is, so two sheets through the same triangle stay
// distinct.
struct crossing_node {
	int dim, desc, sub;
	double param;
	cx w, z;
};

void
compute_crossing_node(const T& t, Simplex(T,2) tri, int r, crossing_node& nd) {
	double l0,l1,l2;
	triangle_bary(t,tri,r,l0,l1,l2);
	double l[3]={l0,l1,l2};
	const double tol=1e-6;
	int zeros=0, zi[3];
	for(int b=0;b<3;++b) if(fabs(l[b])<tol) zi[zeros++]=b;
	triangle_point(t,tri,r,nd.w,nd.z);
	nd.sub=r; nd.param=0;
	if(zeros>=2) {
		int keep=(zeros==3) ? 0 : (3-zi[0]-zi[1]);
		vgtl::array<Vertex(T),3> vs; vertices(t,tri,vs);
		nd.dim=0; nd.desc=vs[keep].desc;
	} else if(zeros==1) {
		int lo=(zi[0]==0)?1:0, hi=(zi[0]==2)?1:2;
		nd.dim=1; nd.desc=face_op(t,tri,zi[0]).desc;
		nd.param=l[hi]/(l[lo]+l[hi]);
	} else {
		nd.dim=2; nd.desc=tri.desc;
	}
}

bool
same_crossing_node(const crossing_node& a, const crossing_node& b) {
	if(a.dim!=b.dim || a.desc!=b.desc) return false;
	if(a.dim==1 && fabs(a.param-b.param)>1e-6) return false;
	if(a.dim==2 && a.sub!=b.sub) return false;
	return true;
}

// Forward stereographic projection, C_infty -> sphere: the inverse of
// from_sphere(). w=0 -> south pole (0,0,-1), |w|->infty -> north pole
// (0,0,1). Used only by the "onion" visualization mode below -- the
// mesh/Newton/refinement machinery only ever goes the other way.
vec<3,double> to_sphere(cx w) {
	double n=std::norm(w); // |w|^2
	double s=n+1.0;
	vec<3,double> p;
	p[0]=2.0*w.real()/s; p[1]=2.0*w.imag()/s; p[2]=(n-1.0)/s;
	return p;
}

// Radial displacement for the "onion" projection: originally a
// fixed-angle projection of w's own stereographic image onto a generic
// direction, relying on every degree-2 catalog curve being of the form
// w^2=f_0(z) so its two roots are always a +-w pair (odd under
// negation, so any single linear probe already separates them). That
// assumption breaks under --generic: F(w,z)=F_raw(GA w+GB,GC z+GD)
// expands to a monic-in-w quadratic with a nonzero linear term whenever
// GB!=0, so the two roots sum to a fixed nonzero constant instead of 0
// (confirmed: for the current --generic GA/GB, w1+w2 ~= -0.57-0.28i,
// not 0) -- their sphere images are no longer antipodal, so probing one
// root's own position said nothing about how it related to the other,
// and the two shells collapsed into each other instead of separating
// (confirmed visually: --generic onion renders as a nearly featureless
// blob, vs. clear nested-shell cutaways without it).
//
// Fixed by computing the *other* root exactly, via Vieta's on the
// mesh-coordinate quadratic, and comparing the two roots' projections
// directly instead of assuming they're antipodal: this generalizes
// correctly whether or not GB=0, and still gives exactly 0 (no
// separation) right at a branch point, where the two roots coincide.
// (n=2-specific, like the rest of the onion projection so far -- see
// onion_other_root.)
const double ONION_ANGLE=0.83; // an arbitrary non-special angle (radians)
double onion_axis_proj(cx w) {
	vec<3,double> s=to_sphere(w);
	return std::cos(ONION_ANGLE)*s[0]+std::sin(ONION_ANGLE)*s[1]; // in [-1,1]
}
// The companion root of w at the same z: F(w,z)=F_raw(GA w+GB,GC z+GD)
// expands (for a degree-2, monic-in-its-first-argument F_raw) to
// GA^2 w^2 + [2 GA GB + GA f_1(Z)] w + [...], Z=GC z+GD -- so by
// Vieta's, w1+w2 = -(2 GB + f_1(Z))/GA (reduces to the plain -f_1(z)
// when GA=1,GB=0, the non-generic default).
cx onion_other_root(cx w, cx z) {
	cx Z=GC*z+GD;
	cx sum=-(cx(2,0)*GB+g_F.f(1,Z))/GA;
	return sum-w;
}
double onion_radial(cx w, cx z) {
	double t1=onion_axis_proj(w), t2=onion_axis_proj(onion_other_root(w,z));
	return (t1-t2)/2.0; // in [-1,1]; exactly 0 when w coincides with its companion (a branch point)
}

bool g_onion=false;
double g_onion_scale=0.3;

// Visualization cutoff: a polygon with at least one vertex whose
// projected (Re w, Im w, Re z) point lies farther than g_cutoff from
// the origin is dropped from the OBJ output entirely, rather than
// clipped -- the far vertices near a pole proxy (up to the from_sphere
// placeholder scale of 1e8) would otherwise dominate the model's
// bounding box. Off (no filtering) unless --cutoff is passed. Only
// applied in flat (non-onion) mode for now -- onion's radius is always
// close to 1 by construction (see onion_radial), so this cutoff isn't
// meaningful there yet.
double g_cutoff=1e300;

// 4D (w,z) in C_infty^2 -> 3D, for a first look at the extracted
// surface. Two modes:
// - flat (default): (Re w, Im w, Re z). Throws away Im(z) entirely, so
//   it's a projection, not an embedding -- it can and will show
//   self-intersections that aren't really there.
// - onion (--onion): place the point by z's position on the sphere
//   (direction), then displace it perpendicular to that sphere -- i.e.
//   radially -- by onion_radial(w). For a fixed z, the (up to n) sheets
//   of the surface then spread into concentric shells instead of
//   overlapping, showing the branch structure as nested spheres that
//   pinch together where sheets meet. Also a projection, not an
//   embedding (w is compressed to one bounded scalar), but a
//   differently-lossy one.
vec<3,double> project_for_viz(cx w, cx z) {
	if(g_onion) {
		vec<3,double> dir=to_sphere(z);
		double r=1.0+g_onion_scale*onion_radial(w,z);
		vec<3,double> p;
		for(int i=0;i<3;++i) p[i]=r*dir[i];
		return p;
	}
	vec<3,double> p;
	p[0]=w.real(); p[1]=w.imag(); p[2]=z.real();
	return p;
}

// Writes an OBJ file: each extracted cycle (polygon, one per sheet
// crossing a 4-simplex) becomes its own OBJ face, with its own
// (unwelded) vertices -- no attempt is made to merge vertices shared
// between polygons from adjacent 4-simplices, since the extraction
// doesn't currently track that adjacency. Fine for a first visual
// check; a real output mesh (with a shared vertex pool and welding)
// is a separate step once the projection choice settles down.
struct obj_writer {
	ofstream out;
	int nverts, nfaces, nclipped;
	obj_writer(const char* path) : out(path), nverts(0), nfaces(0), nclipped(0) {
		if(g_onion)
			out<<"# Riemann surface extraction, onion projection: "
				<<"direction=z on S^2, radius=1+"<<g_onion_scale<<"*onion_radial(w)\n";
		else
			out<<"# Riemann surface extraction, projected via (Re w, Im w, Re z)\n";
		if(!g_onion && g_cutoff<1e299)
			out<<"# cutoff: polygons with a vertex farther than "<<g_cutoff<<" from the origin dropped\n";
	}
	void write_polygon(const vector<crossing_node>& nodes, const vector<int>& cyc) {
		vector<vec<3,double> > pts(cyc.size());
		for(size_t i=0;i<cyc.size();++i)
			pts[i]=project_for_viz(nodes[cyc[i]].w, nodes[cyc[i]].z);
		if(!g_onion) {
			for(size_t i=0;i<pts.size();++i) {
				double d=sqrt(pts[i][0]*pts[i][0]+pts[i][1]*pts[i][1]+pts[i][2]*pts[i][2]);
				if(d>g_cutoff) { ++nclipped; return; }
			}
		}
		for(size_t i=0;i<cyc.size();++i) {
			extern double g_imz_min, g_imz_max;
			double imz=nodes[cyc[i]].z.imag();
			if(imz<g_imz_min) g_imz_min=imz;
			if(imz>g_imz_max) g_imz_max=imz;
			out<<"v "<<pts[i][0]<<" "<<pts[i][1]<<" "<<pts[i][2]<<"\n";
		}
		out<<"f";
		for(size_t i=0;i<cyc.size();++i) out<<" "<<(nverts+(int)i+1);
		out<<"\n";
		nverts+=(int)cyc.size();
		++nfaces;
	}
};

int main(int argc, char* argv[]) {

	int max_depth=12;
	double threshold=0.1;
	int function_index=0;
	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--generic") {
			g_generic_coords=true;
		} else if(arg=="--generic-seed" && i+1<argc) {
			g_generic_coords=true; g_generic_seed=(unsigned)atoi(argv[++i]);
		} else if(arg=="--depth" && i+1<argc) {
			max_depth=atoi(argv[++i]);
		} else if(arg=="--threshold" && i+1<argc) {
			threshold=atof(argv[++i]);
		} else if(arg=="--function" && i+1<argc) {
			function_index=atoi(argv[++i]);
		} else if(arg=="--onion") {
			g_onion=true;
		} else if(arg=="--onion-scale" && i+1<argc) {
			g_onion_scale=atof(argv[++i]);
		} else if(arg=="--cutoff" && i+1<argc) {
			g_cutoff=atof(argv[++i]);
		} else if(arg=="--list-functions") {
			cout<<"available functions:"<<endl;
			print_function_catalog(cout);
			return 0;
		} else {
			cerr<<"unrecognized argument: "<<arg<<endl;
			cerr<<"usage: "<<argv[0]<<" [--function N] [--generic] [--generic-seed N] [--depth N] "
					<<"[--threshold X] [--onion] [--onion-scale X] [--cutoff X] [--list-functions]"<<endl;
			return 1;
		}
	}

	if(g_generic_coords) {
		// GA/GB/GC/GD (the affine transform ON F, see the comment above
		// their declaration) is deliberately left at the identity here --
		// NOT set to a fixed nonzero offset the way it used to be. Reason:
		// it and the rotation below were two separate mechanisms for the
		// same underlying purpose (break the seed mesh's alignment with
		// the curve), added at different times, and only the rotation
		// can reach the w=infty/z=infty vertices (an affine map always
		// fixes infinity). Once the rotation existed, it ALSO already
		// moves every finite landmark (0,1,-1,i,-i) -- not just the
		// poles -- making the affine transform redundant for its
		// original purpose. Confirmed empirically, not just argued:
		// parabola, --depth 10, 4 seeds, rotation-only vs rotation+affine
		// gave statistically indistinguishable bad-cell rates (both
		// ~0.1-1.2%, both >10x better than no --generic at all; neither
		// consistently beat the other across seeds). So this is now done
		// exactly the way examples/top/riemann_PC2/riemann_pc2.cpp does
		// it: ONE mesh-vertex-position transform, not a wrapper re-run on
		// every F evaluation. GA/GB/GC/GD and everything built on them
		// (to_wr/to_zr, onion_other_root) are kept, not deleted -- they
		// already degrade correctly to the plain GA=1,GB=0 formulas at
		// this identity default (that generality was real debugging work,
		// see the comments at to_wr/onion_other_root), and repurposing
		// GA/GB/GC/GD for something else later remains possible.

		// Random, seed-controlled rotation of the seed octahedron itself
		// (independent draws for w and z -- two draws from the same
		// stream, so they generically differ, unlike a shared rotation
		// which would keep the product mesh's own w<->z symmetry).
		std::mt19937 rng(g_generic_seed);
		g_rot_w=random_sphere_rot(rng);
		g_rot_z=random_sphere_rot(rng);
	}

	vector<catalog_entry>& catalog=function_catalog();
	if(function_index<0 || function_index>=(int)catalog.size()) {
		cerr<<"--function "<<function_index<<" out of range; available:"<<endl;
		print_function_catalog(cerr);
		return 1;
	}
	g_F=catalog[function_index].F;
	g_dz=0;
	for(int i=0;i<g_F.n;++i) {
		int deg=(int)g_F.c[i].size()-1;
		if(deg>g_dz) g_dz=deg;
	}
	cout<<"function: "<<function_index<<" ("<<catalog[function_index].name<<") -- "
			<<catalog[function_index].description<<endl;
	cout<<"coordinates: "<<(g_generic_coords?"generic (rotated+translated)":"aligned (original)");
	if(g_generic_coords) cout<<" (seed="<<g_generic_seed<<")";
	cout<<endl;
	cout<<"max_depth="<<max_depth<<" threshold="<<threshold<<endl;
	cout<<"projection: "<<(g_onion?"onion":"flat");
	if(g_onion) cout<<" (scale="<<g_onion_scale<<")";
	cout<<endl;
	if(g_cutoff<1e299)
		cout<<"cutoff: "<<g_cutoff<<(g_onion?" (ignored in onion mode)":"")<<endl;

	// Named after the curve, coordinate mode and projection so runs over
	// different --function/--generic/--onion combinations never silently
	// overwrite each other's output.
	string obj_path_s="riemann_surface_"+catalog[function_index].name
		+(g_generic_coords?"_generic":"")+(g_onion?"_onion":"")+".obj";
	const char* obj_path=obj_path_s.c_str();

	// One combined vertex per (w-label, z-label) pair.
	Vertex(T) V[NLABELS][NLABELS];
	for(int wl=0; wl<NLABELS; ++wl) {
		for(int zl=0; zl<NLABELS; ++zl) {
			Vertex(T) v=add(t);
			w_sphere_set(t,v, g_generic_coords ? g_rot_w.apply(label_sphere(wl)) : label_sphere(wl));
			z_sphere_set(t,v, g_generic_coords ? g_rot_z.apply(label_sphere(zl)) : label_sphere(zl));
			w_label_set(t,v,wl);
			z_label_set(t,v,zl);
			compute_vertex_data(t,v);
			V[wl][zl]=v;
		}
	}

	int ncells=0;
	{
		complex_builder<T> cb(t);
		for(int a=0; a<8; ++a) {
			for(int b=0; b<8; ++b) {
				const octa_tri& sw=octahedron[a]; // w-factor triangle
				const octa_tri& sz=octahedron[b]; // z-factor triangle
				for(int p=0; p<6; ++p) {
					const path6& pa=paths[p];
					vgtl::array<Vertex(T),5> vs;
					for(int k=0; k<5; ++k) {
						int wl=sw.v[pa.i[k]];
						int zl=sz.v[pa.j[k]];
						vs[k]=V[wl][zl];
					}
					// Shuffle-sign formula (eps_w * eps_z * sign(path)); used
					// only as a seed, the BFS pass below re-derives a
					// consistent orientation and overwrites this if needed.
					int ori=sw.sign*sz.sign*pa.sign;
					Cell(T) cv=add(cb,vs);
					orientation_set(t,cv,ori);
					level_set(t,cv,0);
					++ncells;
				}
			}
		}
	}

	cout<<"vertices added: "<<(NLABELS*NLABELS)<<endl;
	cout<<"4-simplices added: "<<ncells<<endl;

	cout<<"incoherent facets from the raw shuffle-sign formula: "
			<<count_incoherent(t)<<" (expected 0 now that the north "
			<<"hemisphere sign is flipped)"<<endl;

	// Kept as a safety net / independent cross-check: re-derive a
	// coherent orientation by BFS over the cell-adjacency graph. If the
	// formula above is already coherent this should change nothing.
	{
		set<Cell(T)> seen;
		queue<Cell(T)> q;
		Cell_it(T) ci,cend;
		simplices(t,ci,cend);
		Cell(T) start=*ci;
		orientation_set(t,start,1);
		seen.insert(start);
		q.push(start);
		int nconflict=0;
		while(!q.empty()) {
			Cell(T) c=q.front(); q.pop();
			for(int j=0; j<=DIM; ++j) {
				Facet(T) f=face_op(t,c,j);
				if(boundary(t,f)) continue;
				Cell(T) c2=adjacent(t,c,j);
				int j2=-1;
				for(int k=0; k<=DIM; ++k) if(face_op(t,c2,k)==f) { j2=k; break; }
				int sign0=(j%2==0)?1:-1;
				int sign1=(j2%2==0)?1:-1;
				int required=-orientation(t,c)*sign0*sign1;
				if(seen.count(c2)) {
					if(orientation(t,c2)!=required) ++nconflict;
					continue;
				}
				orientation_set(t,c2,required);
				seen.insert(c2);
				q.push(c2);
			}
		}
		cout<<"orientation BFS reached "<<seen.size()<<"/"<<ncells<<" cells, "
				<<nconflict<<" conflicts (expected 0 if orientable)"<<endl;
	}

	int nv=0,ne=0,nt=0,ntet=0,n4=0;
	{
		Vertex_it(T) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++nv;
	}
	{
		Simplex_it(T,1) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++ne;
	}
	{
		Simplex_it(T,2) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++nt;
	}
	{
		Simplex_it(T,3) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++ntet;
	}
	{
		Cell_it(T) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++n4;
	}

	cout<<"0-simplices (vertices): "<<nv<<endl;
	cout<<"1-simplices (edges): "<<ne<<endl;
	cout<<"2-simplices (triangles): "<<nt<<endl;
	cout<<"3-simplices (tetrahedra): "<<ntet<<endl;
	cout<<"4-simplices: "<<n4<<endl;

	cout<<"euler characteristic: "<<euler_characteristic(t)
			<<" (expected 4, since C_infty^2 ~ S^2 x S^2)"<<endl;

	int nbnd=0;
	{
		Simplex_it(T,3) i,end;
		for(simplices(t,i,end); i!=end; ++i)
			if(is_current(t,*i) && boundary(t,*i)) ++nbnd;
	}
	cout<<"boundary 3-simplices: "<<nbnd<<" (expected 0, mesh should be closed)"<<endl;

	cout<<"incoherently oriented internal facets (after BFS safety net): "
			<<count_incoherent(t)<<" (expected 0)"<<endl;

	// --- Phase 2: adaptive refinement, F(w,z) = w^2 - z -------------

	cout<<endl<<"--- adaptive refinement ---"<<endl;

	priority_queue<pair<double,Cell(T)> > pq;
	{
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend); ci!=cend; ++ci)
			if(is_current(t,*ci)) pq.push(make_pair(cell_priority(t,*ci),*ci));
	}

	int nsubdivisions=0;
	while(!pq.empty()) {
		pair<double,Cell(T)> top=pq.top(); pq.pop();
		Cell(T) cv=top.second;
		if(!is_current(t,cv)) continue;         // stale queue entry
		if(top.first<threshold) break;          // max-heap: nothing left qualifies
		if(level(t,cv)>=max_depth) continue;     // capped, but keep draining the queue
		refine_app app;
		maubach_subdivide(t,cv,app);
		++nsubdivisions;
		for(size_t i=0;i<app.new_cells.size();++i)
			pq.push(make_pair(cell_priority(t,app.new_cells[i]),app.new_cells[i]));
	}

	cout<<"subdivisions performed: "<<nsubdivisions<<endl;

	int nv2=0,n42=0;
	double tmin=1e300,tmax=0;
	{
		Vertex_it(T) i,end;
		for(simplices(t,i,end); i!=end; ++i) {
			if(!is_current(t,*i)) continue;
			++nv2;
			double tv=tval(t,*i);
			if(tv<tmin) tmin=tv;
			if(tv>tmax) tmax=tv;
		}
	}
	{
		Cell_it(T) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++n42;
	}
	cout<<"vertices after refinement: "<<nv2<<endl;
	cout<<"4-simplices after refinement: "<<n42<<endl;
	cout<<"vertex tval range: ["<<tmin<<", "<<tmax<<"]"<<endl;
	cout<<"euler characteristic after refinement: "<<euler_characteristic(t)
			<<" (expected 4, subdivision must not change the topology)"<<endl;

	int nbnd2=0;
	{
		Simplex_it(T,3) i,end;
		for(simplices(t,i,end); i!=end; ++i)
			if(is_current(t,*i) && boundary(t,*i)) ++nbnd2;
	}
	cout<<"boundary 3-simplices after refinement: "<<nbnd2<<" (expected 0)"<<endl;
	cout<<"incoherent facets after refinement: "<<count_incoherent(t)
			<<" (expected 0; checks that maubach's own orientation rule "
			<<"stays consistent with ours)"<<endl;

	{
		int hist[32]={0};
		Cell_it(T) i,end;
		for(simplices(t,i,end); i!=end; ++i)
			if(is_current(t,*i)) ++hist[level(t,*i)];
		cout<<"final cell level histogram:";
		for(int l=0;l<32;++l) if(hist[l]) cout<<" ["<<l<<"]="<<hist[l];
		cout<<endl;
	}

	// Triangle-shape diagnostic: minimum angle of every current 2-simplex,
	// computed in the ambient R^3 x R^3 = R^6 embedding given by
	// (w_sphere,z_sphere) -- the same space slerp bisection actually
	// operates in (unlike the (w,z) chart, which blows up near the pole
	// proxy). Any 3 points span an affine subspace of dimension <=2, so
	// each triangle is exactly flat there and its angles genuinely sum to
	// 180 deg (min angle <=60 deg). This is an empirical check of the
	// claim (made when conditioning the per-triangle Newton solve) that
	// Maubach bisection keeps simplex shapes within a bounded family
	// rather than degenerating -- not assumed here, measured.
	{
		const double PI=3.14159265358979323846;
		double min_angle=180.0, max_min_angle=0.0, sum_angle=0.0;
		int ntri=0, nsliver5=0, nsliver1=0;
		int hist[18]={0}; // 5-degree buckets, [0,90)
		Simplex_it(T,2) i,end;
		for(simplices(t,i,end); i!=end; ++i) {
			if(!is_current(t,*i)) continue;
			vgtl::array<Vertex(T),3> vs;
			vertices(t,*i,vs);
			double P[3][6];
			for(int k=0;k<3;++k) {
				vec<3,double> ws=w_sphere(t,vs[k]), zs=z_sphere(t,vs[k]);
				for(int c=0;c<3;++c) { P[k][c]=ws[c]; P[k][3+c]=zs[c]; }
			}
			double ang[3];
			for(int k=0;k<3;++k) {
				int a=k, b=(k+1)%3, c=(k+2)%3;
				double u[6],v[6],un=0,vn=0,dot=0;
				for(int d=0;d<6;++d) { u[d]=P[b][d]-P[a][d]; v[d]=P[c][d]-P[a][d]; }
				for(int d=0;d<6;++d) { un+=u[d]*u[d]; vn+=v[d]*v[d]; dot+=u[d]*v[d]; }
				un=sqrt(un); vn=sqrt(vn);
				double cosang=dot/(un*vn);
				if(cosang>1) cosang=1; if(cosang<-1) cosang=-1;
				ang[k]=acos(cosang)*180.0/PI;
			}
			double mn=ang[0]; if(ang[1]<mn) mn=ang[1]; if(ang[2]<mn) mn=ang[2];
			if(mn<min_angle) min_angle=mn;
			if(mn>max_min_angle) max_min_angle=mn;
			sum_angle+=mn;
			++ntri;
			if(mn<5.0) ++nsliver5;
			if(mn<1.0) ++nsliver1;
			int bucket=(int)(mn/5.0); if(bucket>17) bucket=17; if(bucket<0) bucket=0;
			++hist[bucket];
		}
		cout<<endl<<"--- triangle min-angle statistics (w_sphere/z_sphere ambient R^6) ---"<<endl;
		cout<<"triangles: "<<ntri<<endl;
		cout<<"min-angle over all triangles: min="<<min_angle<<" deg, "
				<<"max="<<max_min_angle<<" deg, mean="<<(sum_angle/ntri)<<" deg"<<endl;
		cout<<"slivers: min-angle<5deg: "<<nsliver5<<" ("<<(100.0*nsliver5/ntri)<<"%) "
				<<"min-angle<1deg: "<<nsliver1<<" ("<<(100.0*nsliver1/ntri)<<"%)"<<endl;
		cout<<"histogram (5-degree buckets of the per-triangle min angle):"<<endl;
		for(int b=0;b<18;++b) if(hist[b])
			cout<<"  ["<<(b*5)<<","<<(b*5+5)<<"): "<<hist[b]
					<<" ("<<(100.0*hist[b]/ntri)<<"%)"<<endl;
	}

	// --- Phase 3: surface extraction ---------------------------------
	// Per 4-simplex sigma: find the F=0 crossing(s) of each of its ten
	// 2-faces (cached, so each triangle is solved once even though it
	// is shared by two of sigma's five tetrahedral facets), and
	// identify each crossing point canonically (compute_crossing_node)
	// so that the same physical point found from different triangles
	// -- whether because it lies on a shared edge/vertex, or because a
	// triangle has two roots -- collapses to one node. A tetrahedron
	// contributes an edge of the local polygon between the two nodes
	// on its boundary (or, in the rarer case of 4 nodes -- two sheets
	// crossing the same tetrahedron -- two edges, paired by minimal
	// total (w,z) length); since every crossed triangle of sigma
	// belongs to exactly two of its five tetrahedra, these edges close
	// into cycles of 3 or more points (checked, not assumed -- a sigma
	// with more than one cycle means more than one sheet threads it).

	cout<<endl<<"--- surface extraction ---"<<endl;

	int npoly_out[8]={0,0,0,0,0,0,0,0}; // indexed by polygon size (3..7)
	int ntouching_tets=0;   // tetrahedra grazed at a single point
	int nbad_tets=0;        // tetrahedra with an unhandled node count (not 0,1,2,4)
	int ok_cells=0, bad_cells=0;
	int ok_by_level[32]={0}, bad_by_level[32]={0};
	double fres_min=1e300, fres_max=0;

	obj_writer obj(obj_path);

	{
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend); ci!=cend; ++ci) {
			Cell(T) sigma=*ci;
			if(!is_current(t,sigma)) continue;

			vector<crossing_node> nodes;
			vector<vector<int> > adj;
			bool tet_unhandled=false;

			for(int j=0; j<=DIM; ++j) {
				Simplex(T,3) tet=face_op(t,sigma,j);
				int tnodes[4]; int ntn=0;
				for(int k=0; k<=3; ++k) {
					Simplex(T,2) tri=face_op(t,tet,k);
					triangle_intersection(t,tri);
					int nr=triangle_nroots(t,tri);
					for(int r=0;r<nr;++r) {
						crossing_node nd;
						compute_crossing_node(t,tri,r,nd);
						int idx=-1;
						for(size_t q=0;q<nodes.size();++q)
							if(same_crossing_node(nodes[q],nd)) { idx=(int)q; break; }
						if(idx<0) { idx=(int)nodes.size(); nodes.push_back(nd); adj.push_back(vector<int>()); }
						bool dup=false;
						for(int q=0;q<ntn;++q) if(tnodes[q]==idx) dup=true;
						if(!dup && ntn<4) tnodes[ntn++]=idx;
					}
				}
				if(ntn==2) {
					adj[tnodes[0]].push_back(tnodes[1]);
					adj[tnodes[1]].push_back(tnodes[0]);
				} else if(ntn==4) {
					// Two sheets cross this tetrahedron: pair the 4 points
					// into 2 edges by whichever pairing has the smaller
					// total (w,z) length (the two sheets don't cross each
					// other generically, so the short pairing is correct).
					static const int pairings[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
					int best=0; double bestcost=1e300;
					for(int c=0;c<3;++c) {
						double cost=0;
						for(int e=0;e<2;++e) {
							const crossing_node& A=nodes[tnodes[pairings[c][2*e]]];
							const crossing_node& B=nodes[tnodes[pairings[c][2*e+1]]];
							cost+=sqrt(std::norm(A.w-B.w)+std::norm(A.z-B.z));
						}
						if(cost<bestcost) { bestcost=cost; best=c; }
					}
					for(int e=0;e<2;++e) {
						int a=tnodes[pairings[best][2*e]], b=tnodes[pairings[best][2*e+1]];
						adj[a].push_back(b); adj[b].push_back(a);
					}
				} else if(ntn==1) {
					++ntouching_tets;
				} else if(ntn!=0) {
					++nbad_tets;
					tet_unhandled=true;
				}
			}

			// A node with degree 0 was seen only inside a "touching"
			// tetrahedron (already counted in ntouching_tets) and never
			// paired into an edge by any tetrahedron of sigma -- it isn't
			// part of any polygon and must not affect whether the *other*
			// nodes' cycle(s) are valid, so it's dropped before checking.
			bool deg_ok=true;
			bool any_edge=false;
			for(size_t q=0;q<nodes.size();++q) {
				size_t dg=adj[q].size();
				if(dg==0) continue;
				any_edge=true;
				if(dg!=2) deg_ok=false;
			}
			if(!any_edge) continue;

			vector<vector<int> > cycles_idx;
			bool decompose_ok=deg_ok && !tet_unhandled;
			if(decompose_ok) {
				vector<bool> seen(nodes.size(),false);
				for(size_t q=0;q<nodes.size();++q) {
					if(seen[q] || adj[q].empty()) continue;
					vector<int> cyc;
					int start=(int)q, prev=start, cur=adj[start][0];
					seen[q]=true; cyc.push_back(start);
					while(cur!=start) {
						if(seen[cur]) { decompose_ok=false; break; } // shouldn't happen if deg_ok
						seen[cur]=true; cyc.push_back(cur);
						int nxt=(adj[cur][0]==prev) ? adj[cur][1] : adj[cur][0];
						prev=cur; cur=nxt;
						if(cyc.size()>(size_t)(DIM+1)) { decompose_ok=false; break; } // a 4-simplex has only 5 facets
					}
					if(!decompose_ok) break;
					// A "cycle" of length < 3 is a degenerate digon: two
					// tetrahedra of sigma both connected the same pair of
					// nodes (e.g. every active tetrahedron collapses to the
					// same 2 points because the true crossing is a sliver
					// too thin for this cell to resolve into a real
					// polygon). Every node still has degree exactly 2, so
					// deg_ok alone doesn't catch this -- reject the whole
					// cell rather than emit a degenerate face.
					if(cyc.size()<3) { decompose_ok=false; break; }
					cycles_idx.push_back(cyc);
				}
			}

			if(!decompose_ok) {
				++bad_cells;
				++bad_by_level[level(t,sigma)];
				continue;
			}
			++ok_cells;
			++ok_by_level[level(t,sigma)];

			for(size_t q=0;q<nodes.size();++q) {
				if(adj[q].empty()) continue;
				double r=std::norm(F(nodes[q].w,nodes[q].z));
				if(r<fres_min) fres_min=r;
				if(r>fres_max) fres_max=r;
			}
			for(size_t p=0;p<cycles_idx.size();++p) {
				int sz=(int)cycles_idx[p].size();
				if(sz>=3 && sz<8) ++npoly_out[sz];
				obj.write_polygon(nodes,cycles_idx[p]);
			}
		}
	}

	int ntriangles_out=npoly_out[3];
	int noutput_tris=0;
	for(int sz=3; sz<8; ++sz) noutput_tris+=npoly_out[sz]*(sz-2); // fan triangulation count
	cout<<"extracted cells: ok="<<ok_cells<<" bad="<<bad_cells<<endl;
	cout<<"polygon sizes:";
	for(int sz=3; sz<8; ++sz) if(npoly_out[sz]) cout<<" "<<sz<<"-gon="<<npoly_out[sz];
	cout<<endl;
	cout<<"total output triangles (fan-triangulated): "<<noutput_tris<<endl;
	cout<<"touching tetrahedra (single point, no edge): "<<ntouching_tets<<endl;
	cout<<"tetrahedra with an unhandled node count: "<<nbad_tets<<endl;
	cout<<"[diag] Im(z) range at extracted nodes: ["<<g_imz_min<<", "<<g_imz_max<<"]"<<endl;
	cout<<"wrote "<<obj_path<<": "<<obj.nverts<<" vertices, "<<obj.nfaces<<" faces "
			<<"(projection: "<<(g_onion?"onion, z on S^2 + radial w":"Re(w), Im(w), Re(z)")<<")"<<endl;
	if(obj.nclipped) cout<<"polygons dropped by --cutoff: "<<obj.nclipped<<endl;
	cout<<"|F|^2 residual at extracted nodes, range: ["<<fres_min<<", "<<fres_max<<"]"<<endl;

	cout<<"ok/bad cells by cell level:"<<endl;
	for(int l=0;l<32;++l) if(ok_by_level[l]||bad_by_level[l])
		cout<<"  level "<<l<<": ok="<<ok_by_level[l]<<" bad="<<bad_by_level[l]<<endl;

	return 0;
}
