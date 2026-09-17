#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <queue>
#include <set>
#include <map>
#include <complex>
#include <cmath>
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
bool g_generic_coords=false;
cx GA(1,0), GB(0,0), GC(1,0), GD(0,0);

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

// Inverse stereographic projection, sphere -> C_infty. The north pole
// (point at infinity) is approximated by a large finite value here;
// this is a placeholder, not a projective treatment, and should be
// revisited once refinement actually needs to happen near a pole.
cx from_sphere(const vec<3,double>& p) {
	double denom=1.0-p[2];
	if(fabs(denom)<1e-12) return cx(1e8,0);
	return cx(p[0]/denom,p[1]/denom);
}

// Vertices at/near the 1e8 pole proxy make F and its derivatives
// overflow into meaningless magnitudes (observed: F ~ 1e16, spurious
// "roots" accepted by Newton). Since this example's curve never visits
// a neighborhood of w=infty or z=infty, we simply keep such cells and
// triangles out of both refinement and extraction rather than doing
// real arithmetic on the proxy value. A projective (numerator,
// denominator) treatment would remove the need for this guard.
const double POLE_CUTOFF=1e6;
bool near_pole(cx w) { return std::abs(w)>POLE_CUTOFF; }

void compute_vertex_data(T& t, Vertex(T) v) {
	cx w=from_sphere(w_sphere(t,v));
	cx z=from_sphere(z_sphere(t,v));
	attr(t,v)->fw=Fw(w,z);
	attr(t,v)->fz=Fz(w,z);
	// Cached local branch-gap estimate (small = close to a branch
	// point); cell_priority() below combines this with the cell's own
	// w-extent, since "close to branch point" alone doesn't say the
	// cell's triangles are actually at risk of holding >1 sheet.
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
		array<Vertex(T),2> vs;
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

// t_sigma = (cell's own diameter in w) / (smallest cached branch-gap
// among its vertices). Large means "this cell is wide relative to how
// close the sheets get here" -- the actual risk factor for a triangle
// inside it seeing more than one root -- rather than just "some vertex
// is near a branch point" (which says nothing about the cell's size).
double cell_priority(const T& t, Cell(T) cv) {
	array<Vertex(T),DIM+1> vs;
	vertices(t,cv,vs);
	cx w[DIM+1];
	double mingap=tval(t,vs[0]);
	for(int i=0;i<=DIM;++i) {
		w[i]=from_sphere(w_sphere(t,vs[i]));
		// A cell touching the pole proxy has an unbounded/meaningless
		// w-diameter; refining it further wastes effort and cannot help
		// (see near_pole), so it gets no priority from this criterion.
		if(near_pole(w[i])) return 0;
		if(tval(t,vs[i])<mingap) mingap=tval(t,vs[i]);
	}
	double wdiam=0;
	for(int i=0;i<=DIM;++i)
		for(int j=i+1;j<=DIM;++j) {
			double d=std::abs(w[i]-w[j]);
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

// Newton's method for F(w(x,y),z(x,y))=0 in the triangle's local
// orthogonal frame (see tri_frame above), from a given seed. Converts
// back to barycentric coordinates (w.r.t. the original vertex order) at
// the end, so callers and the domain/boundary check are unaffected by
// the frame used internally.
bool
newton_on_triangle(const tri_frame& fr, double x, double y, double& out_l1, double& out_l2) {
	for(int iter=0; iter<20; ++iter) {
		cx w=fr.Pow+x*fr.Xw+y*fr.Yw;
		cx z=fr.Poz+x*fr.Xz+y*fr.Yz;
		cx Fv=F(w,z);
		if(std::norm(Fv)<1e-24) break;
		cx fw=Fw(w,z), fz=Fz(w,z);
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
	if(std::norm(F(w,z))>1e-20) return false;

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

	array<Vertex(T),3> vs;
	vertices(t,tri,vs);
	cx w0=from_sphere(w_sphere(t,vs[0]));
	cx w1=from_sphere(w_sphere(t,vs[1]));
	cx w2=from_sphere(w_sphere(t,vs[2]));
	cx z0=from_sphere(z_sphere(t,vs[0]));
	cx z1=from_sphere(z_sphere(t,vs[1]));
	cx z2=from_sphere(z_sphere(t,vs[2]));

	if(near_pole(w0)||near_pole(w1)||near_pole(w2)||
		 near_pole(z0)||near_pole(z1)||near_pole(z2)) return;

	double seeds[8][2]; int nseeds=0;
	{
		cx F0=F(w0,z0), F1=F(w1,z1), F2=F(w2,z2);
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
		if(!newton_on_triangle(fr,sx,sy,ol1,ol2)) continue;
		bool dup=false;
		for(int r=0;r<d->nroots;++r)
			if(fabs(d->l1[r]-ol1)<1e-7 && fabs(d->l2[r]-ol2)<1e-7) dup=true;
		if(dup) continue;
		double ol0=1.0-ol1-ol2;
		int r=d->nroots;
		d->l1[r]=ol1; d->l2[r]=ol2;
		d->w_pt[r]=ol0*w0+ol1*w1+ol2*w2;
		d->z_pt[r]=ol0*z0+ol1*z1+ol2*z2;
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
		array<Vertex(T),3> vs; vertices(t,tri,vs);
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

// Radial displacement for the "onion" projection: a fixed-angle
// projection of w's own stereographic image onto a generic direction in
// the (x,y) plane of its sphere, bounded to [-1,1] (points on a unit
// sphere have x^2+y^2<=1). Deliberately not just Re(w)/... along the raw
// axis: it needs to separate the sheets found at a given z, and every
// degree-2 catalog curve so far is of the form w^2=f_0(z), so its two
// roots are always a +-w pair -- any single linear probe of w already
// separates that (odd under negation), but picking a non-axis-aligned
// angle keeps it from being blind to some other curve's symmetry instead
// (e.g. one invariant under w -> -w but not under reflection through the
// real axis).
const double ONION_ANGLE=0.83; // an arbitrary non-special angle (radians)
double onion_radial(cx w) {
	vec<3,double> s=to_sphere(w);
	return std::cos(ONION_ANGLE)*s[0]+std::sin(ONION_ANGLE)*s[1]; // in [-1,1]
}

bool g_onion=false;
double g_onion_scale=0.3;

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
		double r=1.0+g_onion_scale*onion_radial(w);
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
	int nverts, nfaces;
	obj_writer(const char* path) : out(path), nverts(0), nfaces(0) {
		if(g_onion)
			out<<"# Riemann surface extraction, onion projection: "
				<<"direction=z on S^2, radius=1+"<<g_onion_scale<<"*onion_radial(w)\n";
		else
			out<<"# Riemann surface extraction, projected via (Re w, Im w, Re z)\n";
	}
	void write_polygon(const vector<crossing_node>& nodes, const vector<int>& cyc) {
		for(size_t i=0;i<cyc.size();++i) {
			vec<3,double> p=project_for_viz(nodes[cyc[i]].w, nodes[cyc[i]].z);
			out<<"v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<"\n";
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
			GA=std::polar(1.0,0.7); GB=cx(0.13,0.29);
			GC=std::polar(1.0,1.1); GD=cx(-0.21,0.17);
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
		} else if(arg=="--list-functions") {
			cout<<"available functions:"<<endl;
			print_function_catalog(cout);
			return 0;
		} else {
			cerr<<"unrecognized argument: "<<arg<<endl;
			cerr<<"usage: "<<argv[0]<<" [--function N] [--generic] [--depth N] "
					<<"[--threshold X] [--onion] [--onion-scale X] [--list-functions]"<<endl;
			return 1;
		}
	}

	vector<catalog_entry>& catalog=function_catalog();
	if(function_index<0 || function_index>=(int)catalog.size()) {
		cerr<<"--function "<<function_index<<" out of range; available:"<<endl;
		print_function_catalog(cerr);
		return 1;
	}
	g_F=catalog[function_index].F;
	cout<<"function: "<<function_index<<" ("<<catalog[function_index].name<<") -- "
			<<catalog[function_index].description<<endl;
	cout<<"coordinates: "<<(g_generic_coords?"generic (rotated+translated)":"aligned (original)")<<endl;
	cout<<"max_depth="<<max_depth<<" threshold="<<threshold<<endl;
	cout<<"projection: "<<(g_onion?"onion":"flat");
	if(g_onion) cout<<" (scale="<<g_onion_scale<<")";
	cout<<endl;

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
			w_sphere_set(t,v,label_sphere(wl));
			z_sphere_set(t,v,label_sphere(zl));
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
					array<Vertex(T),5> vs;
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
	cout<<"wrote "<<obj_path<<": "<<obj.nverts<<" vertices, "<<obj.nfaces<<" faces "
			<<"(projection: "<<(g_onion?"onion, z on S^2 + radial w":"Re(w), Im(w), Re(z)")<<")"<<endl;
	cout<<"|F|^2 residual at extracted nodes, range: ["<<fres_min<<", "<<fres_max<<"]"<<endl;

	cout<<"ok/bad cells by cell level:"<<endl;
	for(int l=0;l<32;++l) if(ok_by_level[l]||bad_by_level[l])
		cout<<"  level "<<l<<": ok="<<ok_by_level[l]<<" bad="<<bad_by_level[l]<<endl;

	return 0;
}
