#ifndef GLPT_CERTIFY_HPP
#define GLPT_CERTIFY_HPP

/*! \file
 * \brief Prototype (--certify): certified per-cell tests for a smooth
 * curve F=0 in CP^2, measured against the heuristic extraction.
 *
 * A cell is treated as the FLAT 4-simplex spanned by its 5 vertices in
 * the cell's own best-conditioned affine chart (max over charts of the
 * min |coordinate| over the vertices). With P_i the vertices rescaled so
 * that the chart coordinate is 1, a homogeneous polynomial Q of degree m
 * restricted to the simplex is H(l) = Q(sum_i l_i P_i), a homogeneous
 * polynomial of degree m in the barycentric coordinates l. Its Bernstein
 * coefficients b_alpha = c_alpha / multinomial(m; alpha) enclose its
 * values: H(simplex) is inside the convex hull of {b_alpha} in C.
 *
 * Two tests, both "0 is not in the convex hull of the complex Bernstein
 * coefficients" (sharper than separate Re/Im boxes):
 *
 *  EMPTY: Q = F. Then F has no zero on the simplex.
 *
 *  GRAPH: Q = sum_k conj(F_k(c)) F_k, k over the chart's two free
 *         coordinates, c the centroid; i.e. Q(x) = <grad F(x), grad F(c)>
 *         (Hermitian), the complex derivative of F along the fixed
 *         direction v = conj(grad F(c)). If its values lie in an open
 *         half-plane, then (i) grad F != 0 on the simplex, (ii) on every
 *         complex line in direction v, F is injective on the (convex)
 *         slice by the Noshiro-Warschawski theorem, so C meets each such
 *         line at most once, and (iii) by the implicit function theorem
 *         C cap simplex is the graph of a holomorphic function over its
 *         projection onto the tangent line v-perp. This is the complex
 *         analogue of Plantinga-Vegter's small-normal-variation test.
 *
 * Optional subdivision (--certify-level L): the simplex is bisected at
 * the midpoint (in the chart) of its longest edge, L times recursively.
 * EMPTY needs every piece excluded (each on its own). GRAPH needs the
 * UNION of all pieces' coefficients in one half-plane, since the
 * injectivity argument needs one direction for the whole cell.
 *
 * Not certified here (next step): the boundary -- exact crossing counts
 * on 2-faces and arcs on facets -- which together with GRAPH would make
 * C cap cell a disc. Note also that the extraction solves each 2-face in
 * the FACE's best chart, which may differ from the cell's chart used here.
 * Floating-point rounding of the coefficients is not controlled.
 *
 * Since geodesic faces became the default (glpt_crossing.hpp), the cell
 * realization used here (chart-flat) no longer matches the extraction's
 * faces; certifying the geodesic realization needs interval arithmetic,
 * since F composed with it is not a polynomial. The prototype's main
 * result was diagnostic: it exposed the chart mixing that geodesic faces
 * fixed.
 */

#include "glpt_crossing.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <random>
#include <vector>

bool g_certify=false;
int g_certify_level=0;

typedef std::array<int,5> exp5;
typedef std::map<exp5,cx> poly5; // homogeneous polynomial in 5 barycentric vars

inline poly5 mul5(const poly5& A, const poly5& B) {
	poly5 R;
	for(const auto& a: A) for(const auto& b: B) {
		exp5 e; for(int i=0;i<5;++i) e[i]=a.first[i]+b.first[i];
		R[e]+=a.second*b.second;
	}
	return R;
}
inline poly5 one5() { poly5 R; R[exp5{0,0,0,0,0}]=cx(1,0); return R; }

//! Q(sum_i l_i P_i) for homogeneous Q in (X,Y,Z), as a poly5.
inline poly5 compose5(const poly_F3& Q, const pt3 P[5]) {
	poly5 L[3]; // linear forms X,Y,Z in l
	for(int j=0;j<3;++j) for(int i=0;i<5;++i) { exp5 e{0,0,0,0,0}; e[i]=1; L[j][e]=P[i][j]; }
	std::vector<poly5> pw[3]; // pw[j][k] = L_j^k
	for(int j=0;j<3;++j) { pw[j].push_back(one5()); for(int k=1;k<=Q.d;++k) pw[j].push_back(mul5(pw[j].back(),L[j])); }
	poly5 R;
	for(size_t t=0;t<Q.t.size();++t) {
		const term3& tm=Q.t[t];
		poly5 m=mul5(mul5(pw[0][tm.e[0]],pw[1][tm.e[1]]),pw[2][tm.e[2]]);
		for(const auto& kv: m) R[kv.first]+=tm.c*kv.second;
	}
	return R;
}

inline double fact5(int n) { double r=1; for(int i=2;i<=n;++i) r*=i; return r; }

//! Bernstein coefficients of a homogeneous poly5 of degree m: every
//! multi-index of total degree m (zero monomials included -- their
//! Bernstein coefficient is 0 and must be part of the hull).
inline void bernstein5(const poly5& H, int m, std::vector<cx>& out) {
	for(int a=0;a<=m;++a) for(int b=0;a+b<=m;++b) for(int c=0;a+b+c<=m;++c) for(int d=0;a+b+c+d<=m;++d) {
		int e=m-a-b-c-d;
		exp5 k{a,b,c,d,e};
		poly5::const_iterator it=H.find(k);
		cx v = (it==H.end()) ? cx(0,0) : it->second;
		double mult=fact5(m)/(fact5(a)*fact5(b)*fact5(c)*fact5(d)*fact5(e));
		out.push_back(v/mult);
	}
}

//! Is 0 outside the convex hull of S (a finite set in C)? True iff all
//! points are nonzero and lie in an open half-plane through 0, i.e. the
//! largest angular gap between consecutive arguments exceeds pi.
//! `margin` returns gap-pi (radians), a robustness measure.
inline bool zero_outside_hull(const std::vector<cx>& S, double& margin) {
	margin=-M_PI;
	double maxabs=0; for(const cx& s: S) maxabs=std::max(maxabs,std::abs(s));
	if(maxabs==0) return false;
	std::vector<double> ang;
	for(const cx& s: S) {
		if(std::abs(s)<=1e-14*maxabs) return false; // (numerically) zero coefficient
		ang.push_back(std::arg(s));
	}
	std::sort(ang.begin(),ang.end());
	double gap=ang.front()+2*M_PI-ang.back();
	for(size_t i=1;i<ang.size();++i) gap=std::max(gap,ang[i]-ang[i-1]);
	margin=gap-M_PI;
	return margin>1e-12;
}

inline int pick_chart5(const pt3 p[5]) {
	int chart=0; double best=-1;
	for(int c=0;c<3;++c) {
		double m=1e300; for(int i=0;i<5;++i) m=std::min(m,std::abs(p[i][c]));
		if(m>best) { best=m; chart=c; }
	}
	return chart;
}

inline void rescale_to_chart(const pt3 p[5], int chart, pt3 P[5]) {
	for(int i=0;i<5;++i) for(int j=0;j<3;++j) P[i][j]=p[i][j]/p[i][chart];
}

// Recursive longest-edge bisection (in chart coordinates).
inline void split_longest(const pt3 P[5], pt3 A[5], pt3 B[5]) {
	int bi=0,bj=1; double bd=-1;
	for(int i=0;i<5;++i) for(int j=i+1;j<5;++j) {
		double d=0; for(int k=0;k<3;++k) d+=std::norm(P[i][k]-P[j][k]);
		if(d>bd) { bd=d; bi=i; bj=j; }
	}
	pt3 mid; for(int k=0;k<3;++k) mid[k]=0.5*(P[bi][k]+P[bj][k]);
	for(int i=0;i<5;++i) { A[i]=P[i]; B[i]=P[i]; }
	A[bj]=mid; B[bi]=mid;
}

inline bool empty_rec(const pt3 P[5], int level) {
	std::vector<cx> S; bernstein5(compose5(g_F,P), g_F.d, S);
	double mg;
	if(zero_outside_hull(S,mg)) return true;
	if(level<=0) return false;
	pt3 A[5],B[5]; split_longest(P,A,B);
	return empty_rec(A,level-1) && empty_rec(B,level-1);
}

inline void graph_collect(const pt3 P[5], const poly_F3& G, int level, std::vector<cx>& S) {
	if(level<=0) { bernstein5(compose5(G,P), G.d, S); return; }
	pt3 A[5],B[5]; split_longest(P,A,B);
	graph_collect(A,G,level-1,S); graph_collect(B,G,level-1,S);
}

//! Q = sum_k conj(F_k(c)) F_k over the chart's two free coordinates.
inline poly_F3 graph_poly(int chart, const pt3& c) {
	const poly_F3* d[3]={&g_Fx,&g_Fy,&g_Fz};
	poly_F3 Q; Q.d=(g_F.d>0)?g_F.d-1:0;
	for(int k=0;k<3;++k) {
		if(k==chart) continue;
		cx w=std::conj(eval_poly3(*d[k],c[0],c[1],c[2]));
		for(size_t t=0;t<d[k]->t.size();++t) { term3 tm=d[k]->t[t]; tm.c*=w; Q.t.push_back(tm); }
	}
	return Q;
}

enum cert_result { CERT_EMPTY=0, CERT_GRAPH=1, CERT_UNDECIDED=2 };

//! Classify one cell. `graph_margin` is gap-pi of the GRAPH test (or -pi).
inline cert_result certify_cell(const pt3 pts[5], int level, double& graph_margin) {
	graph_margin=-M_PI;
	int chart=pick_chart5(pts);
	pt3 P[5]; rescale_to_chart(pts,chart,P);
	if(empty_rec(P,level)) return CERT_EMPTY;
	pt3 c; for(int k=0;k<3;++k) { c[k]=0; for(int i=0;i<5;++i) c[k]+=P[i][k]; c[k]/=5.0; }
	poly_F3 Q=graph_poly(chart,c);
	std::vector<cx> S; graph_collect(P,Q,level,S);
	if(zero_outside_hull(S,graph_margin)) return CERT_GRAPH;
	return CERT_UNDECIDED;
}

//! Self-test: on random simplices, sampled values of H must lie inside
//! the coefficient hull (checked through its bounding box, a necessary
//! condition), and the hull test must agree with brute force on easy cases.
inline bool certify_selftest() {
	std::mt19937 rng(7); std::uniform_real_distribution<double> U(-1,1);
	bool ok=true; int nchecks=0;
	for(int trial=0; trial<200; ++trial) {
		pt3 P[5];
		for(int i=0;i<5;++i) { P[i][0]=cx(1,0); P[i][1]=cx(U(rng),U(rng)); P[i][2]=cx(U(rng),U(rng)); }
		std::vector<cx> S; bernstein5(compose5(g_F,P), g_F.d, S);
		double rlo=1e300,rhi=-1e300,ilo=1e300,ihi=-1e300;
		for(const cx& s: S) { rlo=std::min(rlo,s.real()); rhi=std::max(rhi,s.real()); ilo=std::min(ilo,s.imag()); ihi=std::max(ihi,s.imag()); }
		for(int smp=0; smp<50; ++smp) {
			double l[5], sum=0; for(int i=0;i<5;++i) { l[i]=-std::log(std::max(1e-12,0.5*(U(rng)+1))); sum+=l[i]; }
			pt3 x; for(int k=0;k<3;++k) { x[k]=0; for(int i=0;i<5;++i) x[k]+=(l[i]/sum)*P[i][k]; }
			cx v=eval_poly3(g_F,x[0],x[1],x[2]);
			double tol=1e-9*(1+std::abs(v));
			if(v.real()<rlo-tol||v.real()>rhi+tol||v.imag()<ilo-tol||v.imag()>ihi+tol) { ok=false; }
			++nchecks;
		}
		// vertices: H(e_i) = F(P_i) must equal the corner coefficient exactly
		for(int i=0;i<5;++i) {
			exp5 e{0,0,0,0,0}; e[i]=g_F.d;
			poly5 H=compose5(g_F,P); cx hv=H.count(e)?H[e]:cx(0,0);
			cx fv=eval_poly3(g_F,P[i][0],P[i][1],P[i][2]);
			if(std::abs(hv-fv)>1e-9*(1+std::abs(fv))) ok=false;
			++nchecks;
		}
	}
	{ double mg; std::vector<cx> S={cx(1,0),cx(0,1),cx(1,1)}; if(!zero_outside_hull(S,mg)) ok=false; }
	{ double mg; std::vector<cx> S={cx(1,0),cx(-1,0.1),cx(0,-1)}; if(zero_outside_hull(S,mg)) ok=false; }
	std::cout<<"certify selftest: "<<nchecks<<" checks, "<<(ok?"OK":"FAILED")<<std::endl;
	return ok;
}

//! Diagnostic: min barycentric coordinate of the point p with respect to
//! the cell's flat simplex in the cell's chart (>= 0 iff p is inside).
inline double min_bary_in_cell(const pt3 pts[5], const pt3& p) {
	int chart=pick_chart5(pts);
	pt3 P[5]; rescale_to_chart(pts,chart,P);
	int f0=(chart+1)%3, f1=(chart+2)%3;
	cx q0=p[f0]/p[chart], q1=p[f1]/p[chart];
	// unknowns l1..l4 (l0 = 1 - sum): real 4x4 system from Re/Im of 2 complex eqs
	double A[4][5];
	cx r0=q0-P[0][f0], r1=q1-P[0][f1];
	for(int j=0;j<4;++j) {
		cx d0=P[j+1][f0]-P[0][f0], d1=P[j+1][f1]-P[0][f1];
		A[0][j]=d0.real(); A[1][j]=d0.imag(); A[2][j]=d1.real(); A[3][j]=d1.imag();
	}
	A[0][4]=r0.real(); A[1][4]=r0.imag(); A[2][4]=r1.real(); A[3][4]=r1.imag();
	for(int c=0;c<4;++c) { // Gaussian elimination, partial pivoting
		int piv=c; for(int r=c+1;r<4;++r) if(std::fabs(A[r][c])>std::fabs(A[piv][c])) piv=r;
		for(int k=0;k<5;++k) std::swap(A[c][k],A[piv][k]);
		if(std::fabs(A[c][c])<1e-300) return -1e300;
		for(int r=0;r<4;++r) if(r!=c) { double f=A[r][c]/A[c][c]; for(int k=c;k<5;++k) A[r][k]-=f*A[c][k]; }
	}
	double l[5]; l[0]=1;
	for(int j=0;j<4;++j) { l[j+1]=A[j][4]/A[j][j]; l[0]-=l[j+1]; }
	double m=l[0]; for(int j=1;j<5;++j) m=std::min(m,l[j]);
	return m;
}

#endif // GLPT_CERTIFY_HPP
