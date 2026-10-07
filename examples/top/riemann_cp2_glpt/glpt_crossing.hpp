#ifndef GLPT_CROSSING_HPP
#define GLPT_CROSSING_HPP

/*! \file
 * \brief Curve-crossing detection on a glpt cell's 2-faces (triangles):
 * the second missing piece for a glpt_tree-based CP^2 mesh (see
 * glpt_points.hpp's own header comment for the first -- real vertex
 * coordinates -- and the wider conversation this continues).
 *
 * ============================================================
 * SCOPE: WHAT THIS PORTS FROM riemann_cp2.cpp
 * ============================================================
 * riemann_cp2.cpp's triangle_intersection() does three genuinely
 * separate things: (1) pure numerics on 3 given points and a curve F
 * (pick a well-conditioned chart, dehomogenize, run Newton from a
 * handful of seeds, deduplicate); (2) a certified Bernstein-Bezier
 * fallback for the rare case every Newton seed misses; (3) nmt<4>-
 * specific caching (extra_data<2>::computed) so a 2-face shared by many
 * 4-cells is only solved once.
 *
 * This file ports (1) essentially verbatim (pick_chart, dehomogenize/
 * rehomogenize, the 3 charts' F/Fa/Fb dispatch, tri_frame/
 * build_tri_frame/bary_to_xy, newton_on_triangle, and
 * triangle_intersection's own seed-generation strategy -- a dynamic
 * linear-system seed plus 4 fixed ones), and REPLACES (3) with a cache
 * keyed by the SORTED GLOBAL VERTEX ID TRIPLE of the face (using
 * glpt_vertex_ids.hpp's ids, not an nmt handle) -- which is exactly
 * the consistency mechanism discussed and tested earlier for this same
 * reason (two different 4-cells sharing a face must compute the SAME
 * crossing, not independently re-solve it and risk disagreeing). (2),
 * the certified Bernstein-Bezier fallback, is ADDED by glpt_bernstein.hpp
 * (--bernstein), via two small function-pointer hooks declared just
 * above compute_face_crossing() below (g_bernstein_maybe_zero,
 * g_bernstein_fallback) rather than included directly here -- see that
 * file's own header comment for why (avoiding a circular dependency).
 * Both are null by default, so --bernstein off costs one branch.
 *
 * pick_chart, dehomogenize/rehomogenize, the F_chart/Fa_chart/Fb_chart
 * family, real_dot, tri_frame/build_tri_frame/bary_to_xy, and
 * newton_on_triangle below are copied verbatim from riemann_cp2.cpp (a
 * monolithic .cpp, not an includable header) -- the same precedent
 * glpt.hpp and glpt_points.hpp already set.
 */

#include "glpt_points.hpp"
#include "glpt_vertex_ids.hpp"
#include "../riemann_cp2/functions_cp2.hpp"
#include <cmath>
#include <vector>

// --- Current curve (set_curve mirrors riemann_cp2.cpp's own global-F
// design: F_chart*/Fa_chart*/Fb_chart* below are plain function
// pointers, cxfun2, so the curve has to live in globals for them to
// close over it without changing that signature).
poly_F3 g_F, g_Fx, g_Fy, g_Fz;
inline void set_curve(const poly_F3& F) {
	g_F=F; g_Fx=diff_poly3(F,0); g_Fy=diff_poly3(F,1); g_Fz=diff_poly3(F,2);
}

inline cx F_chart0(cx a, cx b) { return eval_poly3(g_F, cx(1,0), a, b); }
inline cx Fa_chart0(cx a, cx b){ return eval_poly3(g_Fy, cx(1,0), a, b); }
inline cx Fb_chart0(cx a, cx b){ return eval_poly3(g_Fz, cx(1,0), a, b); }
inline cx F_chart1(cx a, cx b) { return eval_poly3(g_F, a, cx(1,0), b); }
inline cx Fa_chart1(cx a, cx b){ return eval_poly3(g_Fx, a, cx(1,0), b); }
inline cx Fb_chart1(cx a, cx b){ return eval_poly3(g_Fz, a, cx(1,0), b); }
inline cx F_chart2(cx a, cx b) { return eval_poly3(g_F, a, b, cx(1,0)); }
inline cx Fa_chart2(cx a, cx b){ return eval_poly3(g_Fx, a, b, cx(1,0)); }
inline cx Fb_chart2(cx a, cx b){ return eval_poly3(g_Fy, a, b, cx(1,0)); }

inline void dehomogenize(int chart, const pt3& p, cx& a, cx& b) {
	switch(chart) {
		case 0: a=p[1]/p[0]; b=p[2]/p[0]; break;
		case 1: a=p[0]/p[1]; b=p[2]/p[1]; break;
		default: a=p[0]/p[2]; b=p[1]/p[2]; break;
	}
}
inline pt3 rehomogenize(int chart, cx a, cx b) {
	pt3 p;
	switch(chart) {
		case 0: p[0]=cx(1,0); p[1]=a; p[2]=b; break;
		case 1: p[0]=a; p[1]=cx(1,0); p[2]=b; break;
		default: p[0]=a; p[1]=b; p[2]=cx(1,0); break;
	}
	return normalize3(p);
}

//! Best-conditioned chart for 3 given points: max over charts of the
//! min |coordinate| across the 3 vertices.
int g_forced_face_chart = -1; // >=0: every face in this chart (the affine-chart baseline, --box)
inline int pick_chart(const pt3 p[3]) {
	if(g_forced_face_chart>=0) return g_forced_face_chart;
	int chart=0; double bestscore=-1;
	for(int c=0;c<3;++c) {
		double m=std::min(std::abs(p[0][c]), std::min(std::abs(p[1][c]),std::abs(p[2][c])));
		if(m>bestscore) { bestscore=m; chart=c; }
	}
	return chart;
}

inline double real_dot(cx a, cx b) { return a.real()*b.real()+a.imag()*b.imag(); }

struct tri_frame {
	int io, ip, iq;
	cx Pow, Poz;
	cx Xw, Xz, Yw, Yz;
	double xp;
};
inline void build_tri_frame(const cx w[3], const cx z[3], tri_frame& fr) {
	double d01=std::norm(w[1]-w[0])+std::norm(z[1]-z[0]);
	double d12=std::norm(w[2]-w[1])+std::norm(z[2]-z[1]);
	double d02=std::norm(w[2]-w[0])+std::norm(z[2]-z[0]);
	if(d01>=d12 && d01>=d02)      { fr.io=2; fr.ip=0; fr.iq=1; }
	else if(d12>=d01 && d12>=d02) { fr.io=0; fr.ip=1; fr.iq=2; }
	else                          { fr.io=1; fr.ip=0; fr.iq=2; }
	cx dw=w[fr.iq]-w[fr.ip], dz=z[fr.iq]-z[fr.ip];
	cx ow=w[fr.io]-w[fr.ip], oz=z[fr.io]-z[fr.ip];
	double denom=real_dot(dw,dw)+real_dot(dz,dz);
	double tt=(denom>1e-300) ? (real_dot(ow,dw)+real_dot(oz,dz))/denom : 0.5;
	fr.Pow=w[fr.ip]+tt*dw; fr.Poz=z[fr.ip]+tt*dz;
	fr.Xw=w[fr.iq]-fr.Pow; fr.Xz=z[fr.iq]-fr.Poz;
	fr.Yw=w[fr.io]-fr.Pow; fr.Yz=z[fr.io]-fr.Poz;
	double denom2=1.0-tt;
	fr.xp=(fabs(denom2)>1e-9) ? -tt/denom2 : -1e9;
}
inline void bary_to_xy(const tri_frame& fr, double l1, double l2, double& x, double& y) {
	double l[3]={1.0-l1-l2,l1,l2};
	y=l[fr.io];
	x=l[fr.ip]*fr.xp+l[fr.iq];
}

typedef cx (*cxfun2)(cx,cx);

inline bool newton_on_triangle(const tri_frame& fr, double x, double y, double& out_l1, double& out_l2,
		cxfun2 Ffun, cxfun2 Fwfun, cxfun2 Fzfun, double& out_resid, double& out_geom_resid, int& out_iters) {
	int iter=0;
	for(; iter<20; ++iter) {
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
		if(fabs(x)>1e6||fabs(y)>1e6) return false;
	}
	cx w=fr.Pow+x*fr.Xw+y*fr.Yw;
	cx z=fr.Poz+x*fr.Xz+y*fr.Yz;
	cx Fv=Ffun(w,z);
	if(std::norm(Fv)>1e-20) return false;

	cx fw=Fwfun(w,z), fz=Fzfun(w,z);
	cx cX=fw*fr.Xw+fz*fr.Xz;
	cx cY=fw*fr.Yw+fz*fr.Yz;
	double gnorm=std::sqrt(std::norm(cX)+std::norm(cY));
	out_resid=std::abs(Fv);
	out_geom_resid=(gnorm>1e-300) ? out_resid/gnorm : out_resid;
	out_iters=iter;

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

// --bernstein hooks (see glpt_bernstein.hpp for the real implementation
// and the full rationale): both null by default (--bernstein off, one
// branch of overhead), set by glpt_bernstein.hpp's enable_bernstein()
// only when that mode is actually on. Kept as function pointers rather
// than #including glpt_bernstein.hpp here to avoid a circular
// dependency -- glpt_bernstein.hpp itself needs THIS file's pick_chart/
// dehomogenize/build_tri_frame/bary_to_xy/newton_on_triangle/g_F.
// g_bernstein_maybe_zero: certified "can this face contain a zero at
// all?" test, checked by face_crossing_cache::get() BEFORE running any
// Newton search -- a true "certified empty" skip, not a heuristic.
// g_bernstein_fallback: a certified fallback SEED for Newton (not a
// replacement for it), tried by compute_face_crossing only when every
// fixed/dynamic seed already missed.
typedef bool (*bernstein_maybe_zero_fn)(const pt3 p[3], const int id[3]);
typedef bool (*bernstein_fallback_fn)(const pt3 p[3], int chart, double& out_l1, double& out_l2);
bernstein_maybe_zero_fn g_bernstein_maybe_zero = 0;
bernstein_fallback_fn g_bernstein_fallback = 0;
long g_bernstein_pruned_faces = 0; // diagnostic: faces Newton never had to touch
long g_bfallback_tried = 0, g_bfallback_new_root = 0; // fallback diagnostics

//! Core numerics, ported from triangle_intersection() minus the caching
//! (that's face_crossing_cache's job, below): given 3 points, finds up
//! to 2 curve crossings via Newton from a dynamic + 4 fixed seeds,
//! deduplicated, falling back to g_bernstein_fallback (if set) only
//! when every one of those seeds misses. Returns how many were found
//! (0, 1, or 2); out_pts and out_bary (barycentric relative to
//! p[0],p[1],p[2] in the GIVEN order) are filled for that many entries.
inline int compute_face_crossing(const pt3 p[3], pt3 out_pts[2], double out_bary[2][3]) {
	int chart=pick_chart(p);
	cx w0,w1,w2,z0,z1,z2;
	dehomogenize(chart,p[0],w0,z0);
	dehomogenize(chart,p[1],w1,z1);
	dehomogenize(chart,p[2],w2,z2);
	cxfun2 Ffun,Fwfun,Fzfun;
	switch(chart) {
		case 0: Ffun=F_chart0; Fwfun=Fa_chart0; Fzfun=Fb_chart0; break;
		case 1: Ffun=F_chart1; Fwfun=Fa_chart1; Fzfun=Fb_chart1; break;
		default: Ffun=F_chart2; Fwfun=Fa_chart2; Fzfun=Fb_chart2; break;
	}

	double seeds[5][2]; int nseeds=0;
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
	static const double extra_seeds[4][2]={ {1/3.,1/3.},{0.1,0.1},{0.8,0.1},{0.1,0.8} };
	for(int s=0;s<4;++s) { seeds[nseeds][0]=extra_seeds[s][0]; seeds[nseeds][1]=extra_seeds[s][1]; ++nseeds; }

	tri_frame fr;
	{ cx w3[3]={w0,w1,w2}, z3[3]={z0,z1,z2}; build_tri_frame(w3,z3,fr); }

	int nfound=0;
	for(int s=0; s<nseeds && nfound<2; ++s) {
		double sx,sy; bary_to_xy(fr,seeds[s][0],seeds[s][1],sx,sy);
		double ol1,ol2,rF,rG; int riter;
		if(!newton_on_triangle(fr,sx,sy,ol1,ol2,Ffun,Fwfun,Fzfun,rF,rG,riter)) continue;
		bool dup=false;
		for(int r=0;r<nfound;++r)
			if(fabs(out_bary[r][1]-ol1)<1e-7 && fabs(out_bary[r][2]-ol2)<1e-7) dup=true;
		if(dup) continue;
		double ol0=1.0-ol1-ol2;
		out_bary[nfound][0]=ol0; out_bary[nfound][1]=ol1; out_bary[nfound][2]=ol2;
		cx aroot=ol0*w0+ol1*w1+ol2*w2, broot=ol0*z0+ol1*z1+ol2*z2;
		out_pts[nfound]=rehomogenize(chart,aroot,broot);
		++nfound;
	}

	// Certified fallback (see glpt_bernstein.hpp's own comment on
	// bernstein_locate_root): only tried when every fixed/dynamic seed
	// above missed entirely.
	if(nfound==0 && g_bernstein_fallback) {
		++g_bfallback_tried;
		double bl1,bl2;
		if(g_bernstein_fallback(p, chart, bl1, bl2)) {
			double sx,sy; bary_to_xy(fr,bl1,bl2,sx,sy);
			double ol1,ol2,rF,rG; int riter;
			if(newton_on_triangle(fr,sx,sy,ol1,ol2,Ffun,Fwfun,Fzfun,rF,rG,riter)) {
				++g_bfallback_new_root;
				double ol0=1.0-ol1-ol2;
				out_bary[0][0]=ol0; out_bary[0][1]=ol1; out_bary[0][2]=ol2;
				cx aroot=ol0*w0+ol1*w1+ol2*w2, broot=ol0*z0+ol1*z1+ol2*z2;
				out_pts[0]=rehomogenize(chart,aroot,broot);
				nfound=1;
			}
		}
	}
	return nfound;
}

// --- Geodesic faces (the DEFAULT since 2026-10-01; --chart-flat-faces
// restores the old realization): a chart-free, globally consistent
// realization of every 2-face.
//
// The old realization took each face as the FLAT triangle in the
// face's own best chart. "Flat" depends on the chart, so two faces that
// share an edge but use different charts see that edge as two different
// curves, tetrahedral boundaries do not close, and crossing parity on a
// facet breaks. Measured: bad-cell rate 0.1-0.2% when all faces of a cell
// share the cell's chart, 14-21% otherwise. With geodesic faces the
// overall bad-cell rate at depth 14 fell from 3-12% to ~0.13% on four
// smooth curves, and facets with a single crossing all but vanished
// (closed tetrahedral boundaries meet the closed curve C an even number
// of times).
//
// Here a face with vertices A,B,C, sorted by global id, is the geodesic
// cone from its lowest-id vertex A:
//     q(t)   = (1-t) B + t align(B,C)
//     v(s,t) = (1-s) A + s align(A, q(t)),     (s,t) in [0,1]^2,
// with unit representatives. Every edge is a Fubini-Study geodesic,
// parametrized as [(1-u) X + u align(X,Y)] from its lower-id endpoint X,
// whichever face it is seen from, so faces glue exactly along edges and
// fs_midpoint (which is u=1/2 of that same parametrization) lies on its
// edge. Barycentric-like coordinates (sorted order) are
//     l_A = 1-s,  l_B = s(1-t),  l_C = s t,
// which agree with the edge parametrizations on the boundary.
bool g_geodesic_faces=true; // default since 2026-10-01; --chart-flat-faces restores the old realization

struct geo_face { pt3 A,B,C; };
inline geo_face make_geo_face(const pt3& A, const pt3& B, const pt3& C) {
	geo_face g; g.A=normalize3(A); g.B=normalize3(B); g.C=normalize3(C); return g;
}
inline pt3 geo_point(const geo_face& g, double s, double t) {
	pt3 Ct=align_phase(g.B,g.C), q;
	for(int k=0;k<3;++k) q[k]=(1.0-t)*g.B[k]+t*Ct[k];
	pt3 qa=align_phase(g.A,q), v;
	for(int k=0;k<3;++k) v[k]=(1.0-s)*g.A[k]+s*qa[k];
	return v;
}
//! v(s,t) and its two partial derivatives, analytically. With
//! z = <A,q>, u = z/|z| (the alignment phase) and qa = q u:
//!   dv/ds = qa - A,
//!   dv/dt = s (q' u + q u'),  q' = align(B,C) - B,  z' = <A,q'>,
//!   u' = (z' - u Re(conj(u) z')) / |z|.
inline void geo_point_derivs(const geo_face& g, double s, double t, pt3& v, pt3& vs, pt3& vt) {
	pt3 Ct=align_phase(g.B,g.C), q, dq;
	for(int k=0;k<3;++k) { q[k]=(1.0-t)*g.B[k]+t*Ct[k]; dq[k]=Ct[k]-g.B[k]; }
	cx z=hdot(g.A,q), dz=hdot(g.A,dq);
	double az=std::abs(z);
	cx u = (az>1e-14) ? z/az : cx(1,0);
	cx du = (az>1e-14) ? (dz - u*std::real(std::conj(u)*dz))/az : cx(0,0);
	for(int k=0;k<3;++k) {
		cx qa=q[k]*u;
		v[k]=(1.0-s)*g.A[k]+s*qa;
		vs[k]=qa-g.A[k];
		vt[k]=s*(dq[k]*u+q[k]*du);
	}
}
//! Complex gradient of the homogeneous F at v, and the bilinear
//! directional derivative grad F(v) . w.
inline void geo_grad(const pt3& v, cx gF[3]) {
	gF[0]=eval_poly3(g_Fx,v[0],v[1],v[2]); gF[1]=eval_poly3(g_Fy,v[0],v[1],v[2]); gF[2]=eval_poly3(g_Fz,v[0],v[1],v[2]);
}
inline cx geo_dir(const cx gF[3], const pt3& w) { return gF[0]*w[0]+gF[1]*w[1]+gF[2]*w[2]; }
//! |F(v)| / |v|^n: the residual of the normalized point.
inline double geo_resid(const pt3& v) {
	return std::abs(eval_poly3(g_F,v[0],v[1],v[2]))/std::pow(hnorm(v),(double)g_F.d);
}
//! Newton on G(s,t)=F(v(s,t)) as a real 2x2 system with the analytic
//! Jacobian (F homogeneous, so the zeros of F(v) and F(v/|v|) agree).
//! Accepts a root inside [0,1]^2 (up to 1e-6).
inline bool geo_newton(const geo_face& g, double& s, double& t) {
	pt3 v,vs,vt; cx gF[3];
	for(int it=0; it<40; ++it) {
		geo_point_derivs(g,s,t,v,vs,vt);
		cx G=eval_poly3(g_F,v[0],v[1],v[2]);
		if(std::abs(G)/std::pow(hnorm(v),(double)g_F.d)<1e-15) break;
		geo_grad(v,gF);
		cx Gs=geo_dir(gF,vs), Gt=geo_dir(gF,vt);
		double j00=Gs.real(), j01=Gt.real(), j10=Gs.imag(), j11=Gt.imag();
		double det=j00*j11-j01*j10;
		if(std::fabs(det)<1e-300) return false;
		s-=( G.real()*j11-j01*G.imag())/det;
		t-=(-G.real()*j10+j00*G.imag())/det;
		if(std::fabs(s)>10||std::fabs(t)>10) return false;
	}
	if(geo_resid(geo_point(g,s,t))>1e-11) return false;
	const double tol=1e-6;
	return s>=-tol && s<=1+tol && t>=-tol && t<=1+tol;
}
//! Sorted order (by id) of a face's three corners: ord[0] = lowest id.
inline void geo_sort(const int id[3], int ord[3]) {
	ord[0]=0; ord[1]=1; ord[2]=2;
	for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) if(id[ord[j]]<id[ord[i]]) std::swap(ord[i],ord[j]);
}
//! Same contract as compute_face_crossing (bary relative to p[] in the
//! GIVEN order), on the geodesic-cone realization.
inline int compute_face_crossing_geo(const pt3 p[3], const int id[3], pt3 out_pts[2], double out_bary[2][3]) {
	int ord[3]; geo_sort(id,ord);
	geo_face g=make_geo_face(p[ord[0]],p[ord[1]],p[ord[2]]);
	// seeds in sorted barycentrics (lA,lB,lC): the zero of the linear
	// interpolant of F at the phase-aligned corners, plus 4 fixed points
	double seeds[5][2]; int nseeds=0;
	{
		pt3 Bt=align_phase(g.A,g.B), Ct=align_phase(g.A,g.C);
		cx F0=eval_poly3(g_F,g.A[0],g.A[1],g.A[2]), F1=eval_poly3(g_F,Bt[0],Bt[1],Bt[2]), F2=eval_poly3(g_F,Ct[0],Ct[1],Ct[2]);
		cx a1=F1-F0, a2=F2-F0;
		double m00=a1.real(), m01=a2.real(), m10=a1.imag(), m11=a2.imag();
		double det=m00*m11-m01*m10;
		if(std::fabs(det)>1e-14) {
			seeds[nseeds][0]=(-F0.real()*m11+m01*F0.imag())/det;
			seeds[nseeds][1]=(-m00*F0.imag()+F0.real()*m10)/det;
			++nseeds;
		}
	}
	static const double extra_seeds[4][2]={ {1/3.,1/3.},{0.1,0.1},{0.8,0.1},{0.1,0.8} };
	for(int k=0;k<4;++k) { seeds[nseeds][0]=extra_seeds[k][0]; seeds[nseeds][1]=extra_seeds[k][1]; ++nseeds; }
	int nfound=0; double found_l[2][3];
	for(int k=0; k<nseeds && nfound<2; ++k) {
		double lB=seeds[k][0], lC=seeds[k][1];
		double s=lB+lC, t=(s>1e-12)? lC/s : 0.5;
		if(!geo_newton(g,s,t)) continue;
		// the root itself may lie up to 1e-6 outside [0,1]^2: clamp only the
		// coordinates used to classify it, not the point (clamping the point
		// moves it off the curve by ~1e-6 |grad F|)
		pt3 root=geo_point(g,s,t);
		s=std::min(1.0,std::max(0.0,s)); t=std::min(1.0,std::max(0.0,t));
		double l[3]={1.0-s, s*(1.0-t), s*t};
		bool dup=false;
		for(int r=0;r<nfound;++r) if(std::fabs(found_l[r][1]-l[1])<1e-7 && std::fabs(found_l[r][2]-l[2])<1e-7) dup=true;
		if(dup) continue;
		for(int q=0;q<3;++q) found_l[nfound][q]=l[q];
		out_pts[nfound]=normalize3(root);
		for(int q=0;q<3;++q) out_bary[nfound][ord[q]]=l[q]; // back to the caller's order
		++nfound;
	}
	return nfound;
}
//! Transversality score on the geodesic realization: the same quantity
//! as crossing_transversality (Jacobian determinant of (Re F, Im F) on the
//! face, normalized to [0,1]), with the face's tangent vectors at the
//! root projected orthogonally to the root's own line.
inline double crossing_transversality_geo(const pt3 face_pts[3], const int face_ids[3], const double bary[3]) {
	int ord[3]; geo_sort(face_ids,ord);
	geo_face g=make_geo_face(face_pts[ord[0]],face_pts[ord[1]],face_pts[ord[2]]);
	double lA=bary[ord[0]], lC=bary[ord[2]];
	double s=1.0-lA; if(s<1e-9) return 1.0; // at the apex: a vertex node, no tangency signal
	double t=std::min(1.0,std::max(0.0,lC/s));
	pt3 v, Ts, Tt; geo_point_derivs(g,s,t,v,Ts,Tt);
	double vv=hdot(v,v).real();
	cx cs=hdot(Ts,v)/vv, ct=hdot(Tt,v)/vv;
	for(int k=0;k<3;++k) { Ts[k]-=cs*v[k]; Tt[k]-=ct*v[k]; }
	cx gF[3]; geo_grad(v,gF);
	cx w1(0,0), w2(0,0); double g2=0;
	for(int k=0;k<3;++k) { w1+=gF[k]*Ts[k]; w2+=gF[k]*Tt[k]; g2+=std::norm(gF[k]); }
	double n1=hnorm(Ts), n2=hnorm(Tt);
	if(n1<1e-14||n2<1e-14) return 1.0;
	if(g2<1e-300) return 0.0;
	return std::abs(std::imag(std::conj(w1)*w2))/(n1*n2*g2);
}

// --- Face cache, keyed by the sorted GLOBAL vertex id triple (see this
// file's header comment for why: the same physical 2-face is reached
// from potentially many different 4-cells, and they must all agree).
//
// face_key: EXACT (no truncation) identity for a sorted vertex-id
// triple -- hi packs (a,b) as 32+32 bits into one uint64_t, lo holds c
// as a plain uint32_t, so a,b,c can each range over the FULL 32-bit id
// space. REVERTED 2026-09-24 from a single-uint64_t 20-bit-per-id
// packing (see glpt_vertex_ids.hpp's own glpt_edge_cache comment for
// the full story: a real mesh minted more than 2^20-1 vertices and hit
// that scheme's own loud assertion). Since face_crossing_cache/
// bernstein_bounds_cache both need an exact, collision-free key (a
// hash collision here would silently return the WRONG face's cached
// crossing data), face_key stays a genuine two-field key, not a single
// lossy hash -- open addressing (below) still avoids std::map's tree
// overhead, it just needs an explicit `present` flag per slot again
// instead of folding it into an unused key bit.
struct face_key { uint64_t hi; uint32_t lo; };
inline face_key make_face_key(int a, int b, int c) {
	int v[3]={a,b,c};
	for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) if(v[j]<v[i]) { int t=v[i]; v[i]=v[j]; v[j]=t; }
	face_key k;
	k.hi = (uint64_t(uint32_t(v[0]))<<32) | uint32_t(v[1]);
	k.lo = uint32_t(v[2]);
	return k;
}
inline bool operator==(const face_key& x, const face_key& y) { return x.hi==y.hi && x.lo==y.lo; }
inline uint64_t glpt_hash_face_key(const face_key& k) {
	return glpt_hash64(glpt_hash64(k.hi) ^ (uint64_t(k.lo)*0x9E3779B97F4A7C15ULL));
}

// COMPACT on purpose: unsigned char + int is 8 bytes (with padding),
// not the ~100 bytes two full pt3 roots would cost inline -- see
// face_crossing_cache's own comment for why this matters (most cached
// faces, in a realistically refined mesh, have nroots==0: nowhere near
// the curve). Mirrors riemann_cp2.cpp's own extra_data<2> (nroots +
// root_idx into a shared, compact g_roots vector), which already solved
// this exact problem for the SAME reason there -- not a new idea, just
// finally applied here too once it was found to matter (see below).
struct face_result {
	unsigned char nroots;
	int root_idx; // index into face_crossing_cache's own root list; meaningless if nroots==0
};

// Open addressing (glpt_tree.hpp's own design in ~/code/lpt, reusing
// its hash/prime helpers via glpt_hash_util.hpp), NOT std::map -- an
// earlier std::map version was found, by directly measuring
// riemann_cp2_glpt.cpp's peak RSS against riemann_cp2.cpp, to be the
// dominant memory cost at any real refinement depth (2026-09-23; see
// glpt_vertex_ids.hpp's glpt_edge_cache for the identical finding and
// fix, applied here the same way): red-black-tree node overhead (3
// pointers + color) on top of the slot data.
//
// SECOND fix found by the SAME measurement, after the first (map ->
// open addressing) unexpectedly made peak RSS WORSE, not better:
// face_result used to store its up-to-2 root points INLINE (2 pt3 = 96
// bytes) in EVERY slot, whether or not that face actually had a root --
// and in a realistically refined mesh MOST cached faces sit nowhere
// near the curve (nroots==0), so most of that 96 bytes/slot was pure
// waste. At ~982000 distinct faces (--depth 13 on the conic curve) that
// was ~94MB on its own, more than riemann_cp2.cpp's ENTIRE peak RSS at
// the same settings (171MB). Root points now live in a separate,
// compact root_pts_ vector, appended to only when a face actually has
// one -- exactly mirroring riemann_cp2.cpp's own extra_data<2>/g_roots
// split, which already solved this for the identical reason there.
//
// THIRD fix, since superseded (see face_key's own comment above): the
// face_key struct (3 plain ints, 12 bytes) plus a separate bool present
// array were packed into a single uint64_t for a while -- reverted back
// to a genuine (if slightly larger) key + present array once that
// packing's own 20-bit id ceiling proved too small for a real mesh.
//
// No deletion needed here (faces are never forgotten), so this is
// simpler than glpt_edge_cache's own version -- insert/lookup/grow only.
class face_crossing_cache {
	public:
		explicit face_crossing_cache(size_t initial_buckets = 1031)
			: keys_(0), present_(0), results_(0), nbuckets_(0), count_(0)
		{
			alloc_(glpt_next_prime(initial_buckets));
		}
		~face_crossing_cache() { std::free(keys_); std::free(present_); std::free(results_); }

		//! Looks up or computes the crossing(s) of face (p[0],p[1],p[2])
		//! (global ids id[0..2], any order -- canonicalized internally).
		//! Returns a reference valid until the next insertion triggers a
		//! grow (rehash reallocates every slot) -- like glpt_edge_cache,
		//! not meant to be held onto across a later get() call. Use
		//! root_point(fr,r)/root_bary(fr,r) to read back a found root's
		//! actual position and its barycentric coordinates relative to
		//! p[0],p[1],p[2] in the order THIS call was made with (needed by
		//! glpt_extraction.hpp to classify a root as sitting on a vertex,
		//! edge, or the face interior).
		const face_result& get(const pt3 p[3], const int id[3]) {
			face_key key = make_face_key(id[0],id[1],id[2]);
			size_t slot = find_slot_(key);
			if(slot!=size_t(-1)) return results_[slot];

			// --bernstein: certified "can't contain a zero" skip -- Newton
			// never has to run at all here (see glpt_bernstein.hpp's own
			// comment on g_bernstein_maybe_zero).
			if(g_bernstein_maybe_zero && !g_bernstein_maybe_zero(p,id)) {
				++g_bernstein_pruned_faces;
				face_result fr; fr.nroots=0; fr.root_idx=-1;
				return insert_(key, fr);
			}

			face_result fr;
			pt3 out_pts[2]; double out_bary[2][3];
			int nr = g_geodesic_faces ? compute_face_crossing_geo(p, id, out_pts, out_bary)
			                          : compute_face_crossing(p, out_pts, out_bary);
			fr.nroots = (unsigned char)nr;
			fr.root_idx = nr>0 ? (int)root_pts_.size() : -1;
			for(int r=0;r<nr;++r) {
				root_pts_.push_back(out_pts[r]);
				bary3 b; b.l[0]=out_bary[r][0]; b.l[1]=out_bary[r][1]; b.l[2]=out_bary[r][2];
				root_bary_.push_back(b);
			}
			return insert_(key, fr);
		}
		pt3 root_point(const face_result& fr, int r) const { return root_pts_[fr.root_idx+r]; }
		const double* root_bary(const face_result& fr, int r) const { return root_bary_[fr.root_idx+r].l; }
		size_t size() const { return count_; }

	private:
		struct bary3 { double l[3]; };
		face_key* keys_;
		bool* present_;
		face_result* results_;
		size_t nbuckets_, count_;
		std::vector<pt3> root_pts_;
		std::vector<bary3> root_bary_;

		void alloc_(size_t n) {
			keys_ = (face_key*)std::calloc(n,sizeof(face_key));
			present_ = (bool*)std::calloc(n,sizeof(bool));
			results_ = (face_result*)std::calloc(n,sizeof(face_result));
			assert(keys_!=0 && present_!=0 && results_!=0 && "face_crossing_cache: out of memory");
			nbuckets_ = n;
		}
		size_t find_slot_(const face_key& key) const {
			size_t h = glpt_hash_face_key(key) % nbuckets_;
			while(present_[h]) {
				if(keys_[h]==key) return h;
				h=(h+1)%nbuckets_;
			}
			return size_t(-1);
		}
		void grow_if_needed_() {
			if(double(count_+1) <= 0.7*double(nbuckets_)) return;
			face_key* old_k=keys_; bool* old_p=present_; face_result* old_r=results_;
			size_t old_n=nbuckets_;
			alloc_(glpt_next_prime(2*old_n));
			count_=0;
			for(size_t i=0;i<old_n;++i) if(old_p[i]) insert_(old_k[i], old_r[i]);
			std::free(old_k); std::free(old_p); std::free(old_r);
		}
		const face_result& insert_(const face_key& key, const face_result& val) {
			grow_if_needed_();
			size_t h = glpt_hash_face_key(key) % nbuckets_;
			while(present_[h]) h=(h+1)%nbuckets_;
			keys_[h]=key; present_[h]=true; results_[h]=val;
			++count_;
			return results_[h];
		}

		face_crossing_cache(const face_crossing_cache&);
		face_crossing_cache& operator=(const face_crossing_cache&);
};

// --- Enumerates a glpt 4-cell's ten 2-faces (choose 3 of its 5 local
// vertices) as (local vertex index triple), for the caller to extract
// points/ids and query face_crossing_cache.
static const int GLPT_CELL_FACES[10][3] = {
	{0,1,2},{0,1,3},{0,1,4},{0,2,3},{0,2,4},
	{0,3,4},{1,2,3},{1,2,4},{1,3,4},{2,3,4}
};

#endif // GLPT_CROSSING_HPP
