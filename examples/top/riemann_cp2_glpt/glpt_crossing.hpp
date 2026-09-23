#ifndef GLPT_CROSSING_HPP
#define GLPT_CROSSING_HPP

/*! \file
 * \brief Curve-crossing detection on a glpt cell's 2-faces (triangles):
 * the second missing piece for a glpt_tree-based CP^2 mesh (see
 * glpt_points.hpp's own header comment for the first -- real vertex
 * coordinates -- and the wider conversation this continues).
 *
 * ============================================================
 * SCOPE: WHAT THIS PORTS FROM riemann_cp2.cpp, AND WHAT IT DEFERS
 * ============================================================
 * riemann_cp2.cpp's triangle_intersection() does three genuinely
 * separate things: (1) pure numerics on 3 given points and a curve F
 * (pick a well-conditioned chart, dehomogenize, run Newton from a
 * handful of seeds, deduplicate) -- entirely independent of HOW those 3
 * points were obtained; (2) a certified Bernstein-Bezier fallback for
 * the rare case every Newton seed misses; (3) nmt<4>-specific caching
 * (extra_data<2>::computed) so a 2-face shared by many 4-cells is only
 * solved once.
 *
 * This file ports (1) essentially verbatim (pick_chart, dehomogenize/
 * rehomogenize, the 3 charts' F/Fa/Fb dispatch, tri_frame/
 * build_tri_frame/bary_to_xy, newton_on_triangle, and
 * triangle_intersection's own seed-generation strategy -- a dynamic
 * linear-system seed plus 4 fixed ones), REPLACES (3) with a cache
 * keyed by the SORTED GLOBAL VERTEX ID TRIPLE of the face (using
 * glpt_vertex_ids.hpp's ids, not an nmt handle) -- which is exactly
 * the consistency mechanism discussed and tested earlier for this same
 * reason (two different 4-cells sharing a face must compute the SAME
 * crossing, not independently re-solve it and risk disagreeing) -- and
 * DEFERS (2): if every Newton seed misses, this reports 0 roots for now,
 * same as riemann_cp2.cpp's own behavior before the Bernstein fallback
 * was added (see project memory riemann_pc2_bernstein_refinement) --
 * matching an already-working, previously-shipped state of that
 * project, not a new gap. Can be added later the same way it was there.
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
#include <map>
#include <cmath>

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
inline int pick_chart(const pt3 p[3]) {
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

//! Core numerics, ported from triangle_intersection() minus the
//! Bernstein fallback (see this file's header) and minus any caching
//! (that's face_crossing_cache's job, below): given 3 points, finds up
//! to 2 curve crossings via Newton from a dynamic + 4 fixed seeds,
//! deduplicated. Returns how many were found (0, 1, or 2); out_pts and
//! out_bary (barycentric relative to p[0],p[1],p[2] in the GIVEN order)
//! are filled for that many entries.
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
	return nfound;
}

// --- Face cache, keyed by the sorted GLOBAL vertex id triple (see this
// file's header comment for why: the same physical 2-face is reached
// from potentially many different 4-cells, and they must all agree).
struct face_key {
	int a,b,c; // sorted ascending
	face_key(int x,int y,int z) {
		int v[3]={x,y,z};
		for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) if(v[j]<v[i]) { int t=v[i]; v[i]=v[j]; v[j]=t; }
		a=v[0]; b=v[1]; c=v[2];
	}
	bool operator<(const face_key& o) const {
		if(a!=o.a) return a<o.a;
		if(b!=o.b) return b<o.b;
		return c<o.c;
	}
};
struct face_result {
	int nroots;
	pt3 pts[2];
};

class face_crossing_cache {
	public:
		//! Looks up or computes the crossing(s) of face (p[0],p[1],p[2])
		//! (global ids id[0..2], any order -- canonicalized internally).
		//! Returns a reference valid until the next insertion.
		const face_result& get(const pt3 p[3], const int id[3]) {
			face_key key(id[0],id[1],id[2]);
			std::map<face_key,face_result>::iterator it = cache_.find(key);
			if(it!=cache_.end()) return it->second;
			face_result fr;
			pt3 out_pts[2]; double out_bary[2][3];
			fr.nroots = compute_face_crossing(p, out_pts, out_bary);
			for(int r=0;r<fr.nroots;++r) fr.pts[r]=out_pts[r];
			return cache_.insert(std::make_pair(key,fr)).first->second;
		}
		size_t size() const { return cache_.size(); }
	private:
		std::map<face_key,face_result> cache_;
};

// --- Enumerates a glpt 4-cell's ten 2-faces (choose 3 of its 5 local
// vertices) as (local vertex index triple), for the caller to extract
// points/ids and query face_crossing_cache.
static const int GLPT_CELL_FACES[10][3] = {
	{0,1,2},{0,1,3},{0,1,4},{0,2,3},{0,2,4},
	{0,3,4},{1,2,3},{1,2,4},{1,3,4},{2,3,4}
};

#endif // GLPT_CROSSING_HPP
