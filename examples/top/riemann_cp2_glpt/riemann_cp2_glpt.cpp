// riemann_cp2_glpt.cpp -- glpt_tree-based rebuild of riemann_cp2.cpp's
// pipeline (see glpt_points.hpp/glpt_crossing.hpp/glpt_extraction.hpp's
// own header comments for the three pieces this builds on, and the
// conversation that led to this file for the full rationale: an
// O(1)-per-current-cell mesh, ~/code/lpt/glpt_tree.hpp, replacing
// riemann_cp2.cpp's append-only nmt<4>-backed one).
//
// Seeds all 108 Gaifullin cells, then refines by CONTINUATION (default):
// find one root the curve certainly touches (extract_cell()'s own
// any_edge test, not a proxy), then explore outward through only the
// touched cells via glpt_tree::compat_bisect()/recent_leaves() (down to
// --depth) and glpt_tree::neighbor_leaf() (across cells already at
// --depth) -- see Phase 1's own comment below and
// [[riemann_glpt_continuation]] for why this reaches the exact same
// cells as a global scan, for a fraction of the mesh size, on a
// CONNECTED curve. --legacy-priority-refine restores the ORIGINAL
// global gradient-dispersion cell_priority priority-queue refinement
// (ported from riemann_cp2.cpp), kept for comparison. Either way,
// extracts actual connected polygons per cell (glpt_extraction.hpp's
// extract_cell(), a faithful port of riemann_cp2.cpp's own Phase-3
// per-cell graph-building + cycle-decomposition loop) and writes them
// to an OBJ mesh under one of two projections, matching riemann_cp2.cpp's
// own --flat (default)/--alpha-projection choice and semantics exactly:
// flat (per-polygon best-conditioned affine chart, (Re a, Im a, Re b),
// clippable via --cutoff) or alpha-tilde (Dutter, arXiv:2608.04323,
// after Kranich 2015 -- chart-free, bounded by construction, no
// --cutoff needed). Also ports --generic/--generic-seed (random unitary
// change of basis on C^3, applied once to the 15 Gaifullin seed points
// -- see random_unitary()'s own comment below), --proximity (extra
// diam/(mind+diam) factor on cell_priority, biasing refinement toward
// the curve itself -- see g_proximity's own comment below; only has an
// effect under --legacy-priority-refine, since continuation is already
// gated on the curve's actual presence), and
// --bernstein/--bernstein-level/--bernstein-selftest (certified
// Bernstein-Bezier enclosure per 2-face, both pruning cells that
// provably can't contain a zero and seeding Newton when every ordinary
// seed misses -- see glpt_bernstein.hpp's own header comment). Also
// has an experimental --repair-rounds/--repair-extra-depth (not in
// riemann_cp2.cpp): force-refines specifically the cells extraction
// failed on, then re-extracts -- see that flag's own comment below.
// --tangency-threshold/--tangency-extra-depth (also not in
// riemann_cp2.cpp) is the PROACTIVE counterpart: glpt_extraction.hpp's
// crossing_transversality() score, measured (2026-09-27) to predict
// bad-cell risk with a strong, mostly monotonic gradient (22% bad rate
// for score<=0.05, down to ~3% for score in (0.6,0.8]), is checked on
// every cell continuation would otherwise finalize at --depth; a
// suspiciously near-tangent one gets bisected further first, up to
// --tangency-extra-depth beyond --depth, instead of waiting for
// extraction to actually fail before --repair-rounds can act on it.
// FACES (2026-10-01): every 2-face is realized as the Fubini-Study
// geodesic cone from its lowest-id vertex (glpt_crossing.hpp, "Geodesic
// faces"), so faces glue exactly along geodesic edges; this cut the
// bad-cell rate by 20-80x. --chart-flat-faces restores the old
// realization (each face flat in its own best chart), which --bernstein
// still requires. --certify (glpt_certify.hpp) is a measurement
// prototype of certified per-cell tests.
// Deliberately NOT ported: riemann_cp2.cpp's --onion mode (flat+alpha
// are enough to see the extracted surface; onion is deferred, not
// required).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <queue>
#include <fstream>
#include <iostream>
#include <random>
#include <array>
#include <chrono>
#include "glpt_extraction.hpp"
#include "glpt_bernstein.hpp"
#include "glpt_certify.hpp"
#include "glpt_tree.hpp"

using namespace std;

// --alpha-projection (see this file's own header comment): a function
// of the homogeneous point itself, no chart/branch choice at all.
// Copied verbatim from riemann_cp2.cpp's own alpha_projection().
void alpha_projection(const pt3& p, double out[3]) {
	double S=std::norm(p[0])+std::norm(p[1])+std::norm(p[2]);
	cx xy=std::conj(p[0])*p[1];
	out[0]=std::norm(p[0])/S;
	out[1]=xy.real()/S;
	out[2]=xy.imag()/S;
}

// --flat (default): per-polygon best-conditioned affine chart, copied
// verbatim from riemann_cp2.cpp's own pick_chart_polygon() -- picks
// ONE chart for the whole polygon (maximizing the minimum |coordinate|
// across all of its nodes), not per point; that file's own comment
// documents a real bug (0.69% of polygons, "ribbon spray" artifact)
// from picking per-point instead, which this avoids by construction.
// --chart N (flat mode only): force affine chart N (0: X=1, 1: Y=1,
// 2: Z=1) for every polygon instead of the per-polygon best one --
// for figures of a curve as a graph over one chart (e.g. w^2=z^3-z in
// Z=1); combine with --cutoff to drop polygons escaping to infinity.
int g_forced_chart=-1;
// --flat-swap: output (Re b, Im b, Re a) instead of (Re a, Im a, Re b).
bool g_flat_swap=false;
int pick_chart_polygon(const std::vector<crossing_node>& nodes, const std::vector<int>& cyc) {
	if(g_forced_chart>=0) return g_forced_chart;
	int chart=0; double best=-1;
	for(int c=0;c<3;++c) {
		double m=1e300;
		for(size_t i=0;i<cyc.size();++i) {
			double v=std::abs(nodes[cyc[i]].p[c]);
			if(v<m) m=v;
		}
		if(m>best) { best=m; chart=c; }
	}
	return chart;
}
// Same (Re a, Im a, Re b) affine-chart projection as riemann_cp2.cpp's
// own project_for_viz() flat branch, using glpt_crossing.hpp's own
// dehomogenize().
void flat_projection(const pt3& p, int chart, double out[3]) {
	cx a,b; dehomogenize(chart,p,a,b);
	if(g_flat_swap) std::swap(a,b);
	out[0]=a.real(); out[1]=a.imag(); out[2]=b.real();
}

// --generic/--generic-seed: random unitary change of basis on C^3,
// applied to the seed mesh's own vertex positions (NOT to F, which is
// left exactly as catalogued) to break any alignment between the
// curve's fixed X,Y,Z basis and the seed mesh's own -- copied verbatim
// from riemann_cp2.cpp's own random_unitary()/apply_unitary() (see
// that file's comment: 3 of the 15 Gaifullin seed points have TWO
// homogeneous coordinates equal to 0, a real, non-generic resonance
// with that basis, and everything basis-dependent downstream
// (vertex_gradient used by cell_priority, pick_chart/pick_chart_polygon)
// is evaluated in that same fixed basis).
void random_unitary(cx M[3][3], unsigned seed) {
	std::mt19937 rng(seed);
	std::normal_distribution<double> nd(0.0,1.0);
	pt3 cols[3];
	for(int c=0;c<3;++c) {
		pt3 v;
		for(int i=0;i<3;++i) v[i]=cx(nd(rng),nd(rng));
		for(int p=0;p<c;++p) {
			cx proj=hdot(v,cols[p]); // <v,cols[p]>, cols[p] already unit
			for(int i=0;i<3;++i) v[i]-=proj*cols[p][i];
		}
		cols[c]=normalize3(v);
	}
	for(int i=0;i<3;++i) for(int j=0;j<3;++j) M[i][j]=cols[j][i];
}
pt3 apply_unitary(const cx M[3][3], const pt3& v) {
	pt3 r;
	for(int i=0;i<3;++i) {
		r[i]=cx(0,0);
		for(int j=0;j<3;++j) r[i]+=M[i][j]*v[j];
	}
	return r;
}

// --- OBJ output vertex dedup: a crossing_node's identity (see
// glpt_extraction.hpp's own header comment) canonicalized into a
// single, std::map-orderable key -- dim==1's continuous param is
// quantized to CROSSING_MERGE_TOL, matching same_crossing_node()'s own
// tolerance; dim==2's packed face key (60 bits) is split across two
// ints, with `sub` as the tiebreaker distinguishing a face's up-to-2
// roots. Only ever holds one entry per ACTUAL distinct output vertex
// (a small fraction of the mesh, unlike glpt_edge_cache/
// face_crossing_cache), so std::map's own overhead doesn't matter here.
struct node_key {
	int dim, a, b;
	long long qparam;
	bool operator<(const node_key& o) const {
		if(dim!=o.dim) return dim<o.dim;
		if(a!=o.a) return a<o.a;
		if(b!=o.b) return b<o.b;
		return qparam<o.qparam;
	}
};
node_key key_of(const crossing_node& nd) {
	node_key k; k.dim=nd.dim; k.a=0; k.b=0; k.qparam=0;
	if(nd.dim==0) { k.a=nd.a; }
	else if(nd.dim==1) { k.a=nd.a; k.b=nd.b; k.qparam=(long long)(nd.param/CROSSING_MERGE_TOL+0.5); }
	else {
		// nd.fkey is the exact (hi,lo) face_key (glpt_crossing.hpp) --
		// hi packs the sorted triple's first two ids as 32+32 bits, lo
		// holds the third; folding lo and sub (0 or 1) together into
		// qparam via a bijective 2x+sub keeps all four quantities
		// (a,b,c,sub) exactly distinguished across k.a/k.b/k.qparam.
		k.a=int(nd.fkey.hi>>32); k.b=int(nd.fkey.hi&0xFFFFFFFFu);
		k.qparam=(long long)nd.fkey.lo*2+nd.sub;
	}
	return k;
}

// --- cell_priority: gradient-dispersion refinement criterion, ported
// from riemann_cp2.cpp's own cell_priority() (see that file's comment)
// -- unchanged in spirit, just reading a glpt cell's own vertex points
// and gradients (computed on demand, not cached, same as riemann_cp2.cpp
// itself does per its "on-demand replacements for the per-vertex
// Fx/Fy/Fz/Fval cache" comment) instead of an nmt<4> cell's.
pt3 vertex_gradient(const pt3& p) {
	pt3 g;
	g[0]=eval_poly3(g_Fx,p[0],p[1],p[2]);
	g[1]=eval_poly3(g_Fy,p[0],p[1],p[2]);
	g[2]=eval_poly3(g_Fz,p[0],p[1],p[2]);
	return g;
}
cx vertex_Fval(const pt3& p) { return eval_poly3(g_F,p[0],p[1],p[2]); }

// Shared by cell_priority's own --proximity branch and
// cell_priority_bernstein (below): the cell's own Fubini-Study
// diameter, and (proximity_factor) the diam/(mind+diam) factor itself,
// mind=min over the cell's own vertices of the first-order distance-
// to-curve estimate |F(v)|/|gradF(v)| -- copied verbatim in spirit
// from riemann_cp2.cpp's own cell_diam()/proximity_factor().
double cell_diam(const pt3 pts[glpt::DIM+1]) {
	double diam=0;
	for(int i=0;i<=glpt::DIM;++i)
		for(int j=i+1;j<=glpt::DIM;++j)
			diam=std::max(diam,fs_dist(pts[i],pts[j]));
	return diam;
}
double proximity_factor(const pt3 pts[glpt::DIM+1]) {
	double mind=1e300;
	for(int i=0;i<=glpt::DIM;++i) {
		pt3 g=vertex_gradient(pts[i]);
		double gn=hnorm(g);
		if(gn<1e-12) continue;
		double dv=std::abs(vertex_Fval(pts[i]))/gn;
		if(dv<mind) mind=dv;
	}
	double diam=cell_diam(pts);
	return diam/(mind+diam);
}

// --proximity: optional extra factor biasing refinement toward the
// curve itself (gradient dispersion alone is a global signal, evaluated
// identically whether or not a cell is anywhere near F=0) -- copied
// verbatim in spirit from riemann_cp2.cpp's own cell_priority() (see
// that file's own comment on g_proximity for the full derivation of
// diam/(mind+diam)).
bool g_proximity=false;
double cell_priority(const pt3 pts[glpt::DIM+1]) {
	pt3 ghat[glpt::DIM+1];
	int n=0;
	double mind=1e300;
	for(int i=0;i<=glpt::DIM;++i) {
		pt3 g = vertex_gradient(pts[i]);
		double gn = hnorm(g);
		if(gn<1e-12) continue; // near a singular point: can't normalize
		if(g_proximity) {
			double dv=std::abs(vertex_Fval(pts[i]))/gn;
			if(dv<mind) mind=dv;
		}
		for(int c=0;c<3;++c) g[c]/=gn;
		ghat[n++]=g;
	}
	if(n<2) return 1e18; // most of the cell sits right at a singular point: force refinement
	double S=0;
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			S += std::norm(hdot(ghat[i],ghat[j]));
	double kappa = S/(double(n)*double(n));
	double priority = 1.0-kappa;
	if(g_proximity) {
		double diam=cell_diam(pts);
		priority *= diam/(mind+diam);
	}
	return priority;
}

// --- Bernstein-certified refinement/pruning (--bernstein), a second,
// independent refinement criterion alongside gradient dispersion above
// -- ported from riemann_cp2.cpp's own cell_priority_bernstein() (see
// that file's comment): the widest candidate face box among a cell's
// 10 distinct 2-faces (glpt_crossing.hpp's own GLPT_CELL_FACES) that
// bernstein_bounds_cache can't certify empty, or -1 (a sentinel meaning
// "certified empty -- drop this cell, never refine it again") if none
// of them can contain a zero. Divides by the cell's own diameter for
// the same reason --proximity does (raw box width isn't comparable
// across refinement levels -- see that file's own comment on
// cell_priority_bernstein), and, like ordinary cell_priority, can be
// further scaled by --proximity's own diam/(mind+diam) factor.
double cell_priority_bernstein(const pt3 pts[glpt::DIM+1], const int ids[glpt::DIM+1], bernstein_bounds_cache& bcache) {
	double best=-1.0;
	for(int f=0; f<10; ++f) {
		pt3 face_pts[3]; int face_ids[3];
		for(int k=0;k<3;++k) { face_pts[k]=pts[GLPT_CELL_FACES[f][k]]; face_ids[k]=ids[GLPT_CELL_FACES[f][k]]; }
		const bernstein_bounds_cache::bbounds& b = bcache.get(face_pts, face_ids);
		if(bernstein_bounds_cache::maybe_zero(b)) {
			double w=(b.reHi-b.reLo)+(b.imHi-b.imLo);
			if(w>best) best=w;
		}
	}
	if(best>=0) {
		double diam=cell_diam(pts);
		if(diam>1e-12) best/=diam;
		if(g_proximity) best*=proximity_factor(pts);
	}
	return best;
}

// --bernstein-selftest: direct numerical verification, independent of
// the mesh, that (a) level-0 and level-1 enclosures both actually
// contain the TRUE range of Re(F)/Im(F) over a triangle's own CHART
// domain (the same one compute_face_crossing searches), sampled by
// brute force, and (b) level-1's box is a subset of level-0's
// (mathematically required: subdivision only tightens). Copied
// verbatim in spirit from riemann_cp2.cpp's own bernstein_selftest().
void bernstein_selftest() {
	pt3 P0=normalize3(mkpt(cx(1,0),cx(0.3,0.1),cx(-0.2,0.4)));
	pt3 P1=normalize3(mkpt(cx(0.2,-0.5),cx(1,0),cx(0.1,0.3)));
	pt3 P2=normalize3(mkpt(cx(-0.3,0.2),cx(0.4,-0.1),cx(1,0)));
	pt3 p[3]={P0,P1,P2};
	int chart=pick_chart(p);
	cpt Q0,Q1,Q2;
	dehomogenize(chart,P0,Q0.a,Q0.b);
	dehomogenize(chart,P1,Q1.a,Q1.b);
	dehomogenize(chart,P2,Q2.a,Q2.b);

	double l0Lo,l0Hi,i0Lo,i0Hi, l1Lo,l1Hi,i1Lo,i1Hi;
	bernstein_bounds_recursive(g_F,chart,Q0,Q1,Q2,0,l0Lo,l0Hi,i0Lo,i0Hi);
	bernstein_bounds_recursive(g_F,chart,Q0,Q1,Q2,1,l1Lo,l1Hi,i1Lo,i1Hi);

	double trueReLo=1e300,trueReHi=-1e300,trueImLo=1e300,trueImHi=-1e300;
	mt19937 rng(1);
	uniform_real_distribution<double> ud(0.0,1.0);
	int N=2000000;
	for(int s=0;s<N;++s) {
		double u=ud(rng), v=ud(rng);
		if(u+v>1.0) { u=1.0-u; v=1.0-v; }
		double L0=1.0-u-v, L1=u, L2=v;
		cx a=L0*Q0.a+L1*Q1.a+L2*Q2.a;
		cx b=L0*Q0.b+L1*Q1.b+L2*Q2.b;
		pt3 P=rehomogenize(chart,a,b);
		cx val=eval_poly3(g_F,P[0],P[1],P[2]);
		double re=val.real(), im=val.imag();
		if(re<trueReLo) trueReLo=re; if(re>trueReHi) trueReHi=re;
		if(im<trueImLo) trueImLo=im; if(im>trueImHi) trueImHi=im;
	}

	cout<<"chart: "<<chart<<endl;
	cout<<"level 0 box:  Re["<<l0Lo<<","<<l0Hi<<"]  Im["<<i0Lo<<","<<i0Hi<<"]"<<endl;
	cout<<"level 1 box:  Re["<<l1Lo<<","<<l1Hi<<"]  Im["<<i1Lo<<","<<i1Hi<<"]"<<endl;
	cout<<"true range (brute force, "<<N<<" samples): Re["<<trueReLo<<","<<trueReHi
			<<"]  Im["<<trueImLo<<","<<trueImHi<<"]"<<endl;

	bool ok=true;
	if(!(l0Lo<=trueReLo+1e-9 && l0Hi>=trueReHi-1e-9)) { cout<<"FAIL: level-0 Re doesn't enclose true range"<<endl; ok=false; }
	if(!(i0Lo<=trueImLo+1e-9 && i0Hi>=trueImHi-1e-9)) { cout<<"FAIL: level-0 Im doesn't enclose true range"<<endl; ok=false; }
	if(!(l1Lo<=trueReLo+1e-9 && l1Hi>=trueReHi-1e-9)) { cout<<"FAIL: level-1 Re doesn't enclose true range"<<endl; ok=false; }
	if(!(i1Lo<=trueImLo+1e-9 && i1Hi>=trueImHi-1e-9)) { cout<<"FAIL: level-1 Im doesn't enclose true range"<<endl; ok=false; }
	if(!(l1Lo>=l0Lo-1e-9 && l1Hi<=l0Hi+1e-9)) { cout<<"FAIL: level-1 Re box is not a subset of level-0's"<<endl; ok=false; }
	if(!(i1Lo>=i0Lo-1e-9 && i1Hi<=i0Hi+1e-9)) { cout<<"FAIL: level-1 Im box is not a subset of level-0's"<<endl; ok=false; }
	cout<<(ok?"ALL CHECKS PASSED":"CHECKS FAILED")<<endl;
}

// --- Cell vertex points/ids, computed from scratch via the cell's own
// bisection path (see glpt_points.hpp/glpt_vertex_ids.hpp -- both are
// inherently path-dependent, not derivable from a bare glpt code).
void path_from_root(const glpt& c, int path[], int& path_len) {
	path_len = 0; glpt cur = c;
	while(!cur.is_root()) { path[path_len++] = cur.is_child0()?0:1; cur = cur.parent(); }
	for(int i=0;i<path_len/2;++i) { int t=path[i]; path[i]=path[path_len-1-i]; path[path_len-1-i]=t; }
}
void cell_points_and_ids(const glpt& c, const vector<pt3>& gp, glpt_edge_cache& cache,
		pt3 points[glpt::DIM+1], int ids[glpt::DIM+1]) {
	int path[64], path_len;
	path_from_root(c, path, path_len);
	glpt_root_vertex_points(c.seed_index(), gp, points);
	glpt_root_vertex_ids(c.seed_index(), ids);
	glpt cur(c.seed_index());
	for(int i=0;i<path_len;++i) {
		pt3 child_points[glpt::DIM+1]; int child_ids[glpt::DIM+1];
		glpt_child_vertex_points(cur, points, path[i], child_points);
		glpt_child_vertex_ids(cur, ids, path[i], cache, child_ids);
		cur = cur.child(path[i]);
		for(int k=0;k<=glpt::DIM;++k) { points[k]=child_points[k]; ids[k]=child_ids[k]; }
	}
}

// --genus-check: combinatorial Euler characteristic of the extracted
// complex IN CP^2 -- vertices identified by crossing-node identity
// (node_key), independent of any projection. Bad cells leave holes;
// pinched vertices (link with k>1 components) are split into k copies;
// then, if every hole is a disc, g = sum g_i = (2c - b - chi)/2 with c
// components and b boundary loops.
bool g_genus_check=false;
static int dsu_find(std::vector<int>& p, int x) { while(p[x]!=x) { p[x]=p[p[x]]; x=p[x]; } return x; }
// Fubini-Study area of the extracted complex (each polygon fanned into
// triangles, each triangle measured in the horizontal tangent space at
// its first vertex; second-order accurate) against Wirtinger's theorem:
// with d(p,q)=arccos|<p,q>| a projective line is a sphere of radius 1/2,
// so a smooth curve of degree n has area n*pi. Also |F| at the vertices
// (unit representatives).
void surface_stats(const std::vector<pt3>& P, const std::vector<std::vector<int> >& faces, int degree) {
	double area=0;
	for(const std::vector<int>& f: faces) {
		const pt3& a=P[f[0]];
		for(size_t i=1;i+1<f.size();++i) {
			pt3 b=align_phase(a,P[f[i]]), c=align_phase(a,P[f[i+1]]), e1, e2;
			for(int k=0;k<3;++k) { e1[k]=b[k]-a[k]; e2[k]=c[k]-a[k]; }
			cx p1=hdot(e1,a), p2=hdot(e2,a);
			for(int k=0;k<3;++k) { e1[k]-=p1*a[k]; e2[k]-=p2*a[k]; }
			double n1=hdot(e1,e1).real(), n2=hdot(e2,e2).real(), r=hdot(e1,e2).real();
			area+=0.5*std::sqrt(std::max(0.0,n1*n2-r*r));
		}
	}
	std::vector<double> res;
	for(const pt3& p: P) { pt3 u=normalize3(p); res.push_back(std::abs(eval_poly3(g_F,u[0],u[1],u[2]))); }
	std::sort(res.begin(),res.end());
	std::cout<<"FS area = "<<area<<"  (Wirtinger n*pi = "<<degree*M_PI<<", ratio "<<area/(degree*M_PI)<<")"<<std::endl;
	if(!res.empty())
		std::cout<<"|F| at "<<res.size()<<" vertices: median="<<res[res.size()/2]<<" p99="<<res[(99*res.size())/100]<<" max="<<res.back()<<std::endl;
}
void genus_check(const std::vector<std::vector<int> >& faces_in, int nverts, int expected_genus) {
	// corners grouped by link component of their vertex
	std::vector<std::vector<std::pair<int,int> > > corners(nverts); // (face, pos)
	for(size_t f=0; f<faces_in.size(); ++f)
		for(size_t i=0; i<faces_in[f].size(); ++i) corners[faces_in[f][i]].push_back(std::make_pair((int)f,(int)i));
	std::vector<std::vector<int> > faces(faces_in);
	int nid=0, pinched=0;
	for(int v=0; v<nverts; ++v) {
		const std::vector<std::pair<int,int> >& cs=corners[v];
		if(cs.empty()) continue;
		// link vertices of v: union-find over neighbour ids
		std::map<int,int> li; std::vector<int> par;
		auto idx=[&](int x){ std::map<int,int>::iterator it=li.find(x); if(it!=li.end()) return it->second; int k=(int)par.size(); li[x]=k; par.push_back(k); return k; };
		std::vector<int> ca(cs.size());
		for(size_t j=0;j<cs.size();++j) {
			const std::vector<int>& F=faces_in[cs[j].first]; int n=(int)F.size(), i=cs[j].second;
			int a=idx(F[(i+n-1)%n]), b=idx(F[(i+1)%n]);
			int ra=dsu_find(par,a), rb=dsu_find(par,b); if(ra!=rb) par[ra]=rb;
			ca[j]=a;
		}
		std::map<int,int> comp;
		for(size_t j=0;j<cs.size();++j) {
			int r=dsu_find(par,ca[j]);
			if(!comp.count(r)) { int k=(int)comp.size(); comp[r]=k; }
		}
		if(comp.size()>1) ++pinched;
		for(size_t j=0;j<cs.size();++j) faces[cs[j].first][cs[j].second]=nid+comp[dsu_find(par,ca[j])];
		nid+=(int)comp.size();
	}
	std::map<std::pair<int,int>,int> edges;
	for(size_t f=0; f<faces.size(); ++f) { int n=(int)faces[f].size();
		for(int i=0;i<n;++i) { int a=faces[f][i], b=faces[f][(i+1)%n]; edges[std::make_pair(std::min(a,b),std::max(a,b))]++; } }
	long V=nid, E=(long)edges.size(), F=(long)faces.size(); long chi=V-E+F;
	std::vector<int> pb(nid), ps(nid);
	for(int i=0;i<nid;++i) pb[i]=ps[i]=i;
	std::vector<char> onb(nid,0); long nonman=0;
	for(auto& e: edges) {
		if(e.second==1) { onb[e.first.first]=onb[e.first.second]=1;
			int r1=dsu_find(pb,e.first.first), r2=dsu_find(pb,e.first.second); if(r1!=r2) pb[r1]=r2; }
		if(e.second>2) ++nonman;
		int s1=dsu_find(ps,e.first.first), s2=dsu_find(ps,e.first.second); if(s1!=s2) ps[s1]=s2;
	}
	std::set<int> loops, comps;
	for(int i=0;i<nid;++i) { if(onb[i]) loops.insert(dsu_find(pb,i)); comps.insert(dsu_find(ps,i)); }
	long b=(long)loops.size(), c=(long)comps.size();
	cout<<endl<<"--- genus check (complex in CP^2, projection-independent) ---"<<endl;
	cout<<"pinched vertices split: "<<pinched<<", non-manifold edges: "<<nonman<<", components c="<<c<<endl;
	cout<<"V="<<V<<" E="<<E<<" F="<<F<<" chi="<<chi<<" boundary loops b="<<b<<endl;
	cout<<"genus estimate (2c-b-chi)/2 = "<<(2.0*c-b-chi)/2.0;
	if(expected_genus>=0) cout<<"   (expected "<<expected_genus<<" = (n-1)(n-2)/2)";
	cout<<endl;
}


// --close: topological post-processing of the extracted complex in CP^2.
// (1) split pinched vertices (one copy per link component); (2) keep the
// main connected component (tiny islands cut off by the pinches are
// dropped); (3) cap every boundary loop -- a hole left by bad cells --
// with a cone from a NEW vertex (a fan from an existing vertex could
// duplicate an edge), placed at the loop's Fubini-Study barycentre and
// moved onto the curve by Newton's method in the best affine chart.
// The result is verified by genus_check() (closed, connected, chi=2-2g).
bool g_close=false;
bool g_cause_stats=false;
bool g_global_scan=false; // continuation from all 108 roots instead of the first touched one
bool g_timing=false;      // phase timings + peak RSS
bool g_leaf_hash=false;   // order-independent hash of the final leaf set
bool g_neighbor_test=false; // every leaf's facet neighbours must carry the facet's 4 vertices, at the same points
// --kahler-filter X: omega/area of an extracted polygon, omega the
// Fubini-Study Kahler form, on the fan at its first vertex with tangent
// vectors in the horizontal space. By Wirtinger it is +-1 exactly on
// complex lines, so a polygon close to the curve has |ratio| near 1; a cell
// with a polygon below X is treated as bad (examples/top/riemann_s2s2_glpt
// has the same test, where it removes polygons cutting across a cusp).
double g_kahler_filter=0;
long g_kahler_rejected=0;
double polygon_kahler_ratio(const std::vector<crossing_node>& nodes, const std::vector<int>& cyc) {
	auto tang=[](const pt3& a, const pt3& b){ pt3 e=align_phase(a,b); for(int k=0;k<3;++k) e[k]-=a[k];
		cx c=hdot(e,a); for(int k=0;k<3;++k) e[k]-=c*a[k]; return e; };
	pt3 a=normalize3(nodes[cyc[0]].p); double om=0, ar=0;
	for(size_t i=1;i+1<cyc.size();++i) {
		pt3 e1=tang(a,normalize3(nodes[cyc[i]].p)), e2=tang(a,normalize3(nodes[cyc[i+1]].p));
		cx h=hdot(e2,e1);
		double n1=hdot(e1,e1).real(), n2=hdot(e2,e2).real();
		om+=0.5*h.imag(); ar+=0.5*std::sqrt(std::max(0.0,n1*n2-h.real()*h.real()));
	}
	return ar>0 ? om/ar : 1.0;
}
// --box R: the affine-chart baseline. The seed is Kuhn's triangulation of
// the cube [-R,R]^4 in the coordinates (Re x, Im x, Re y, Im y) of the
// chart Z=1: 16 corners (id = bitmask of the coordinates at +R) and 24
// cells, one per permutation, each ordered by the number of coordinates
// at +R -- a balanced, colour-ordered seed with boundary. New vertices
// are affine midpoints in the chart and every face is flat in it, so
// everything is consistent within the one chart.
double g_box_R = 0;
// F'(p) = F(M p): the curve seen from a seed whose points were NOT moved
// by M. --generic moves the Gaifullin seed by M; with --box the seed is
// kept axis-aligned in its chart and the curve is moved instead, which
// gives the same relative position of curve and mesh (and keeps the
// curve from meeting the box grid in non-generic ways: with real
// coefficients its real locus would lie in a union of mesh faces).
inline poly_F3 compose_linear(const poly_F3& F, const cx M[3][3]) {
	typedef std::map<std::array<int,3>,cx> P3;
	auto mul=[](const P3& A, const P3& B){ P3 R; for(auto& a: A) for(auto& b: B) { std::array<int,3> e{{a.first[0]+b.first[0],a.first[1]+b.first[1],a.first[2]+b.first[2]}}; R[e]+=a.second*b.second; } return R; };
	P3 L[3];
	for(int i=0;i<3;++i) for(int j=0;j<3;++j) { std::array<int,3> e{{0,0,0}}; e[j]=1; L[i][e]=M[i][j]; }
	P3 out;
	for(const term3& t: F.t) {
		P3 m; m[std::array<int,3>{{0,0,0}}]=t.c;
		for(int i=0;i<3;++i) for(int k=0;k<t.e[i];++k) m=mul(m,L[i]);
		for(auto& kv: m) out[kv.first]+=kv.second;
	}
	poly_F3 G; G.d=F.d;
	for(auto& kv: out) if(std::abs(kv.second)>1e-15) { term3 t; t.e[0]=kv.first[0]; t.e[1]=kv.first[1]; t.e[2]=kv.first[2]; t.c=kv.second; G.t.push_back(t); }
	return G;
}
static int g_box_cells[24][5];
inline void make_box_seed(double R, vector<pt3>& gp) {
	gp.assign(16, pt3());
	for(int m=0;m<16;++m) {
		double x[4]; for(int k=0;k<4;++k) x[k] = (m>>k & 1) ? R : -R;
		pt3 p; p[0]=cx(x[0],x[1]); p[1]=cx(x[2],x[3]); p[2]=cx(1,0);
		gp[m]=normalize3(p);
	}
	int perm[4]={0,1,2,3}, c=0;
	do {
		int v=0; g_box_cells[c][0]=0;
		for(int k=0;k<4;++k) { v|=1<<perm[k]; g_box_cells[c][k+1]=v; }
		++c;
	} while(std::next_permutation(perm,perm+4));
	glpt_set_seed(g_box_cells, 24, 16, true);
	g_affine_midpoint_chart = 2;
	g_forced_face_chart = 2;
	g_geodesic_faces = false;
}
// --fs-diam h: refine every cell that meets C until its Fubini-Study
// diameter is below h (instead of, or on top of, a fixed --depth);
// gives the intrinsic and the affine-chart meshes the same resolution.
double g_fs_diam = 0;
inline double cell_fs_diam(const pt3 pts[glpt::DIM+1]) {
	double d=0; for(int a=0;a<=glpt::DIM;++a) for(int b=a+1;b<=glpt::DIM;++b) d=std::max(d,fs_dist(pts[a],pts[b]));
	return d;
}
inline double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
inline long peak_rss_kb() { // VmHWM from /proc/self/status (Linux)
	std::ifstream f("/proc/self/status"); std::string line;
	while(std::getline(f,line)) if(line.rfind("VmHWM:",0)==0) return std::atol(line.c_str()+6);
	return -1;
}
static pt3 newton_onto_curve(pt3 p, bool& ok) {
	int k=0; for(int i=1;i<3;++i) if(std::abs(p[i])>std::abs(p[k])) k=i;
	int i1=(k+1)%3, i2=(k+2)%3; if(i1>i2) std::swap(i1,i2);
	const poly_F3* d[3]={&g_Fx,&g_Fy,&g_Fz};
	pt3 q=p; for(int t=0;t<3;++t) q[t]=p[t]/p[k];
	ok=false;
	for(int it=0; it<30; ++it) {
		cx f=eval_poly3(g_F,q[0],q[1],q[2]);
		cx ga=eval_poly3(*d[i1],q[0],q[1],q[2]), gb=eval_poly3(*d[i2],q[0],q[1],q[2]);
		double n2=std::norm(ga)+std::norm(gb);
		if(n2<1e-300) break;
		// minimum-norm Newton step for one complex equation in C^2
		q[i1]-=f*std::conj(ga)/n2; q[i2]-=f*std::conj(gb)/n2;
		double nq=hnorm(q);
		if(std::abs(eval_poly3(g_F,q[0],q[1],q[2]))/std::pow(nq,(double)g_F.d)<1e-13) { ok=true; break; }
	}
	return normalize3(q);
}
void close_surface(std::vector<pt3>& P, std::vector<std::vector<int> >& faces) {
	int nv=(int)P.size();
	// (0) drop every polygon on a non-manifold edge (one shared by more than
	// two polygons): each cell passed its own checks, but their polygons do
	// not form a surface there -- seen at near-tangencies, where a 2-face
	// carries two close crossings and several cells join them. The hole is
	// capped in (3) like the hole of a bad cell.
	int nonman_dropped=0;
	{
		std::map<std::pair<int,int>,int> ec;
		for(auto& f: faces) for(size_t i=0;i<f.size();++i) { int a=f[i], b=f[(i+1)%f.size()]; ec[std::make_pair(std::min(a,b),std::max(a,b))]++; }
		std::vector<std::vector<int> > keep;
		for(auto& f: faces) {
			bool bad=false;
			for(size_t i=0;i<f.size() && !bad;++i) { int a=f[i], b=f[(i+1)%f.size()]; if(ec[std::make_pair(std::min(a,b),std::max(a,b))]>2) bad=true; }
			if(bad) ++nonman_dropped; else keep.push_back(f);
		}
		faces.swap(keep);
	}
	// (1) split pinched vertices
	std::vector<std::vector<std::pair<int,int> > > corners(nv);
	for(size_t f=0; f<faces.size(); ++f)
		for(size_t i=0; i<faces[f].size(); ++i) corners[faces[f][i]].push_back(std::make_pair((int)f,(int)i));
	std::vector<std::vector<int> > nf(faces);
	std::vector<pt3> NP; int pinched=0;
	for(int v=0; v<nv; ++v) {
		const std::vector<std::pair<int,int> >& cs=corners[v];
		if(cs.empty()) continue;
		std::map<int,int> li; std::vector<int> par;
		auto idx=[&](int x){ std::map<int,int>::iterator it=li.find(x); if(it!=li.end()) return it->second; int k=(int)par.size(); li[x]=k; par.push_back(k); return k; };
		std::vector<int> ca(cs.size());
		for(size_t j=0;j<cs.size();++j) {
			const std::vector<int>& F=faces[cs[j].first]; int n=(int)F.size(), i=cs[j].second;
			int a=idx(F[(i+n-1)%n]), b=idx(F[(i+1)%n]);
			int ra=dsu_find(par,a), rb=dsu_find(par,b); if(ra!=rb) par[ra]=rb;
			ca[j]=a;
		}
		std::map<int,int> comp;
		for(size_t j=0;j<cs.size();++j) { int r=dsu_find(par,ca[j]); if(!comp.count(r)) { int k=(int)comp.size(); comp[r]=k; } }
		if(comp.size()>1) ++pinched;
		int base=(int)NP.size();
		for(size_t k=0;k<comp.size();++k) NP.push_back(P[v]);
		for(size_t j=0;j<cs.size();++j) nf[cs[j].first][cs[j].second]=base+comp[dsu_find(par,ca[j])];
	}
	// (2) keep the main component
	int n2=(int)NP.size();
	std::vector<int> ps(n2); for(int i=0;i<n2;++i) ps[i]=i;
	for(auto& f: nf) for(size_t i=1;i<f.size();++i) { int a=dsu_find(ps,f[0]), b=dsu_find(ps,f[i]); if(a!=b) ps[a]=b; }
	std::map<int,int> csize; for(auto& f: nf) csize[dsu_find(ps,f[0])]+=1;
	int mainr=-1, best=-1; for(auto& c: csize) if(c.second>best) { best=c.second; mainr=c.first; }
	std::vector<std::vector<int> > kept; int dropped=0;
	for(auto& f: nf) { if(dsu_find(ps,f[0])==mainr) kept.push_back(f); else ++dropped; }
	// (3) cap boundary loops with cones
	std::map<std::pair<int,int>,int> ec;
	for(auto& f: kept) for(size_t i=0;i<f.size();++i) { int a=f[i], b=f[(i+1)%f.size()]; ec[std::make_pair(std::min(a,b),std::max(a,b))]++; }
	std::vector<std::pair<int,int> > bnd;
	for(auto& f: kept) for(size_t i=0;i<f.size();++i) { int a=f[i], b=f[(i+1)%f.size()]; if(ec[std::make_pair(std::min(a,b),std::max(a,b))]==1) bnd.push_back(std::make_pair(b,a)); }
	std::vector<int> pb(n2); for(int i=0;i<n2;++i) pb[i]=i;
	for(auto& e: bnd) { int a=dsu_find(pb,e.first), b=dsu_find(pb,e.second); if(a!=b) pb[a]=b; }
	std::map<int,std::vector<int> > loopv;
	for(auto& e: bnd) loopv[dsu_find(pb,e.first)].push_back(e.first);
	std::map<int,int> centre; int onc=0;
	for(auto& L: loopv) {
		pt3 ref=NP[L.second[0]], sum; for(int t=0;t<3;++t) sum[t]=cx(0,0);
		for(int x: L.second) { pt3 q=align_phase(ref,NP[x]); for(int t=0;t<3;++t) sum[t]+=q[t]; }
		pt3 c=normalize3(sum);
		bool ok; pt3 cc=newton_onto_curve(c, ok);
		double r=0; for(int x: L.second) r=std::max(r,fs_dist(c,NP[x]));
		if(ok && fs_dist(c,cc)<=2*r+1e-12) { c=cc; ++onc; }
		centre[L.first]=(int)NP.size(); NP.push_back(c);
	}
	for(auto& e: bnd) kept.push_back(std::vector<int>{e.first, e.second, centre[dsu_find(pb,e.first)]});
	// compact vertex ids
	std::vector<int> remap(NP.size(),-1); std::vector<pt3> out;
	for(auto& f: kept) for(int& x: f) { if(remap[x]<0) { remap[x]=(int)out.size(); out.push_back(NP[x]); } x=remap[x]; }
	P.swap(out); faces.swap(kept);
	cout<<endl<<"--- close (topological post-processing in CP^2) ---"<<endl;
	cout<<"polygons dropped at non-manifold edges: "<<nonman_dropped<<endl;
	cout<<"pinched vertices split: "<<pinched<<", island faces dropped: "<<dropped
		<<", holes capped: "<<loopv.size()<<" (cone apex moved onto the curve: "<<onc<<")"<<endl;
}

int main(int argc, char* argv[]) {
	int function_idx=0;
	int max_depth=8;
	double threshold=0.05;
	string obj_path="riemann_cp2_glpt.obj";
	// projection: same default/flags as riemann_cp2.cpp's own --flat
	// (implicit default)/--alpha-projection/--cutoff (onion deliberately
	// not ported here -- see this file's own header comment).
	bool alpha_mode=false;
	double cutoff=1e300;
	bool generic=false;
	unsigned generic_seed=12345;
	bool bernstein_selftest_flag=false;
	bool certify_selftest_flag=false;
	int repair_rounds=0;
	int repair_extra_depth=4;
	bool legacy_priority_refine=false;
	double tangency_threshold=0.0; // 0 = disabled
	int tangency_extra_depth=4;

	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--function" && i+1<argc) { function_idx=atoi(argv[++i]); }
		else if(arg=="--depth" && i+1<argc) { max_depth=atoi(argv[++i]); }
		else if(arg=="--threshold" && i+1<argc) { threshold=atof(argv[++i]); }
		else if(arg=="--obj" && i+1<argc) { obj_path=argv[++i]; }
		else if(arg=="--flat") { alpha_mode=false; }
		else if(arg=="--alpha-projection") { alpha_mode=true; }
		else if(arg=="--cutoff" && i+1<argc) { cutoff=atof(argv[++i]); }
		else if(arg=="--chart" && i+1<argc) { g_forced_chart=atoi(argv[++i]); }
		else if(arg=="--flat-swap") { g_flat_swap=true; }
		else if(arg=="--genus-check") { g_genus_check=true; }
		else if(arg=="--close") { g_close=true; g_genus_check=true; }
		else if(arg=="--generic") { generic=true; }
		else if(arg=="--generic-seed" && i+1<argc) { generic=true; generic_seed=(unsigned)atoi(argv[++i]); }
		else if(arg=="--proximity") { g_proximity=true; }
		else if(arg=="--bernstein") { g_bernstein=true; }
		else if(arg=="--bernstein-level" && i+1<argc) { g_bernstein=true; g_bernstein_level=atoi(argv[++i]); }
		else if(arg=="--bernstein-selftest") { bernstein_selftest_flag=true; }
		else if(arg=="--certify") { g_certify=true; }
		else if(arg=="--geodesic-faces") { g_geodesic_faces=true; } // the default; accepted for old command lines
		else if(arg=="--chart-flat-faces") { g_geodesic_faces=false; }
		else if(arg=="--cause-stats") { g_cause_stats=true; }
		else if(arg=="--global-scan") { g_global_scan=true; }
		else if(arg=="--box" && i+1<argc) { g_box_R=atof(argv[++i]); }
		else if(arg=="--neighbor-test") { g_neighbor_test=true; }
		else if(arg=="--kahler-filter" && i+1<argc) { g_kahler_filter=atof(argv[++i]); }
		else if(arg=="--fs-diam" && i+1<argc) { g_fs_diam=atof(argv[++i]); }
		else if(arg=="--timing") { g_timing=true; }
		else if(arg=="--leaf-hash") { g_leaf_hash=true; }
		else if(arg=="--certify-level" && i+1<argc) { g_certify=true; g_certify_level=atoi(argv[++i]); }
		else if(arg=="--certify-selftest") { g_certify=true; certify_selftest_flag=true; }
		else if(arg=="--repair-rounds" && i+1<argc) { repair_rounds=atoi(argv[++i]); }
		else if(arg=="--repair-extra-depth" && i+1<argc) { repair_extra_depth=atoi(argv[++i]); }
		else if(arg=="--tangency-threshold" && i+1<argc) { tangency_threshold=atof(argv[++i]); }
		else if(arg=="--tangency-extra-depth" && i+1<argc) { tangency_extra_depth=atoi(argv[++i]); }
		else if(arg=="--legacy-priority-refine") { legacy_priority_refine=true; }
		else if(arg=="--list") { print_function_catalog_cp2(cout); return 0; }
		else {
			cerr<<"usage: "<<argv[0]<<" [--function N] [--depth N] [--threshold X] [--obj PATH] "
				<<"[--flat] [--alpha-projection] [--cutoff X] [--chart N] [--flat-swap] [--genus-check] [--close] [--generic] [--generic-seed N] [--proximity] "
				<<"[--bernstein] [--bernstein-level N] [--bernstein-selftest] [--certify] [--certify-level N] [--certify-selftest] [--geodesic-faces | --chart-flat-faces] [--cause-stats] [--global-scan] [--timing] [--leaf-hash] [--box R] [--fs-diam h] [--neighbor-test] "
				<<"[--repair-rounds N] [--repair-extra-depth N] "
				<<"[--tangency-threshold X] [--tangency-extra-depth N] [--legacy-priority-refine] [--list]"<<endl;
			if(arg!="--help" && arg!="-h") return 1;
			return 0;
		}
	}
	if(g_proximity && !legacy_priority_refine)
		cerr<<"note: --proximity has no effect under the default continuation refinement "
			<<"(no priority queue) -- pass --legacy-priority-refine to use it"<<endl;

	vector<catalog_entry_cp2>& cat = function_catalog_cp2();
	if(function_idx<0 || function_idx>=(int)cat.size()) {
		cerr<<"bad --function index "<<function_idx<<" (--list to see the catalog)"<<endl;
		return 1;
	}
	set_curve(cat[function_idx].F);

	if(certify_selftest_flag) return certify_selftest() ? 0 : 1;
	if(g_geodesic_faces && g_bernstein) {
		cerr<<"--bernstein assumes chart-flat faces: pass --chart-flat-faces with it (geodesic faces are the default)"<<endl;
		return 1;
	}
	if(g_geodesic_faces && g_certify)
		cerr<<"note: --certify tests the cell's chart-flat simplex, which does not match the default geodesic faces"<<endl;
	if(bernstein_selftest_flag) {
		bernstein_selftest();
		return 0;
	}

	cout<<"curve: "<<cat[function_idx].name<<" -- "<<cat[function_idx].description<<endl;
	cout<<"generic="<<(generic?"on":"off");
	if(generic) cout<<" (seed="<<generic_seed<<")";
	cout<<"  proximity="<<(g_proximity?"on":"off");
	cout<<"  bernstein="<<(g_bernstein?"on":"off");
	cout<<"  faces="<<(g_box_R>0 ? "chart-flat (box)" : (g_geodesic_faces?"geodesic":"chart-flat"));
	if(g_bernstein) cout<<" (level="<<g_bernstein_level<<")";
	cout<<endl;

	bernstein_bounds_cache bcache;
	if(g_bernstein) enable_bernstein(bcache);

	if(g_box_R>0 && (g_bernstein || g_close)) {
		cerr<<"--box is incompatible with --bernstein and --close"<<endl;
		return 1;
	}
	vector<pt3> gp = gaifullin_points();
	if(g_box_R>0) {
		make_box_seed(g_box_R, gp);
		cout<<"seed: Kuhn box [-"<<g_box_R<<","<<g_box_R<<"]^4 in the chart Z=1 (affine midpoints, chart-flat faces)"<<endl;
	}
	if(generic) {
		cx M[3][3]; random_unitary(M,generic_seed);
		if(g_box_R>0) set_curve(compose_linear(cat[function_idx].F, M)); // see compose_linear()
		else for(size_t i=0;i<gp.size();++i) gp[i]=apply_unitary(M,gp[i]);
	}
	glpt_edge_cache id_cache;
	glpt_tree tree;
	tree.seed_all_roots();
	cout<<"seeded "<<tree.leaf_count()<<" root cells"<<endl;

	// --- Phase 1: adaptive refinement --
	//
	// Default (continuation): find ONE root cell the curve certainly
	// passes through (extract_cell()'s own any_edge, tested directly --
	// not a proxy like gradient dispersion), then explore outward
	// through ONLY the cells actually touched by the curve: bisect a
	// touched cell down to max_depth via compat_bisect() (recursing into
	// every cell that cascade creates, since a graded-mesh cascade can
	// touch cells elsewhere too -- each independently re-tested by
	// any_edge, so an irrelevant collateral cell is dropped immediately);
	// once at max_depth, propagate to glpt_tree::neighbor_leaf() across
	// all 5 facets. Never looks at the other 107 roots at all unless
	// reached this way. For a CONNECTED curve this was measured
	// (2026-09-25) to recover the EXACT SAME final leaf set as testing
	// every one of the 108 roots directly, while building a mesh 2.5-6x
	// smaller than this file's own gradient-dispersion Phase 1 at the
	// same depth/--threshold, and never a WORSE bad-cell rate (only ever
	// better, on curves where the old criterion stopped refining a
	// touched cell too early) -- see [[riemann_glpt_continuation]].
	// face_crossing_cache is shared with Phase 2 below on purpose: every
	// 2-face this phase already solved via any_edge testing is reused
	// there, not recomputed.
	//
	// --legacy-priority-refine: the ORIGINAL global gradient-dispersion
	// priority queue (cell_priority/cell_priority_bernstein, modulated
	// by --proximity), kept only for comparison -- see this file's own
	// git history and [[riemann_pc2_gradient_dispersion]].
	double t_refine0=now_s();
	cout<<endl<<"--- adaptive refinement ---"<<endl;
	face_crossing_cache fcache;
	int nsubdivisions=0;
	if(legacy_priority_refine) {
		priority_queue<pair<double,glpt> > pq;
		int ndropped=0;
		for(int s=0;s<glpt_seed_count();++s) {
			pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
			glpt c(s);
			cell_points_and_ids(c, gp, id_cache, pts, ids);
			double p = g_bernstein ? cell_priority_bernstein(pts,ids,bcache) : cell_priority(pts);
			if(g_bernstein && p<0) { ++ndropped; continue; } // certified empty
			pq.push(make_pair(p, c));
		}
		while(!pq.empty()) {
			pair<double,glpt> top = pq.top(); pq.pop();
			glpt cv = top.second;
			if(!tree.exists(cv)) continue;
			if(top.first<threshold) break;
			if(cv.simplex_level()>=max_depth) continue;
			tree.clear_recent();
			tree.compat_bisect(cv);
			++nsubdivisions;
			const vector<glpt>& recent = tree.recent_leaves();
			for(size_t i=0;i<recent.size();++i) {
				if(!tree.exists(recent[i])) continue; // superseded later in this same cascade -- see recent_leaves()'s own comment
				pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
				cell_points_and_ids(recent[i], gp, id_cache, pts, ids);
				double p = g_bernstein ? cell_priority_bernstein(pts,ids,bcache) : cell_priority(pts);
				if(g_bernstein && p<0) { ++ndropped; continue; }
				pq.push(make_pair(p, recent[i]));
			}
		}
		cout<<"subdivisions performed: "<<nsubdivisions<<endl;
		if(g_bernstein) cout<<"cells certified empty (dropped, never refined): "<<ndropped<<endl;
	} else {
		int seed_root=-1;
		long n_cells_tested=0;
		long n_tangency_refined=0;
		if(g_global_scan) seed_root=0; // every root is pushed below
		else for(int s=0;s<glpt_seed_count();++s) {
			pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
			cell_points_and_ids(glpt(s), gp, id_cache, pts, ids);
			cell_extraction_result res;
			extract_cell(pts, ids, fcache, res);
			++n_cells_tested;
			if(res.any_edge) { seed_root=s; break; }
		}
		if(seed_root<0) {
			cout<<"curve doesn't touch any root cell -- nothing to mesh"<<endl;
		} else {
			if(g_global_scan) cout<<"global scan: continuation from all "<<glpt_seed_count()<<" roots"<<endl;
			else cout<<"seed root: "<<seed_root<<" (scanned "<<(seed_root+1)<<"/"<<glpt_seed_count()<<")"<<endl;
			set<glpt> visited;
			vector<glpt> frontier(1, glpt(seed_root));
			if(g_global_scan) { frontier.clear(); for(int r=glpt_seed_count()-1;r>=0;--r) frontier.push_back(glpt(r)); }
			while(!frontier.empty()) {
				glpt c = frontier.back(); frontier.pop_back();
				if(!tree.exists(c)) continue;      // superseded by an earlier cascade
				if(visited.count(c)) continue;
				visited.insert(c);
				pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
				cell_points_and_ids(c, gp, id_cache, pts, ids);
				cell_extraction_result res;
				extract_cell(pts, ids, fcache, res);
				++n_cells_tested;
				if(!res.any_edge) continue;        // curve doesn't reach here -- stop
				// --tangency-threshold: a cell whose min_transversality (see
				// glpt_extraction.hpp's own comment) is suspiciously close to
				// 0 -- measured to correlate strongly with bad-cell risk --
				// gets bisected further even past --depth, up to a separate
				// --tangency-extra-depth budget, exactly mirroring
				// --repair-rounds' extra-depth philosophy but applied
				// PROACTIVELY here, before extraction, rather than only after
				// decompose_ok already failed.
				bool needs_tangency_refine = tangency_threshold>0.0
					&& res.min_transversality<tangency_threshold
					&& c.simplex_level()<max_depth+tangency_extra_depth;
				if(needs_tangency_refine && c.simplex_level()>=max_depth) ++n_tangency_refined;
				bool needs_fs_refine = g_fs_diam>0 && c.simplex_level()<40 && cell_fs_diam(pts)>g_fs_diam;
				if(c.simplex_level()<max_depth || needs_tangency_refine || needs_fs_refine) {
					tree.clear_recent();
					tree.compat_bisect(c);
					++nsubdivisions;
					const vector<glpt>& recent = tree.recent_leaves();
					for(size_t i=0;i<recent.size();++i) frontier.push_back(recent[i]);
				} else {
					for(int i=0;i<=glpt::DIM;++i) {
						glpt nb;
						glpt_tree::neighbor_status ns = tree.neighbor_leaf(c, i, nb);
						if(ns==glpt_tree::FOUND && !visited.count(nb)) frontier.push_back(nb);
					}
				}
			}
		}
		cout<<"subdivisions performed: "<<nsubdivisions<<", cells tested: "<<n_cells_tested<<endl;
		if(tangency_threshold>0.0) cout<<"cells past --depth forced deeper by --tangency-threshold: "<<n_tangency_refined<<endl;
	}
	cout<<"final leaf count: "<<tree.leaf_count()<<endl;
	cout<<"distinct vertices minted: "<<id_cache.next_id()<<" ("<<glpt_seed_vertex_count()<<" base + "
		<<(id_cache.next_id()-glpt_seed_vertex_count())<<" from bisection)"<<endl;

	// --- Phase 2: extraction -- actual connected polygons, not just
	// points (glpt_extraction.hpp's extract_cell(), one call per final
	// leaf; face_crossing_cache solves each DISTINCT 2-face -- by
	// vertex-id triple -- only once, however many cells' facets touch
	// it).
	//
	// Output vertex handling differs by projection mode:
	//   alpha (--alpha-projection): chart-free (see alpha_projection()'s
	//   own comment) -- a given point always projects the same way no
	//   matter which polygon references it, so vertices are deduped by
	//   crossing_node IDENTITY (node_key), not coordinate proximity: the
	//   same node reached from different cells gets the SAME output
	//   vertex index, and shared edges between adjacent cells' polygons
	//   connect exactly.
	//   flat (default, --flat): NOT chart-free -- riemann_cp2.cpp's own
	//   history (pick_chart_polygon's comment) found a real bug from
	//   letting the same point be dehomogenized in different (a,b)
	//   frames depending on which polygon it's part of ("ribbon spray").
	//   Its fix -- one best-conditioned chart per WHOLE polygon -- is
	//   only correct if vertices are NOT globally shared across
	//   polygons with different chart choices, so flat mode writes
	//   independent vertices per polygon, exactly like riemann_cp2.cpp's
	//   own OBJ output does (verified there directly: every output edge
	//   at multiplicity 1, no global dedup at all).
	if(g_neighbor_test) { // facet-neighbour consistency of the final mesh
		struct T { glpt_tree* tr; const vector<pt3>* gp; glpt_edge_cache* ic; long *nq,*nbad,*nposbad,*nbnd;
			void operator()(const glpt& c) const {
				pt3 P[5]; int I[5]; cell_points_and_ids(c,*gp,*ic,P,I);
				for(int i=0;i<5;++i) {
					glpt nb; glpt_tree::neighbor_status st=tr->neighbor_leaf(c,i,nb);
					if(st!=glpt_tree::FOUND) { if(st==glpt_tree::UNRESOLVED) ++*nbnd; continue; }
					++*nq;
					pt3 Q[5]; int J[5]; cell_points_and_ids(nb,*gp,*ic,Q,J);
					if(nb.simplex_level()!=c.simplex_level()) continue; // coarser neighbour: facet not shared exactly
					for(int k=0;k<5;++k) { if(k==i) continue; int m=-1; for(int l=0;l<5;++l) if(J[l]==I[k]) m=l;
						if(m<0) { ++*nbad; break; }
						if(fs_dist(P[k],Q[m])>1e-6) { ++*nposbad; break; } }
				}
			} };
		long nq=0,nbad=0,nposbad=0,nbnd=0; T t; t.tr=&tree; t.gp=&gp; t.ic=&id_cache; t.nq=&nq; t.nbad=&nbad; t.nposbad=&nposbad; t.nbnd=&nbnd;
		tree.for_each_leaf(t);
		cout<<"neighbor test: "<<nq<<" neighbour queries, "<<nbad<<" with a facet vertex missing, "<<nposbad<<" with a position mismatch, "<<nbnd<<" unresolved (boundary)"<<endl;
	}
	double t_extract0=now_s();
	cout<<endl<<"--- surface extraction ---"<<endl;
	int npoly_out[8]={0,0,0,0,0,0,0,0};
	int ntouching_tets=0, nbad_tets=0, ok_cells=0, bad_cells=0, nclipped=0;
	long n_cells_visited=0;

	// alpha mode output (globally deduped)
	vector<pt3> vert_pts;
	map<node_key,int> vert_index;
	vector<vector<int> > faces_out;
	// flat mode output (independent per-polygon vertices, already projected)
	struct vec3d { double v[3]; };
	vector<vec3d> flat_verts;
	vector<vector<int> > flat_faces;

	// --certify (glpt_certify.hpp): certified per-cell tests, tabulated
	// against the heuristic extraction's own good/bad verdict.
	struct CertStats {
		long arc[2][3];      // [extraction ok?][cert_result] for cells with an arc
		long noarc[3];       // [cert_result] for leaves without an arc
		long chartmix[2][2]; // [extraction ok?][some face solved in another chart?]
		std::vector<double> margin[2]; // GRAPH-test margin, [extraction ok?]
		void clear() { for(int a=0;a<2;++a) for(int b=0;b<3;++b) arc[a][b]=0; for(int b=0;b<3;++b) noarc[b]=0; for(int a=0;a<2;++a) for(int b=0;b<2;++b) chartmix[a][b]=0; margin[0].clear(); margin[1].clear(); }
	};
	CertStats cstats; cstats.clear();
	// --cause-stats: why bad cells fail, and the transversality score of
	// good vs bad cells.
	struct CauseStats {
		long reason[4]; std::vector<double> tv[2]; // tv[extraction ok?]
		std::vector<double> diam; // FS diameter of every cell with an arc
		void clear() { for(int i=0;i<4;++i) reason[i]=0; tv[0].clear(); tv[1].clear(); diam.clear(); }
	};
	CauseStats causes; causes.clear();
	struct Extractor {
		CauseStats* cz;
		CertStats* cs;
		const vector<pt3>* gp; glpt_edge_cache* idc; face_crossing_cache* fc;
		long* n_cells; int *ok, *bad, *ntouch, *nbad_t, *npoly, *nclip;
		bool alpha_mode; double cutoff;
		vector<pt3>* vpts; map<node_key,int>* vidx; vector<vector<int> >* faces;
		vector<vec3d>* fverts; vector<vector<int> >* ffaces;
		vector<glpt>* bad_list; // --repair-rounds: which leaves failed to extract, for possible re-refinement

		int emit(const crossing_node& nd) const {
			node_key k = key_of(nd);
			map<node_key,int>::iterator it = vidx->find(k);
			if(it!=vidx->end()) return it->second;
			int idx=(int)vpts->size();
			vpts->push_back(nd.p);
			(*vidx)[k]=idx;
			return idx;
		}
		void operator()(const glpt& c) const {
			++(*n_cells);
			pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
			cell_points_and_ids(c, *gp, *idc, pts, ids);
			cell_extraction_result res;
			extract_cell(pts, ids, *fc, res);
			if(g_certify) {
				double mg; cert_result cr=certify_cell(pts, g_certify_level, mg);
				if(!res.any_edge) ++cs->noarc[cr];
				else {
					int okx=res.decompose_ok?1:0; ++cs->arc[okx][cr]; cs->margin[okx].push_back(mg);
					{ // chart consistency: faces whose best chart differs from the cell's
						int cc=pick_chart5(pts), fdiff=0;
						for(int f=0;f<10;++f) { pt3 fp[3]; for(int k=0;k<3;++k) fp[k]=pts[GLPT_CELL_FACES[f][k]]; if(pick_chart(fp)!=cc) ++fdiff; }
						++cs->chartmix[okx][fdiff>0?1:0];
					}
				}
			}
			if(!res.any_edge) return; // curve doesn't cross this cell at all -- not a failure, nothing to count
			if(g_cause_stats) { ++cz->reason[res.fail_reason]; cz->tv[res.decompose_ok?1:0].push_back(res.min_transversality); cz->diam.push_back(cell_fs_diam(pts)); }
			*ntouch += res.ntouching_tets;
			*nbad_t += res.nbad_tets;
			// Orient every polygon by the complex orientation: reverse the cycle
			// when omega/area < 0 (see polygon_kahler_ratio). Output faces then
			// agree with the curve's orientation, whichever way extract_cell()
			// happened to walk the cycle.
			if(res.decompose_ok)
				for(size_t p=0;p<res.cycles.size();++p) {
					double r=polygon_kahler_ratio(res.nodes,res.cycles[p]);
					if(g_kahler_filter>0 && std::fabs(r)<g_kahler_filter) { res.decompose_ok=false; ++g_kahler_rejected; break; }
					if(r<0) std::reverse(res.cycles[p].begin()+1,res.cycles[p].end()); // keeps cyc[0], the fan apex
				}
			if(!res.decompose_ok) { ++(*bad); if(bad_list) bad_list->push_back(c); return; }
			++(*ok);
			for(size_t p=0;p<res.cycles.size();++p) {
				const vector<int>& cyc = res.cycles[p];
				int sz=(int)cyc.size();
				if(sz>=3 && sz<8) ++npoly[sz];
				if(alpha_mode) {
					vector<int> face;
					for(size_t i=0;i<cyc.size();++i) face.push_back(emit(res.nodes[cyc[i]]));
					faces->push_back(face);
				} else {
					if(g_genus_check) { // combinatorial complex in CP^2 (genus check, --close)
						vector<int> face;
						for(size_t i=0;i<cyc.size();++i) face.push_back(emit(res.nodes[cyc[i]]));
						faces->push_back(face);
					}
					int chart = pick_chart_polygon(res.nodes, cyc);
					vector<vec3d> proj(cyc.size());
					bool clip=false;
					for(size_t i=0;i<cyc.size();++i) {
						double q[3]; flat_projection(res.nodes[cyc[i]].p, chart, q);
						proj[i].v[0]=q[0]; proj[i].v[1]=q[1]; proj[i].v[2]=q[2];
						double dd=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
						if(dd>cutoff) clip=true;
					}
					if(clip) { ++(*nclip); continue; }
					int base=(int)fverts->size();
					vector<int> face(cyc.size());
					for(size_t i=0;i<proj.size();++i) { fverts->push_back(proj[i]); face[i]=base+(int)i; }
					ffaces->push_back(face);
				}
			}
		}
	};
	vector<glpt> bad_list;
	Extractor ext;
	ext.cs=&cstats; ext.cz=&causes;
	ext.gp=&gp; ext.idc=&id_cache; ext.fc=&fcache; ext.n_cells=&n_cells_visited;
	ext.ok=&ok_cells; ext.bad=&bad_cells; ext.ntouch=&ntouching_tets; ext.nbad_t=&nbad_tets; ext.npoly=npoly_out;
	ext.nclip=&nclipped; ext.alpha_mode=alpha_mode; ext.cutoff=cutoff;
	ext.vpts=&vert_pts; ext.vidx=&vert_index; ext.faces=&faces_out;
	ext.fverts=&flat_verts; ext.ffaces=&flat_faces; ext.bad_list=&bad_list;
	tree.for_each_leaf(ext);

	// --repair-rounds: an experiment (user's own idea, revisited --
	// "já tentamos isso antes, mas vale tentar de novo" 2026-09-24):
	// bad cells are never specifically targeted for MORE refinement --
	// Phase 1 only refines by cell_priority/proximity/threshold, so a
	// cell whose LOCAL crossing pattern doesn't fit extract_cell()'s
	// pairing assumptions (near-tangency, multiple close crossings on
	// one facet, ...) just stays bad forever, however deep its
	// neighbors go. This forces exactly those cells (and only those) to
	// bisect further -- up to max_depth+repair_extra_depth, a separate,
	// deliberately generous budget, since a cell already at max_depth
	// after Phase 1 needs genuine EXTRA room to have any chance of
	// resolving -- then re-extracts the WHOLE mesh fresh (simplest
	// correct approach; not the cheapest, but repair mode is opt-in and
	// meant for exploring whether this helps at all, not production
	// use). Known NOT to be a complete fix by itself: some bad
	// configurations (genuine tangency, a crossing sitting exactly on a
	// shared lower-dimensional stratum) don't resolve at ANY depth --
	// see glpt_extraction.hpp's own accepted-degeneracy discussion and
	// [[riemann_glpt_pointerless_mesh]] memory.
	for(int round=0; round<repair_rounds && !bad_list.empty(); ++round) {
		int nrepaired=0;
		for(size_t i=0;i<bad_list.size();++i) {
			glpt c = bad_list[i];
			if(!tree.exists(c)) continue; // already refined via a different path, or superseded
			if(c.simplex_level()>=max_depth+repair_extra_depth) continue;
			tree.compat_bisect(c);
			++nrepaired;
		}
		cout<<"repair round "<<(round+1)<<": "<<bad_list.size()<<" bad cells, "<<nrepaired<<" bisected further"<<endl;
		if(nrepaired==0) break; // every bad cell already maxed out -- no point re-extracting

		n_cells_visited=0; ok_cells=0; bad_cells=0; ntouching_tets=0; nbad_tets=0; nclipped=0;
		cstats.clear(); causes.clear();
		for(int sz=0;sz<8;++sz) npoly_out[sz]=0;
		vert_pts.clear(); vert_index.clear(); faces_out.clear();
		flat_verts.clear(); flat_faces.clear();
		bad_list.clear();
		tree.for_each_leaf(ext);
	}
	if(repair_rounds>0) {
		cout<<"after repair: bad cells remaining: "<<bad_cells<<endl;
		cout<<"final leaf count: "<<tree.leaf_count()<<" (after repair rounds)"<<endl;
	}

	cout<<"cells visited: "<<n_cells_visited<<", distinct faces solved: "<<fcache.size()<<endl;
	if(g_bernstein) cout<<"faces Bernstein-pruned (Newton skipped, certified no root): "<<g_bernstein_pruned_faces<<endl;
	cout<<"extracted cells: ok="<<ok_cells<<" bad="<<bad_cells<<endl;
	if(g_kahler_filter>0) cout<<"cells rejected by --kahler-filter (all rounds): "<<g_kahler_rejected<<endl;
	cout<<"polygon sizes:";
	for(int sz=3;sz<8;++sz) if(npoly_out[sz]) cout<<" "<<sz<<"-gon="<<npoly_out[sz];
	cout<<endl;
	cout<<"touching tetrahedra: "<<ntouching_tets<<"  unhandled node count: "<<nbad_tets<<endl;
	if(g_cause_stats) {
		cout<<endl<<"--- cause stats ---"<<endl;
		if(!causes.diam.empty()) {
			vector<double>& d=causes.diam; sort(d.begin(),d.end());
			cout<<"FS diameter of cells meeting C: median="<<d[d.size()/2]<<" p90="<<d[(9*d.size())/10]<<" max="<<d.back()<<endl;
		}
		cout<<"bad cells: irregular facet="<<causes.reason[1]<<" node degree="<<causes.reason[2]<<" bad cycle="<<causes.reason[3]<<endl;
		for(int okx=1;okx>=0;--okx) {
			vector<double>& m=causes.tv[okx];
			if(m.empty()) continue;
			sort(m.begin(),m.end());
			long low=0; for(double x: m) if(x<0.05) ++low;
			cout<<(okx?"good":"bad ")<<" cells: transversality p10="<<m[m.size()/10]<<" median="<<m[m.size()/2]
				<<" fraction<0.05="<<(double)low/m.size()<<endl;
		}
	}
	if(g_certify) {
		const char* nm[3]={"EMPTY","GRAPH","UNDECIDED"};
		cout<<endl<<"--- certify (level "<<g_certify_level<<") ---"<<endl;
		cout<<"leaves without an arc:";
		for(int b=0;b<3;++b) cout<<" "<<nm[b]<<"="<<cstats.noarc[b];
		cout<<endl;
		for(int mix=0;mix<2;++mix) {
			long g=cstats.chartmix[1][mix], bd=cstats.chartmix[0][mix];
			cout<<"cells with arc, "<<(mix?"some face in another chart":"all faces in cell chart")<<": good="<<g<<" bad="<<bd
				<<" bad rate="<<(g+bd?100.0*bd/(g+bd):0.0)<<"%"<<endl;
		}
		for(int okx=1;okx>=0;--okx) {
			long tot=0; for(int b=0;b<3;++b) tot+=cstats.arc[okx][b];
			cout<<(okx?"good":"bad ")<<" cells (with arc), "<<tot<<":";
			for(int b=0;b<3;++b) cout<<" "<<nm[b]<<"="<<cstats.arc[okx][b]<<" ("<<(tot?100.0*cstats.arc[okx][b]/tot:0.0)<<"%)";
			cout<<endl;
			vector<double>& m=cstats.margin[okx];
			if(!m.empty()) {
				sort(m.begin(),m.end());
				cout<<"   GRAPH margin (gap-pi, rad): p10="<<m[m.size()/10]<<" median="<<m[m.size()/2]<<" p90="<<m[(9*m.size())/10]<<endl;
			}
		}
	}
	if(!alpha_mode && cutoff<1e299) cout<<"polygons dropped by --cutoff: "<<nclipped<<endl;
	if(g_bernstein) cout<<"Bernstein fallback seeds used: "<<g_bfallback_tried<<" faces tried, "
		<<g_bfallback_new_root<<" found a root the fixed/dynamic seeds missed"<<endl;

	if(g_genus_check && getenv("GLPT_NONMANIFOLD_DEBUG")) { // which nodes form a non-manifold edge
		std::map<std::pair<int,int>,std::vector<int> > ef;
		for(size_t f=0;f<faces_out.size();++f) { int n=(int)faces_out[f].size();
			for(int i=0;i<n;++i) { int a=faces_out[f][i], b=faces_out[f][(i+1)%n]; ef[std::make_pair(std::min(a,b),std::max(a,b))].push_back((int)f); } }
		std::vector<node_key> inv(vert_pts.size());
		for(auto& kv: vert_index) inv[kv.second]=kv.first;
		for(auto& e: ef) if(e.second.size()>2) {
			cout<<"non-manifold edge "<<e.first.first<<"-"<<e.first.second<<" in "<<e.second.size()<<" faces, FS length "
				<<fs_dist(vert_pts[e.first.first],vert_pts[e.first.second])<<endl;
			for(int v: {e.first.first,e.first.second}) { const node_key& k=inv[v];
				cout<<"  node "<<v<<": dim="<<k.dim<<" a="<<k.a<<" b="<<k.b<<" q="<<k.qparam<<endl; }
			for(int f: e.second) { cout<<"  face "<<f<<":"; for(int x: faces_out[f]) cout<<" "<<x; cout<<endl; }
		}
		int worst=-1; double wr=-1; // vertex with the largest |F|
		for(size_t v=0;v<vert_pts.size();++v) { pt3 u=normalize3(vert_pts[v]); double r=std::abs(eval_poly3(g_F,u[0],u[1],u[2])); if(r>wr) { wr=r; worst=(int)v; } }
		if(worst>=0) { const node_key& k=inv[worst];
			cout<<"largest |F| = "<<wr<<" at node "<<worst<<": dim="<<k.dim<<" a="<<k.a<<" b="<<k.b<<" q="<<k.qparam<<endl; }
	}
	if(g_genus_check) {
		int deg=g_F.d;
		genus_check(faces_out, (int)vert_pts.size(), (deg-1)*(deg-2)/2);
		surface_stats(vert_pts, faces_out, deg);
		if(g_close) {
			close_surface(vert_pts, faces_out);
			genus_check(faces_out, (int)vert_pts.size(), (deg-1)*(deg-2)/2);
		surface_stats(vert_pts, faces_out, deg);
			if(!alpha_mode) { // re-project the closed complex (per polygon, as in flat mode)
				flat_verts.clear(); flat_faces.clear(); nclipped=0;
				for(auto& f: faces_out) {
					int chart=g_forced_chart;
					if(chart<0) { double bestm=-1; for(int c=0;c<3;++c) { double m=1e300; for(int x: f) m=std::min(m,std::abs(vert_pts[x][c])); if(m>bestm) { bestm=m; chart=c; } } }
					std::vector<int> face; bool clip=false; int base=(int)flat_verts.size();
					for(int x: f) { double q[3]; flat_projection(vert_pts[x],chart,q); vec3d w; w.v[0]=q[0]; w.v[1]=q[1]; w.v[2]=q[2];
						if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2])>cutoff) clip=true; flat_verts.push_back(w); face.push_back((int)flat_verts.size()-1); }
					if(clip) { flat_verts.resize(base); ++nclipped; continue; }
					flat_faces.push_back(face);
				}
			}
		}
	}
	if(g_leaf_hash) {
		struct H { uint64_t* h; long* n; void operator()(const glpt& c) const {
			uint64_t z=c.raw()+0x9e3779b97f4a7c15ULL; z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; z^=z>>31;
			*h+=z; ++*n; } };
		uint64_t h=0; long n=0; H f; f.h=&h; f.n=&n; tree.for_each_leaf(f);
		cout<<"leaf set: "<<n<<" leaves, hash "<<std::hex<<h<<std::dec<<endl;
	}
	if(g_timing) {
		double t_end=now_s();
		cout<<"timing: refinement "<<(t_extract0-t_refine0)<<" s, extraction (incl. repair) + post-processing "
			<<(t_end-t_extract0)<<" s, peak RSS "<<peak_rss_kb()/1024.0<<" MB"<<endl;
	}
	// --- Phase 3: output (OBJ mesh) --
	ofstream out(obj_path.c_str());
	if(!out) { cerr<<"couldn't open "<<obj_path<<" for writing"<<endl; return 1; }
	out<<"# riemann_cp2_glpt surface extraction: curve="<<cat[function_idx].name
		<<" depth="<<max_depth<<" threshold="<<threshold<<"\n";
	if(alpha_mode) {
		out<<"# projected via alpha-tilde (Dutter arXiv:2608.04323, after Kranich 2015)\n";
		for(size_t i=0;i<vert_pts.size();++i) {
			double p[3]; alpha_projection(vert_pts[i], p);
			out<<"v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<"\n";
		}
		for(size_t i=0;i<faces_out.size();++i) {
			out<<"f";
			for(size_t k=0;k<faces_out[i].size();++k) out<<" "<<(faces_out[i][k]+1);
			out<<"\n";
		}
		cout<<"wrote "<<obj_path<<": "<<vert_pts.size()<<" vertices, "<<faces_out.size()<<" faces"<<endl;
	} else {
		out<<"# projected via per-polygon best affine chart (Re a, Im a, Re b)\n";
		if(cutoff<1e299) out<<"# cutoff: polygons with a vertex farther than "<<cutoff<<" from the origin dropped\n";
		for(size_t i=0;i<flat_verts.size();++i)
			out<<"v "<<flat_verts[i].v[0]<<" "<<flat_verts[i].v[1]<<" "<<flat_verts[i].v[2]<<"\n";
		for(size_t i=0;i<flat_faces.size();++i) {
			out<<"f";
			for(size_t k=0;k<flat_faces[i].size();++k) out<<" "<<(flat_faces[i][k]+1);
			out<<"\n";
		}
		cout<<"wrote "<<obj_path<<": "<<flat_verts.size()<<" vertices, "<<flat_faces.size()<<" faces"<<endl;
	}
	return 0;
}
