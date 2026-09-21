// slice_probe -- a controlled stress test for the codim-2 crossing-graph
// extraction algorithm used by examples/top/riemann/riemann.cpp and
// examples/top/riemann_cp2/riemann_cp2.cpp, isolated from everything
// algebraic (Newton iteration, poles, branch points, ...).
//
// Setup: take the standard 4-simplex in R^4 (vertices 0, e1, e2, e3, e4).
// Draw a uniformly random interior point P and two random unit
// directions u, v. The affine subspace
//
//   L = { x : dot(x-P,u)=0, dot(x-P,v)=0 }
//
// is codimension 2 in R^4 and passes through P, which is interior to
// the simplex, so L always meets it. Maubach-refine the single seed
//4-simplex (biased toward the region L actually crosses -- a filter
// that is EXACT here, not a heuristic proxy, since the field is
// literally affine: an affine map's extrema over a simplex are at its
// vertices), then run the *same* per-tetrahedron/per-4-simplex
// crossing-graph closure used by riemann.cpp's Phase 3 (find crossing
// nodes on 2-faces, canonically identify them, connect them into
// edges via each of a cell's 5 tetrahedral facets, and check the
// result closes into cycles).
//
// Why this is a strong test, not just a simpler one: every individual
// 4-simplex cell is convex, so its intersection with the affine plane
// L is provably EITHER empty, a point, OR a single convex polygon --
// never more than one polygon, never an open (non-cycle) fan. So here,
// unlike the algebraic-curve case, a "bad cell" (crossing graph doesn't
// close) or a cell producing more than one cycle is not a possible
// geometric outcome -- it can only be a bug or a tolerance/precision
// issue in the extraction code itself (canonical node identity,
// degree/cycle bookkeeping, the barycentric "on this face" tolerance).
// Repeating this over thousands of random (P,u,v) draws, at increasing
// refinement depth, is a clean-room way to find those cases.
//
// Per triangle (2-simplex), the crossing (if any) is found by an EXACT
// 2x2 linear solve -- not Newton -- since the field is affine on the
// triangle: this is the same closed-form "linear seed" riemann.cpp
// computes as a Newton starting guess (triangle_intersection there),
// except here it IS the answer, not just a seed.

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/euler.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/alg/vec.hpp>

#define DIM 4

using namespace std;
using namespace vgtl;

typedef vgtl::nmt<DIM> T;
typedef vec<DIM,double> pt4;

// --- the probe: a random point P and two random unit directions u,v ----

pt4 g_P, g_u, g_v;

struct field2 { double f1, f2; };

field2 field_at(const pt4& x) {
	pt4 d=x-g_P;
	field2 r; r.f1=dot(d,g_u); r.f2=dot(d,g_v);
	return r;
}

namespace vgtl {

	template <>
	struct extra_data<0> {
		pt4 pos;
	};

	// Cache of the (at most one) L-crossing of a 2-simplex. Unlike
	// riemann.cpp's extra_data<2> (up to 2 roots, since F has degree n
	// in w), a codim-2 AFFINE field crosses a 2-simplex in at most one
	// point generically -- exactly one linear system, exactly one
	// answer (or none, if it falls outside the triangle, or the system
	// is degenerate).
	template <>
	struct extra_data<2> {
		unsigned char computed, has_root;
		double l1,l2; // barycentric coords of the crossing (l0=1-l1-l2)
		extra_data() : computed(0), has_root(0), l1(0), l2(0) {}
	};

}

pt4 pos(const T& t, Vertex(T) v) { return attr(t,v)->pos; }
void pos_set(T& t, Vertex(T) v, const pt4& p) { attr(t,v)->pos=p; }

// --- refinement: bisect only cells L could possibly cross ---------------

// EXACT necessary condition: f1 and f2 are each affine over the cell, so
// each attains its extrema at the cell's own vertices. If 0 isn't
// bracketed by a vertex's min/max in BOTH f1 and f2, no point of the
// cell (let alone of L) lies in it. This never discards a cell that
// could actually contain a crossing (it's a superset test, just not
// tight against the joint (f1,f2) convex hull) -- a scheduling filter
// only, extraction correctness below never depends on it.
bool cell_relevant(const T& t, Cell(T) cv) {
	vgtl::array<Vertex(T),DIM+1> vs;
	vertices(t,cv,vs);
	double f1lo=1e300,f1hi=-1e300,f2lo=1e300,f2hi=-1e300;
	for(int i=0;i<=DIM;++i) {
		field2 f=field_at(pos(t,vs[i]));
		if(f.f1<f1lo) f1lo=f.f1; if(f.f1>f1hi) f1hi=f.f1;
		if(f.f2<f2lo) f2lo=f.f2; if(f.f2>f2hi) f2hi=f.f2;
	}
	return f1lo<=0 && f1hi>=0 && f2lo<=0 && f2hi>=0;
}

struct refine_app : do_nothing {
	using do_nothing::apply;
	vector<Cell(T)> new_cells;

	Vertex(T) add_edge_vertex(T& t, Edge(T) e) {
		vgtl::array<Vertex(T),2> vs;
		vertices(t,e,vs);
		Vertex(T) v=add(t);
		pos_set(t,v, 0.5*(pos(t,vs[0])+pos(t,vs[1])) );
		return v;
	}

	void apply(T& t, Cell(T) s) { new_cells.push_back(s); }
};

// --- Phase 3 analogue: per-triangle exact linear crossing --------------

long long g_degenerate_tris=0;   // rank-deficient (f1,f2) image on this triangle
long long g_tri_in_plane=0;      // degenerate AND field ~0 at all 3 vertices (whole triangle in L)

void triangle_intersection(T& t, Simplex(T,2) tri) {
	extra_data<2>* d=attr(t,tri);
	if(d->computed) return;
	d->computed=true;

	vgtl::array<Vertex(T),3> vs;
	vertices(t,tri,vs);
	pt4 X0=pos(t,vs[0]), X1=pos(t,vs[1]), X2=pos(t,vs[2]);
	field2 F0=field_at(X0), F1=field_at(X1), F2=field_at(X2);

	double a1f1=F1.f1-F0.f1, a2f1=F2.f1-F0.f1;
	double a1f2=F1.f2-F0.f2, a2f2=F2.f2-F0.f2;
	double na1=sqrt(a1f1*a1f1+a1f2*a1f2), na2=sqrt(a2f1*a2f1+a2f2*a2f2);
	double det=a1f1*a2f2-a2f1*a1f2;
	// Scale-invariant degeneracy test: det/(na1*na2) is sin(angle) between
	// the two edge images in field space, so this stays meaningful however
	// small the triangle has become under refinement (a fixed absolute
	// threshold on det itself would not).
	double reldet = (na1>1e-300 && na2>1e-300) ? det/(na1*na2) : 0.0;
	if(fabs(reldet)<1e-9) {
		++g_degenerate_tris;
		double m=max(max(fabs(F0.f1),fabs(F0.f2)),max(max(fabs(F1.f1),fabs(F1.f2)),max(fabs(F2.f1),fabs(F2.f2))));
		if(m<1e-9) ++g_tri_in_plane;
		return;
	}

	double l1=(-F0.f1*a2f2+a2f1*F0.f2)/det;
	double l2=(-a1f1*F0.f2+F0.f1*a1f2)/det;
	double l0=1.0-l1-l2;
	const double dom_tol=1e-9;
	if(l0<-dom_tol || l1<-dom_tol || l2<-dom_tol) return; // crossing outside this triangle

	d->has_root=1; d->l1=l1; d->l2=l2;
}

bool triangle_has_root(const T& t, Simplex(T,2) tri) { return attr(t,tri)->has_root!=0; }

void triangle_bary(const T& t, Simplex(T,2) tri, double& l0, double& l1, double& l2) {
	l1=attr(t,tri)->l1; l2=attr(t,tri)->l2; l0=1.0-l1-l2;
}

// EXPERIMENT 1 (reverted): tried converting the fixed 1e-6 barycentric
// tolerance into a real-space-absolute one (loosening it in ill-shaped
// cells), both via the triangle's own diameter and via each specific
// edge's own length. Both made --corner-bias MEASURABLY WORSE (62 ->
// 479, then 514, bad cells /500 trials), not better.
//
// The per-tetrahedron trace this file's failure dumps now include
// explains why, and identifies the REAL bug -- the opposite of what a
// looser tolerance would fix. Example (fail_t22_c58): two DIFFERENT
// triangles (desc 347 and 377), each with its OWN genuinely distinct
// crossing point (~4e-7 apart -- far above float noise, and comparable
// to how far each point actually sits off the shared edge, i.e. two
// real, separate near-edge points, not one point found twice), both
// independently land within the SAME 1e-6 window of shared edge 147 and
// get MERGED into one node by same_crossing_node. That collapses what
// should be a proper 3-node cycle (0 -> node_A -> node_B -> 0, a
// triangle) into a spurious 0<->1 digon: the tetrahedron that should
// contribute the (node_A, node_B) edge instead sees them as the SAME
// node and reports ntn=1 ("touching", no edge at all).
//
// So the classification tolerance (is this point close enough to a
// face to special-case it) and the MERGE tolerance (are two
// already-classified points actually the same one) are answering two
// different questions and were sharing one constant by accident. A
// genuinely-identical point found from two different incident
// triangles agrees far more tightly than 1e-6 (their formula reduces to
// the same edge-endpoint-only computation when the point is truly ON
// the edge, so the only disagreement left is ordinary floating-point
// roundoff) -- so MERGE_TOL can be, and needs to be, far tighter than
// SNAP_TOL to stop conflating "close to a shared edge" with "the same
// point on it".
const double SNAP_TOL=1e-6;   // classification: is this point on a lower-dim face?
const double MERGE_TOL=1e-9;  // identity: are two same-face points the same point?

// Canonical identity of a crossing point, exactly as riemann.cpp's
// compute_crossing_node: which mesh element (vertex/edge/triangle
// interior) it lands on, in barycentric (scale-invariant) terms.
struct crossing_node {
	int dim, desc;
	double param;
	pt4 x;
};

void compute_crossing_node(const T& t, Simplex(T,2) tri, crossing_node& nd) {
	double l0,l1,l2; triangle_bary(t,tri,l0,l1,l2);
	double l[3]={l0,l1,l2};
	int zeros=0, zi[3];
	for(int b=0;b<3;++b) if(fabs(l[b])<SNAP_TOL) zi[zeros++]=b;
	vgtl::array<Vertex(T),3> vs; vertices(t,tri,vs);
	nd.x = l0*pos(t,vs[0]) + l1*pos(t,vs[1]) + l2*pos(t,vs[2]);
	nd.param=0;
	if(zeros>=2) {
		int keep=(zeros==3) ? 0 : (3-zi[0]-zi[1]);
		nd.dim=0; nd.desc=vs[keep].desc;
	} else if(zeros==1) {
		int lo=(zi[0]==0)?1:0, hi=(zi[0]==2)?1:2;
		nd.dim=1; nd.desc=face_op(t,tri,zi[0]).desc;
		nd.param=l[hi]/(l[lo]+l[hi]);
	} else {
		nd.dim=2; nd.desc=tri.desc;
	}
}

bool same_crossing_node(const crossing_node& a, const crossing_node& b) {
	if(a.dim!=b.dim || a.desc!=b.desc) return false;
	if(a.dim==1 && fabs(a.param-b.param)>MERGE_TOL) return false;
	return true;
}

// --- orthonormal frame for L, used only to write debug wireframes ------

struct frame4 { pt4 e1,e2,e3,e4; };

frame4 build_frame() {
	frame4 fr;
	fr.e1=g_u;
	pt4 w=g_v-dot(g_v,fr.e1)*fr.e1;
	double n=sqrt(dot(w,w));
	fr.e2 = (n>1e-9) ? (1.0/n)*w : pt4();
	pt4 cand[2]; int nc=0;
	for(int i=0;i<4 && nc<2;++i) {
		pt4 b; b[i]=1.0;
		b -= dot(b,fr.e1)*fr.e1;
		b -= dot(b,fr.e2)*fr.e2;
		for(int k=0;k<nc;++k) b -= dot(b,cand[k])*cand[k];
		double bn=sqrt(dot(b,b));
		if(bn>1e-6) { cand[nc]=(1.0/bn)*b; ++nc; }
	}
	fr.e3 = (nc>0)?cand[0]:pt4();
	fr.e4 = (nc>1)?cand[1]:pt4();
	return fr;
}

// project to (e2,e3,e4): a crossing node genuinely on L has e1=e2=0, so
// this view still shows a nonzero "e2 residual" for anything that
// ISN'T really on L -- useful to spot a canonicalization bug visually.
void proj3(const frame4& fr, const pt4& x, double out[3]) {
	pt4 d=x-g_P;
	out[0]=dot(d,fr.e2); out[1]=dot(d,fr.e3); out[2]=dot(d,fr.e4);
}

// --- per-trial run --------------------------------------------------

struct trial_result {
	int ok_cells, bad_cells, multi_cycle_cells;
	int touching_tets, unhandled_tets;
	long long n_subdivisions, n_leaf_cells;
};

pt4 sample_unit(mt19937& rng) {
	normal_distribution<double> nd(0.0,1.0);
	pt4 r; double n2=0;
	for(int i=0;i<DIM;++i) { r[i]=nd(rng); n2+=r[i]*r[i]; }
	double n=sqrt(n2);
	for(int i=0;i<DIM;++i) r[i]/=n;
	return r;
}

void write_obj_wireframe(const string& path, const frame4& fr,
                          const pt4 V[DIM+1], const vector<crossing_node>& nodes,
                          const vector<vector<int> >& adj) {
	ofstream o(path.c_str());
	o<<"# slice_probe debug wireframe: seed 4-simplex edges + crossing-node graph\n";
	o<<"# axes are (e2,e3,e4) of the L-adapted orthonormal frame: a node truly\n";
	o<<"# on L has zero e1,e2 residual by construction, so any visible gap\n";
	o<<"# between two edges meant to meet at 'the same' node is real.\n";
	vector<double> pts; // flattened xyz, in write order
	int base=1;
	for(int i=0;i<=DIM;++i) {
		double p[3]; proj3(fr,V[i],p);
		o<<"v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" # simplex vertex "<<i<<"\n";
	}
	for(int i=0;i<=DIM;++i)
		for(int j=i+1;j<=DIM;++j)
			o<<"l "<<(base+i)<<" "<<(base+j)<<"\n";
	base+=DIM+1;
	for(size_t q=0;q<nodes.size();++q) {
		double p[3]; proj3(fr,nodes[q].x,p);
		o<<"v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" # crossing node "<<q
			<<" dim="<<nodes[q].dim<<" deg="<<adj[q].size()<<"\n";
	}
	for(size_t q=0;q<nodes.size();++q)
		for(size_t k=0;k<adj[q].size();++k)
			if(adj[q][k]>(int)q) o<<"l "<<(base+q)<<" "<<(base+adj[q][k])<<"\n";
}

// Adversarial knobs: the pure-random draw below essentially never hits
// a floating-point-relevant degenerate alignment (measure zero, and
// numerically the shrinking cell size keeps crossings well clear of
// exact ties at any depth we can afford -- confirmed empirically: zero
// bad cells over ~700k extracted cells, depths 8-20). These bias the
// sampling toward the kinds of alignment that, in riemann.cpp, DID
// cause real failures (a probe landing exactly on/near a low-level
// mesh feature that Maubach bisection never perturbs).
struct probe_config {
	double corner_bias;   // 0=off; else P is pulled this fraction of the way
	                      // toward a random ORIGINAL simplex vertex (the
	                      // codim-2 analogue of riemann.cpp's branch point
	                      // sitting exactly on a seed vertex forever)
	double parallel_eps;  // 0=off; else v is forced within this many radians
	                      // of u, stressing the per-triangle rank test
	double axis_align_p;  // probability u (independently, v) is replaced by
	                      // a coordinate axis or a simplex edge direction,
	                      // instead of a fully generic random unit vector
	probe_config() : corner_bias(0), parallel_eps(0), axis_align_p(0) {}
};

pt4 maybe_axis_align(mt19937& rng, const pt4& fallback, double p, const pt4 V[DIM+1]) {
	uniform_real_distribution<double> U01(0.0,1.0);
	if(U01(rng)>=p) return fallback;
	uniform_int_distribution<int> pick(0,2*DIM-1);
	int k=pick(rng);
	if(k<DIM) { pt4 r; r[k]=1.0; return r; } // coordinate axis
	// a simplex edge direction (Vj-Vi), normalized
	int i=k-DIM, j=(i+1)%(DIM+1);
	pt4 d=V[j]-V[i];
	double n=sqrt(dot(d,d));
	return n>1e-12 ? (1.0/n)*d : fallback;
}

trial_result run_trial(unsigned seed, int max_depth, int trial_idx,
                        const string& dump_dir, int& dump_budget,
                        int ok_by_level[32], int bad_by_level[32],
                        const probe_config& cfg) {
	mt19937 rng(seed);
	uniform_real_distribution<double> U(0.0,1.0);

	// standard 4-simplex: 0, e1, e2, e3, e4
	pt4 V[DIM+1];
	for(int i=1;i<=DIM;++i) V[i][i-1]=1.0;

	double b[DIM+1], bsum=0;
	for(int i=0;i<=DIM;++i) { double u=U(rng); b[i]=-log(1.0-u); bsum+=b[i]; }
	for(int i=0;i<=DIM;++i) b[i]/=bsum;
	pt4 P; for(int i=0;i<=DIM;++i) P += b[i]*V[i];

	if(cfg.corner_bias>0) {
		uniform_int_distribution<int> pickv(0,DIM);
		int cvx=pickv(rng);
		P = (1.0-cfg.corner_bias)*P + cfg.corner_bias*V[cvx];
	}

	pt4 u=sample_unit(rng), v=sample_unit(rng);
	for(int tries=0; tries<20 && fabs(dot(u,v))>0.999; ++tries) v=sample_unit(rng);
	u=maybe_axis_align(rng,u,cfg.axis_align_p,V);
	v=maybe_axis_align(rng,v,cfg.axis_align_p,V);

	if(cfg.parallel_eps>0) {
		// v = cos(eps)*u + sin(eps)*(component of the random draw
		// orthogonal to u), i.e. exactly eps radians from u -- forces a
		// small, controlled, nonzero angle instead of a fully random one.
		pt4 w=v-dot(v,u)*u;
		double wn=sqrt(dot(w,w));
		if(wn>1e-12) {
			w=(1.0/wn)*w;
			v = cos(cfg.parallel_eps)*u + sin(cfg.parallel_eps)*w;
		}
	}

	g_P=P; g_u=u; g_v=v;

	T t;
	vgtl::array<Vertex(T),DIM+1> vs;
	for(int i=0;i<=DIM;++i) { vs[i]=add(t); pos_set(t,vs[i],V[i]); }
	Cell(T) seed_cell;
	{
		vgtl::complex_builder<T> cb(t);
		seed_cell=add(cb,vs);
	}
	orientation_set(t,seed_cell,1);
	level_set(t,seed_cell,0);

	// --- refine: bisect every cell that could possibly contain a
	// crossing, down to max_depth ---
	long long nsub=0;
	{
		vector<Cell(T)> stack;
		{
			Cell_it(T) ci,cend;
			for(simplices(t,ci,cend); ci!=cend; ++ci)
				if(is_current(t,*ci) && cell_relevant(t,*ci)) stack.push_back(*ci);
		}
		while(!stack.empty()) {
			Cell(T) cv=stack.back(); stack.pop_back();
			if(!is_current(t,cv)) continue;
			if(level(t,cv)>=max_depth) continue;
			refine_app app;
			maubach_subdivide(t,cv,app);
			++nsub;
			for(size_t i=0;i<app.new_cells.size();++i)
				if(cell_relevant(t,app.new_cells[i])) stack.push_back(app.new_cells[i]);
		}
	}

	// --- extraction: identical crossing-graph closure to riemann.cpp's
	// Phase 3, minus the multi-sheet (up to 4 nodes/tet) pairing case,
	// which cannot legitimately arise for a single affine plane: a
	// tetrahedron meeting >2 nodes here is treated as unhandled, exactly
	// the diagnostic we want (see file header).
	trial_result res; memset(&res,0,sizeof(res));
	res.n_subdivisions=nsub;

	frame4 fr; bool fr_built=false;

	Cell_it(T) ci,cend;
	for(simplices(t,ci,cend); ci!=cend; ++ci) {
		Cell(T) sigma=*ci;
		if(!is_current(t,sigma)) continue;
		++res.n_leaf_cells;

		vector<crossing_node> nodes;
		vector<vector<int> > adj;
		bool tet_unhandled=false;

		// Diagnostic-only: full per-tetrahedron, per-triangle-face trace
		// (raw barycentric l0,l1,l2, classification, matched node index),
		// so a failure dump can show not just the final (deduped) node
		// list but WHY each of the cell's 5 tetrahedra ended up with the
		// node count it did -- e.g. whether a facet contributing only 1
		// node ("touching") is genuinely grazing the surface there, or
		// silently lost what should have been its 2nd node to a
		// same_crossing_node match against the wrong existing node.
		ostringstream tet_trace;

		for(int j=0;j<=DIM;++j) {
			Simplex(T,3) tet=face_op(t,sigma,j);
			int tnodes[4]; int ntn=0;
			tet_trace<<"  tet "<<j<<":\n";
			for(int k=0;k<=3;++k) {
				Simplex(T,2) tri=face_op(t,tet,k);
				triangle_intersection(t,tri);
				if(!triangle_has_root(t,tri)) { tet_trace<<"    face "<<k<<": no root\n"; continue; }
				double rl0,rl1,rl2; triangle_bary(t,tri,rl0,rl1,rl2);
				crossing_node nd; compute_crossing_node(t,tri,nd);
				int idx=-1;
				for(size_t q=0;q<nodes.size();++q)
					if(same_crossing_node(nodes[q],nd)) { idx=(int)q; break; }
				if(idx<0) { idx=(int)nodes.size(); nodes.push_back(nd); adj.push_back(vector<int>()); }
				tet_trace<<"    face "<<k<<" tri_desc="<<tri.desc<<": raw l=("<<rl0<<","<<rl1<<","<<rl2<<")"
					<<" -> dim="<<nd.dim<<" desc="<<nd.desc<<" param="<<nd.param<<" -> node "<<idx<<"\n";
				bool dup=false;
				for(int q=0;q<ntn;++q) if(tnodes[q]==idx) dup=true;
				if(!dup && ntn<4) tnodes[ntn++]=idx;
			}
			tet_trace<<"    ntn="<<ntn<<" nodes={"; for(int q=0;q<ntn;++q) tet_trace<<(q?",":"")<<tnodes[q]; tet_trace<<"}\n";
			if(ntn==2) {
				adj[tnodes[0]].push_back(tnodes[1]);
				adj[tnodes[1]].push_back(tnodes[0]);
			} else if(ntn==1) {
				++res.touching_tets;
			} else if(ntn>=3) {
				++res.unhandled_tets; tet_unhandled=true;
			}
		}

		bool deg_ok=true, any_edge=false;
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
					if(seen[cur]) { decompose_ok=false; break; }
					seen[cur]=true; cyc.push_back(cur);
					int nxt=(adj[cur][0]==prev)?adj[cur][1]:adj[cur][0];
					prev=cur; cur=nxt;
					if(cyc.size()>(size_t)(DIM+1)) { decompose_ok=false; break; }
				}
				if(!decompose_ok) break;
				if(cyc.size()<3) { decompose_ok=false; break; }
				cycles_idx.push_back(cyc);
			}
		}

		int lvl=level(t,sigma); if(lvl<0) lvl=0; if(lvl>31) lvl=31;

		if(!decompose_ok) {
			++res.bad_cells; ++bad_by_level[lvl];
			if(dump_budget>0 && !dump_dir.empty()) {
				if(!fr_built) { fr=build_frame(); fr_built=true; }
				ostringstream name; name<<dump_dir<<"/fail_t"<<trial_idx<<"_c"<<res.n_leaf_cells<<".obj";
				write_obj_wireframe(name.str(),fr,V,nodes,adj);
				ostringstream tname; tname<<dump_dir<<"/fail_t"<<trial_idx<<"_c"<<res.n_leaf_cells<<".txt";
				ofstream tf(tname.str().c_str());
				tf<<setprecision(17);
				tf<<"trial="<<trial_idx<<" seed="<<seed<<" depth="<<max_depth<<" cell_level="<<lvl<<"\n";
				tf<<"P="<<P[0]<<" "<<P[1]<<" "<<P[2]<<" "<<P[3]<<"\n";
				tf<<"u="<<u[0]<<" "<<u[1]<<" "<<u[2]<<" "<<u[3]<<"\n";
				tf<<"v="<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<"\n";
				tf<<"deg_ok="<<deg_ok<<" tet_unhandled="<<tet_unhandled<<"\n";
				tf<<"nodes="<<nodes.size()<<"\n";
				for(size_t q=0;q<nodes.size();++q) {
					tf<<"  node "<<q<<": dim="<<nodes[q].dim<<" desc="<<nodes[q].desc
						<<" param="<<nodes[q].param<<" deg="<<adj[q].size()
						<<" x=("<<nodes[q].x[0]<<","<<nodes[q].x[1]<<","<<nodes[q].x[2]<<","<<nodes[q].x[3]<<")"
						<<" adj={"; for(size_t k=0;k<adj[q].size();++k) tf<<(k?",":"")<<adj[q][k]; tf<<"}\n";
				}
				tf<<"per-tetrahedron trace:\n"<<tet_trace.str();
				--dump_budget;
			}
			continue;
		}
		++res.ok_cells; ++ok_by_level[lvl];
		if(cycles_idx.size()>1) ++res.multi_cycle_cells;
	}

	return res;
}

int main(int argc, char* argv[]) {
	int trials=2000;
	int max_depth=8;
	unsigned base_seed=12345;
	string dump_dir;
	int dump_budget=20;
	probe_config cfg;

	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--trials" && i+1<argc) trials=atoi(argv[++i]);
		else if(arg=="--depth" && i+1<argc) max_depth=atoi(argv[++i]);
		else if(arg=="--seed" && i+1<argc) base_seed=(unsigned)atoi(argv[++i]);
		else if(arg=="--dump-dir" && i+1<argc) dump_dir=argv[++i];
		else if(arg=="--dump-budget" && i+1<argc) dump_budget=atoi(argv[++i]);
		else if(arg=="--corner-bias" && i+1<argc) cfg.corner_bias=atof(argv[++i]);
		else if(arg=="--parallel-eps" && i+1<argc) cfg.parallel_eps=atof(argv[++i]);
		else if(arg=="--axis-align-p" && i+1<argc) cfg.axis_align_p=atof(argv[++i]);
		else {
			cerr<<"unrecognized argument: "<<arg<<endl;
			cerr<<"usage: "<<argv[0]<<" [--trials N] [--depth N] [--seed N] "
					<<"[--dump-dir DIR] [--dump-budget N] [--corner-bias X] "
					<<"[--parallel-eps RAD] [--axis-align-p X]"<<endl;
			return 1;
		}
	}
	if(!dump_dir.empty()) {
		string cmd="mkdir -p "+dump_dir;
		system(cmd.c_str());
	}

	cout<<"slice_probe: "<<trials<<" trials, max depth "<<max_depth
			<<", base seed "<<base_seed<<endl;

	long long tot_ok=0, tot_bad=0, tot_multi=0, tot_touch=0, tot_unhandled=0;
	long long tot_leaf=0, tot_sub=0;
	int trials_with_bad=0, trials_with_multi=0;
	int worst_bad=0, worst_bad_trial=-1;
	int ok_by_level[32]={0}, bad_by_level[32]={0};

	for(int trial=0; trial<trials; ++trial) {
		unsigned seed=base_seed+(unsigned)trial*2654435761u; // Knuth multiplicative spread
		trial_result r=run_trial(seed,max_depth,trial,dump_dir,dump_budget,ok_by_level,bad_by_level,cfg);
		tot_ok+=r.ok_cells; tot_bad+=r.bad_cells; tot_multi+=r.multi_cycle_cells;
		tot_touch+=r.touching_tets; tot_unhandled+=r.unhandled_tets;
		tot_leaf+=r.n_leaf_cells; tot_sub+=r.n_subdivisions;
		if(r.bad_cells>0) ++trials_with_bad;
		if(r.multi_cycle_cells>0) ++trials_with_multi;
		if(r.bad_cells>worst_bad) { worst_bad=r.bad_cells; worst_bad_trial=trial; }

		if(trials<=50 || trial%max(1,trials/20)==0)
			cout<<"trial "<<trial<<": leaf_cells="<<r.n_leaf_cells
					<<" ok="<<r.ok_cells<<" bad="<<r.bad_cells
					<<" multi_cycle="<<r.multi_cycle_cells<<endl;
	}

	cout<<endl<<"=== summary over "<<trials<<" trials (depth "<<max_depth<<") ==="<<endl;
	cout<<"total leaf cells: "<<tot_leaf<<" (avg "<<(double)tot_leaf/trials<<"/trial, "
			<<tot_sub<<" total subdivisions)"<<endl;
	cout<<"ok_cells="<<tot_ok<<" bad_cells="<<tot_bad
			<<" ("<<(100.0*tot_bad/max(1LL,tot_ok+tot_bad))<<"%)"<<endl;
	cout<<"multi_cycle_cells="<<tot_multi<<" (geometrically impossible for a convex "
			<<"cell -- any nonzero count here is a confirmed bug, not a tolerance issue)"<<endl;
	cout<<"touching_tets="<<tot_touch<<" unhandled_tets="<<tot_unhandled<<endl;
	cout<<"degenerate triangles (rank-deficient field image): "<<g_degenerate_tris
			<<", of which whole-triangle-in-L: "<<g_tri_in_plane<<endl;
	cout<<"trials with >=1 bad cell: "<<trials_with_bad<<"/"<<trials
			<<" ("<<(100.0*trials_with_bad/trials)<<"%)"<<endl;
	cout<<"trials with >=1 multi-cycle cell: "<<trials_with_multi<<"/"<<trials<<endl;
	if(worst_bad_trial>=0)
		cout<<"worst single trial: #"<<worst_bad_trial<<" with "<<worst_bad<<" bad cells"<<endl;

	cout<<"bad cells by refinement level:";
	for(int l=0;l<32;++l) if(bad_by_level[l]) cout<<" ["<<l<<"]="<<bad_by_level[l];
	cout<<endl;
	cout<<"ok cells by refinement level:";
	for(int l=0;l<32;++l) if(ok_by_level[l]) cout<<" ["<<l<<"]="<<ok_by_level[l];
	cout<<endl;

	if(!dump_dir.empty())
		cout<<"failing-cell dumps (up to budget) written under: "<<dump_dir<<endl;

	return 0;
}
