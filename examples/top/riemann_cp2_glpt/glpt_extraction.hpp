#ifndef GLPT_EXTRACTION_HPP
#define GLPT_EXTRACTION_HPP

/*! \file
 * \brief Per-cell cycle extraction: turns face_crossing_cache's raw
 * root points into actual connected polygons -- the third (and, for
 * now, last) missing piece of a glpt_tree-based riemann_cp2.cpp
 * (glpt_points.hpp: real coordinates; glpt_crossing.hpp: root-finding
 * on 2-faces). Ports riemann_cp2.cpp's own crossing_node/
 * compute_crossing_node/same_crossing_node and its Phase-3 per-cell
 * graph-building + cycle-decomposition loop, adapted to glpt's own
 * facet/2-face combinatorics and vertex ids in place of nmt<4>'s
 * Simplex(T,3)/Simplex(T,2) handles and Vertex(T).desc.
 *
 * ============================================================
 * NODE IDENTITY: WHY THIS NEEDS TO BE EXPLICITLY CONSTRUCTED HERE
 * ============================================================
 * riemann_cp2.cpp's crossing_node identity (nd.desc) is just whichever
 * nmt<4> topological handle (a Vertex/Edge/Simplex(T,2)) the crossing
 * sits on -- nmt already guarantees two different cells reaching the
 * SAME vertex/edge/face get the SAME handle, for free, since it's one
 * shared topological structure. glpt has no such shared structure (by
 * design -- see glpt_tree.hpp's own header comment on why): a "vertex"
 * is just a small integer (glpt_vertex_ids.hpp), an "edge" has no
 * stored identity at all, and a "face" is only ever a vertex-id-triple
 * key into face_crossing_cache (glpt_crossing.hpp). So THIS file builds
 * the equivalent identity by hand: a dim-0 node's identity is its
 * global vertex id; a dim-1 node's is a canonically-ordered (sorted)
 * pair of global vertex ids plus a parameter along that edge; a dim-2
 * node's is the SAME packed face_key face_crossing_cache itself already
 * uses, plus which of its (up to 2) roots. All three are exactly the
 * quantities already available from glpt_vertex_ids.hpp/glpt_crossing.hpp
 * -- nothing new needs to be invented, just assembled per cell.
 */

#include "glpt_crossing.hpp"
#include <vector>
#include <cmath>

const double CROSSING_SNAP_TOL=1e-6;  // classification: on this vertex/edge or not?
const double CROSSING_MERGE_TOL=1e-9; // identity: same edge-parameter as another node?

struct crossing_node {
	int dim; // 0=vertex, 1=edge, 2=face
	int a,b; // dim==0: a=global vertex id (b unused); dim==1: sorted global vertex id pair
	double param; // dim==1 only: fraction of the way from a to b
	face_key fkey; int sub; // dim==2 only: exact packed face key + which root (0/1)
	pt3 p; // the actual point (for display, and the ntn==4 pairing heuristic)
};

inline bool same_crossing_node(const crossing_node& x, const crossing_node& y) {
	if(x.dim!=y.dim) return false;
	if(x.dim==0) return x.a==y.a;
	if(x.dim==1) return x.a==y.a && x.b==y.b && std::fabs(x.param-y.param)<CROSSING_MERGE_TOL;
	return x.fkey==y.fkey && x.sub==y.sub;
}

//! face_ids/bary are relative to the SAME face, in the SAME order, that
//! produced this root (face_crossing_cache::get()'s own p[]/id[]
//! order) -- see this file's own header comment for the identity
//! construction.
inline void compute_crossing_node(const int face_ids[3], const double bary[3], const face_key& fkey, int sub, const pt3& p, crossing_node& nd) {
	double l[3]={bary[0],bary[1],bary[2]};
	int zeros=0, zi[3];
	for(int b=0;b<3;++b) if(std::fabs(l[b])<CROSSING_SNAP_TOL) zi[zeros++]=b;
	nd.p=p;
	if(zeros>=2) {
		int keep = (zeros==3) ? 0 : (3-zi[0]-zi[1]); // the one index not in zi[]
		nd.dim=0; nd.a=face_ids[keep];
	} else if(zeros==1) {
		int lo=(zi[0]==0)?1:0, hi=(zi[0]==2)?1:2;
		int va=face_ids[lo], vb=face_ids[hi];
		double t = l[hi]/(l[lo]+l[hi]);
		if(va>vb) { int tmp=va; va=vb; vb=tmp; t=1.0-t; } // canonicalize (a<b), keep param meaning "fraction toward b"
		nd.dim=1; nd.a=va; nd.b=vb; nd.param=t;
	} else {
		nd.dim=2; nd.fkey=fkey; nd.sub=sub;
	}
}

// --- glpt 4-cell facet -> 2-face mapping: facet j (the tetrahedron
// opposite local vertex j) is bounded by exactly the 4 of
// GLPT_CELL_FACES's 10 entries that don't involve j either -- computed
// by hand once (DIM=4 fixed), matching GLPT_CELL_FACES's own order.
static const int GLPT_FACET_SUBFACES[5][4] = {
	{6,7,8,9}, // facet 0 (excl. vertex 0): {1,2,3},{1,2,4},{1,3,4},{2,3,4}
	{3,4,5,9}, // facet 1 (excl. vertex 1): {0,2,3},{0,2,4},{0,3,4},{2,3,4}
	{1,2,5,8}, // facet 2 (excl. vertex 2): {0,1,3},{0,1,4},{0,3,4},{1,3,4}
	{0,2,4,7}, // facet 3 (excl. vertex 3): {0,1,2},{0,1,4},{0,2,4},{1,2,4}
	{0,1,3,6}, // facet 4 (excl. vertex 4): {0,1,2},{0,1,3},{0,2,3},{1,2,3}
};

// --- Transversality score of a single crossing root: w_k := Fa*da_k +
// Fb*db_k is the ordinary complex directional derivative of F along a
// 2-face edge (da_k,db_k) (chart-local, from dehomogenizing the face's
// 3 vertices), since F is holomorphic. The score
//   |Im(conj(w1)*w2)| / (|e1| |e2| (|Fa|^2+|Fb|^2))
// is in [0,1] by Cauchy-Schwarz, and is exactly 0 when the 2-face's own
// (flat) tangent plane is tangent to the curve at this root -- every
// real direction in it would then have zero directional derivative.
// Equivalently: the face's two edge vectors are both close to
// perpendicular to grad(Re F) AND grad(Im F) -- which collapse into
// this single complex test because grad(Re F) is always perpendicular
// to grad(Im F), with equal norm (Cauchy-Riemann). Measured
// (2026-09-27, all 6 catalog curves, depth 14) to predict extraction
// bad-cell risk with a strong, mostly monotonic gradient: 22% bad rate
// for score<=0.05, down to ~3% for score in (0.6,0.8] -- see
// --tangency-threshold in riemann_cp2_glpt.cpp, which uses this to
// force extra refinement on suspiciously near-tangent cells BEFORE
// extraction, rather than only repairing bad ones after the fact.
inline double crossing_transversality(const pt3 face_pts[3], const double bary[3]) {
	int chart = pick_chart(face_pts);
	cx a[3], b[3];
	for(int k=0;k<3;++k) dehomogenize(chart, face_pts[k], a[k], b[k]);
	cx da1=a[1]-a[0], db1=b[1]-b[0];
	cx da2=a[2]-a[0], db2=b[2]-b[0];
	double e1n = std::sqrt(std::norm(da1)+std::norm(db1));
	double e2n = std::sqrt(std::norm(da2)+std::norm(db2));
	if(e1n<1e-14 || e2n<1e-14) return 1.0; // degenerate (near-zero-length) edge -- not a tangency signal
	cx ra = bary[0]*a[0]+bary[1]*a[1]+bary[2]*a[2];
	cx rb = bary[0]*b[0]+bary[1]*b[1]+bary[2]*b[2];
	cx Fa = (chart==0)?Fa_chart0(ra,rb):(chart==1)?Fa_chart1(ra,rb):Fa_chart2(ra,rb);
	cx Fb = (chart==0)?Fb_chart0(ra,rb):(chart==1)?Fb_chart1(ra,rb):Fb_chart2(ra,rb);
	double gnorm2 = std::norm(Fa)+std::norm(Fb);
	if(gnorm2<1e-24) return 0.0; // at/near a singular point of the curve -- treat as maximally degenerate
	cx w1 = Fa*da1 + Fb*db1;
	cx w2 = Fa*da2 + Fb*db2;
	return std::abs(std::imag(std::conj(w1)*w2)) / (e1n*e2n*gnorm2);
}

struct cell_extraction_result {
	std::vector<crossing_node> nodes;
	std::vector<std::vector<int> > cycles; // each a closed walk of indices into nodes
	int ntouching_tets, nbad_tets;
	bool any_edge;     // false: the curve doesn't cross this cell at all -- not a failure, just nothing to extract
	bool decompose_ok; // true only when any_edge and the graph cleanly decomposed into cycles
	double min_transversality; // min crossing_transversality() over every root examined; 1.0 (safe default) if none
};

//! Builds the crossing-node graph for one glpt 4-cell (its 5 facets'
//! worth of 2-face root queries, via fcache) and decomposes it into
//! closed cycles -- one per boundary component of the curve's
//! intersection with this cell. Faithful port of riemann_cp2.cpp's own
//! Phase-3 per-cell loop (see that file, "for(int j=0;j<=DIM;++j)"
//! onward) -- same ntn==1/2/4/other cases, same degree check, same
//! cycle-walk.
inline void extract_cell(const pt3 pts[glpt::DIM+1], const int ids[glpt::DIM+1], face_crossing_cache& fcache,
		cell_extraction_result& out) {
	const int DIM = glpt::DIM;
	out.nodes.clear();
	out.cycles.clear();
	out.ntouching_tets=0; out.nbad_tets=0;
	out.min_transversality=1.0;
	std::vector<std::vector<int> > adj;
	bool tet_unhandled=false;

	for(int j=0;j<=DIM;++j) {
		int tnodes[4]; int ntn=0;
		for(int s=0;s<4;++s) {
			int f = GLPT_FACET_SUBFACES[j][s];
			pt3 face_pts[3]; int face_ids[3];
			for(int k=0;k<3;++k) { face_pts[k]=pts[GLPT_CELL_FACES[f][k]]; face_ids[k]=ids[GLPT_CELL_FACES[f][k]]; }
			const face_result& fr = fcache.get(face_pts, face_ids);
			for(int r=0;r<fr.nroots;++r) {
				double tv = g_geodesic_faces ? crossing_transversality_geo(face_pts, face_ids, fcache.root_bary(fr,r))
				                             : crossing_transversality(face_pts, fcache.root_bary(fr,r));
				if(tv<out.min_transversality) out.min_transversality=tv;
				crossing_node nd;
				compute_crossing_node(face_ids, fcache.root_bary(fr,r), make_face_key(face_ids[0],face_ids[1],face_ids[2]), r,
					fcache.root_point(fr,r), nd);
				int idx=-1;
				for(size_t q=0;q<out.nodes.size();++q) if(same_crossing_node(out.nodes[q],nd)) { idx=(int)q; break; }
				if(idx<0) { idx=(int)out.nodes.size(); out.nodes.push_back(nd); adj.push_back(std::vector<int>()); }
				bool dup=false;
				for(int q=0;q<ntn;++q) if(tnodes[q]==idx) dup=true;
				if(!dup && ntn<4) tnodes[ntn++]=idx;
			}
		}
		if(ntn==2) {
			adj[tnodes[0]].push_back(tnodes[1]);
			adj[tnodes[1]].push_back(tnodes[0]);
		} else if(ntn==4) {
			static const int pairings[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
			int best=0; double bestcost=1e300;
			for(int c=0;c<3;++c) {
				double cost=0;
				for(int e=0;e<2;++e) {
					const crossing_node& A=out.nodes[tnodes[pairings[c][2*e]]];
					const crossing_node& B=out.nodes[tnodes[pairings[c][2*e+1]]];
					cost+=fs_dist(A.p,B.p);
				}
				if(cost<bestcost) { bestcost=cost; best=c; }
			}
			for(int e=0;e<2;++e) {
				int a=tnodes[pairings[best][2*e]], b=tnodes[pairings[best][2*e+1]];
				adj[a].push_back(b); adj[b].push_back(a);
			}
		} else if(ntn==1) {
			++out.ntouching_tets;
		} else if(ntn!=0) {
			++out.nbad_tets;
			tet_unhandled=true;
		}
	}

	bool deg_ok=true;
	out.any_edge=false;
	for(size_t q=0;q<out.nodes.size();++q) {
		size_t dg=adj[q].size();
		if(dg==0) continue;
		out.any_edge=true;
		if(dg!=2) deg_ok=false;
	}
	out.decompose_ok = out.any_edge && deg_ok && !tet_unhandled;
	if(!out.decompose_ok) return;

	std::vector<bool> seen(out.nodes.size(),false);
	for(size_t q=0;q<out.nodes.size();++q) {
		if(seen[q] || adj[q].empty()) continue;
		std::vector<int> cyc;
		int start=(int)q, prev=start, cur=adj[start][0];
		seen[q]=true; cyc.push_back(start);
		while(cur!=start) {
			if(seen[cur]) { out.decompose_ok=false; return; }
			seen[cur]=true; cyc.push_back(cur);
			int nxt=(adj[cur][0]==prev) ? adj[cur][1] : adj[cur][0];
			prev=cur; cur=nxt;
			if(cyc.size()>(size_t)(DIM+1)) { out.decompose_ok=false; return; }
		}
		if(cyc.size()<3) { out.decompose_ok=false; return; }
		out.cycles.push_back(cyc);
	}
}

#endif // GLPT_EXTRACTION_HPP
