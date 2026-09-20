#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <queue>
#include <set>
#include <map>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <random>
#include <limits>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/euler.hpp>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/top/do_nothing.hpp>
#include <vgtl/alg/vec.hpp>
#include "functions_pc2.hpp"

// Seed mesh: A. Gaifullin's 15-vertex/108-cell triangulation of CP^2
// (arXiv:0904.4222, also GAP simpcomp SCLib entry 397) -- Maubach-
// compatible and coherently orientable DIRECTLY, no subdivision step
// (see gaifullin_points()/gaifullin_cells() below for the full story,
// including why this replaced the project's earlier Kuhnel-9-vertex +
// barycentric-subdivision seed). Unlike examples/top/riemann/ (built on
// C_infty x C_infty = P^1 x P^1, w and z treated as two separate
// factors with all the asymmetric-pole handling that entailed),
// everything here lives in a single CP^2: one homogeneous point
// [X:Y:Z] per vertex, one curve F(X,Y,Z)=0.

#define DIM 4

using namespace std;
using namespace vgtl;

typedef vgtl::nmt<DIM> T;

T t;

typedef vec<3,cx> pt3;

pt3 mkpt(cx a, cx b, cx c) { pt3 p; p[0]=a; p[1]=b; p[2]=c; return p; }

// Storage-only half-precision twin of pt3/cx (see narrow_pt/widen_pt
// below) -- isolated test of extra_data<2>::pt[] alone (the Newton root
// cache, read back only for chart-selection heuristics and final mesh
// export), NOT extra_data<0>::p (vertex positions, which feed
// fs_barycenter's recursive averaging across refinement depth).
typedef complex<float> cxf;
typedef vec<3,cxf> pt3f;
pt3f narrow_pt(const pt3& p) {
	pt3f q;
	for(int i=0;i<3;++i) q[i]=cxf((float)p[i].real(),(float)p[i].imag());
	return q;
}
pt3 widen_pt(const pt3f& p) {
	pt3 q;
	for(int i=0;i<3;++i) q[i]=cx((double)p[i].real(),(double)p[i].imag());
	return q;
}

// Pool for the Newton root cache: the overwhelming majority of
// triangles have ZERO roots (measured: 7083822/7211426 = 98.2% for the
// parabola baseline, still 98.2% for the quartic), so reserving 2
// slots' worth of storage (pt3f pt[2] + double l1[2],l2[2]) INSIDE
// every extra_data<2> paid for a near-empty array on ~98% of triangles.
// Instead, extra_data<2> holds a single index into this flat pool
// (root_idx, -1 if empty); a triangle with nroots>0 owns the
// contiguous run g_roots[root_idx .. root_idx+nroots-1]. Never shrinks
// (matches nmt's own no-real-deletion semantics -- see del() for k!=Dim
// in model/nmt.hpp), same lifetime pattern extra_data itself already
// has.
struct root_data { pt3f pt; double l1,l2; };
vector<root_data> g_roots;

namespace vgtl {

	template <>
	struct extra_data<0> {
		// Gradient/F(p) used to live here too (cx Fx,Fy,Fz,Fval), cached
		// once at vertex creation. But cell_priority()/proximity_factor()
		// (the only readers) are each called exactly once per CELL, at
		// push time into the refinement priority queue -- never
		// recomputed on pop -- so caching per VERTEX instead only pays
		// off if a vertex's average cell-incidence is below the
		// gradient's own recompute cost, which it isn't here (F is a low-
		// degree polynomial, eval_poly3 is cheap). Recomputing on demand
		// (see vertex_gradient()/vertex_Fval() below) trades a bit of
		// redundant polynomial evaluation for dropping this struct from
		// 112 to 48 bytes/vertex.
		pt3 p;           // homogeneous point, kept unit-normalized in C^3
	};

	pt3 point(const T& t, Vertex(T) v) { return attr(t,v)->p; }
	void point_set(T& t, Vertex(T) v, const pt3& p) { attr(t,v)->p=p; }

	// Cache of the F=0 intersection point(s) of a 2-simplex, exactly as in
	// examples/top/riemann/riemann.cpp -- up to 2 roots per triangle,
	// identified by barycentric coordinates (l1,l2), l0=1-l1-l2. Only ONE
	// homogeneous point per root now (not a separate w/z pair).
	template <>
	struct extra_data<2> {
		// Index into the global g_roots pool (see above), not the roots
		// themselves -- most triangles have none (measured: 98.2% for
		// both the parabola and quartic catalog curves), so paying for
		// 2 root slots inline here was wasted on nearly every triangle.
		// -1 = no roots. A triangle with nroots>0 owns the contiguous
		// run g_roots[root_idx .. root_idx+nroots-1].
		int root_idx;
		// --proximity/--bernstein (see compute_bernstein_bounds): a
		// certified real-interval enclosure of Re(F) and Im(F) over this
		// 2-simplex, from its Bernstein-Bezier coefficients -- immutable
		// once the triangle's 3 vertices exist, cached the same way the
		// Newton roots above are.
		//
		// float, not double: this is a certified enclosure used only for
		// pruning/priority, never fed back into the geometry, so halving
		// its footprint (the field that dominates the per-triangle memory
		// budget at deep --bernstein refinement) is a clean win as long as
		// the stored interval still contains the double-precision one --
		// see the outward rounding in compute_bernstein_bounds() below.
		float reLo,reHi,imLo,imHi;
		// unsigned char, not bool/int: nroots is always 0..2 (see comment
		// above), and grouping the three 1-byte flags here (after every
		// 4/8-byte-aligned member) means they cost only their own 3 bytes
		// instead of also forcing padding around themselves.
		unsigned char computed;
		unsigned char nroots;
		unsigned char bernstein_computed;
		extra_data() : root_idx(-1), computed(0), nroots(0), bernstein_computed(0) {}
	};

}

// --- Fubini-Study geometry on CP^2 (unit representatives in C^3) -------
// Hermitian inner product <a,b> = sum a_i * conj(b_i).
cx hdot(const pt3& a, const pt3& b) {
	cx r(0,0);
	for(int i=0;i<3;++i) r += a[i]*std::conj(b[i]);
	return r;
}
double hnorm(const pt3& a) { return std::sqrt(std::max(0.0, hdot(a,a).real())); }
pt3 normalize3(pt3 a) {
	double n=hnorm(a);
	for(int i=0;i<3;++i) a[i]/=n;
	return a;
}
// Fubini-Study distance between two unit representatives: arccos|<a,b>|,
// the direct analogue of the old sphere_dist (chordal S^2 distance).
double fs_dist(const pt3& a, const pt3& b) {
	double ip=std::abs(hdot(a,b));
	if(ip>1.0) ip=1.0;
	return std::acos(ip);
}
// Rephase b so that <a,b> is real and non-negative -- resolves the U(1)
// ambiguity in b's homogeneous representative relative to a, needed
// before b can be meaningfully averaged with a (see fs_barycenter).
pt3 align_phase(const pt3& a, pt3 b) {
	cx ip=hdot(a,b);
	double m=std::abs(ip);
	if(m<1e-14) return b; // near-orthogonal (FS-distance ~pi/2): no sensible phase to lock onto
	cx phase=ip/m;
	for(int i=0;i<3;++i) b[i]*=phase;
	return b;
}
// Geodesic barycenter of k+1 points on CP^2 (Fubini-Study): rephase each
// against the running sum and renormalize -- the direct C^3 analogue of
// slerp_midpoint's "sum unit vectors, renormalize" on the real sphere
// S^2 (which is itself the Fubini-Study geometry of CP^1). For k=1 this
// is exactly the geodesic midpoint used by Maubach edge-bisection; for
// k+1>2 it is what barycentric_vertex() below uses to place a new vertex
// at the barycenter of a face during the seed's barycentric subdivision.
pt3 fs_barycenter(const vector<pt3>& pts) {
	pt3 acc=pts[0];
	for(size_t i=1;i<pts.size();++i) {
		pt3 a=align_phase(acc,pts[i]);
		for(int c=0;c<3;++c) acc[c]+=a[c];
	}
	return normalize3(acc);
}

// --- Random unitary change of basis on C^3 (--generic), applied to the
// seed mesh's own vertex positions (NOT to F, which is left exactly as
// catalogued) to break any alignment between the curve's fixed X,Y,Z
// basis and the seed mesh's own: 3 of the 15 seed points (the V4\{e}
// ones, see gaifullin_points() below) have TWO homogeneous coordinates
// equal to 0, a real, non-generic resonance with that basis, and
// everything basis-dependent downstream (Fx,Fy,Fz used by
// cell_priority, per-triangle chart selection in
// triangle_intersection) is evaluated in that same fixed basis. This is
// the direct CP^2 analogue of examples/top/riemann's --generic seed
// rotation (see riemann_pole_rotation_and_viz_gaps memory), done here
// by rotating the MESH rather than the polynomial: cheaper (no
// symbolic linear-substitution/expansion of F needed), and every
// already-verified topological property (Euler char, orientation BFS)
// is purely combinatorial, so it's untouched by this.
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
bool g_generic=false;
unsigned g_generic_seed=12345;

// --- The curve and its derivatives, generic (see functions_pc2.hpp) ---
// Only the gradient is needed now: refinement priority is gradient
// dispersion across a cell's own vertices (see cell_priority below), not
// a Hessian-based per-point signal, so no second derivatives are kept.
poly_F3 g_F, g_Fx, g_Fy, g_Fz;
int g_d=0;

void set_curve(const poly_F3& F) {
	g_F=F; g_d=F.d;
	g_Fx=diff_poly3(F,0); g_Fy=diff_poly3(F,1); g_Fz=diff_poly3(F,2);
}

// On-demand replacements for the per-vertex Fx/Fy/Fz/Fval cache: same
// eval_poly3 calls, just done at cell_priority/proximity_factor time
// (each vertex is looked at only a handful of times per call, and each
// call happens only once per cell -- see the extra_data<0> comment).
pt3 vertex_gradient(const T& t, Vertex(T) v) {
	pt3 p=point(t,v);
	pt3 g;
	g[0]=eval_poly3(g_Fx,p[0],p[1],p[2]);
	g[1]=eval_poly3(g_Fy,p[0],p[1],p[2]);
	g[2]=eval_poly3(g_Fz,p[0],p[1],p[2]);
	return g;
}
cx vertex_Fval(const T& t, Vertex(T) v) {
	pt3 p=point(t,v);
	return eval_poly3(g_F,p[0],p[1],p[2]);
}

// --- Seed: Gaifullin's 15-vertex triangulation X of CP^2 (A. Gaifullin,
// arXiv:0904.4222, Construction 1.1 + Sec. 3's explicit coordinates;
// independently cross-checked against GAP's simpcomp SCLib entry 397,
// "Gaifullin CP^2" -- identical f-vector (15,90,240,270,108) and facet
// list up to relabeling). REPLACES Kuhnel's 9-vertex/36-cell CP^2_9 +
// barycentric subdivision used earlier this project: Kuhnel's raw
// 36-cell complex has a PROVEN bipartiteness obstruction to coherent
// orientation under any Maubach-compatible reordering (see
// riemann_pc2_new_approach memory / find_reordering.cpp), forcing the
// barycentric-subdivision workaround -- which in turn introduced its
// own bug (3 of the 9 Hesse points are exact vector sums of two others,
// P[i]+P[6]=P[i+3], so 3 subdivision-created flag-centroids collapse
// exactly onto existing vertices, producing 48 zero-length-edge cells
// that Maubach's longest-edge rule can never split). Gaifullin's
// triangulation's dual graph IS bipartite (verified: BFS 2-coloring,
// 0 conflicts) and admits a proper 5-edge-coloring (verified: exact
// backtracking search) -- so it gets BOTH Maubach compatibility AND
// coherent orientation directly, with no subdivision step and no
// vertex-coincidence risk (there's no vertex-averaging construction
// left to coincide). Also directly verified: none of its 240 actual
// 2-faces are projectively collinear (all 240 determinants nonzero).
//
// Vertex set V = (V4\{e}) x-union ({1,2,3,4}x{1,2,3}), V4 the Klein
// four-group in S4; omega = primitive cube root of unity. Labels 0-2
// are V4\{e} = (12)(34),(13)(24),(14)(23); labels 3-14 are (a,b) for
// a=1..4 (4 blocks of 3), b=1..3, in that order -- exactly the order
// used below and in gaifullin_cells[].
vector<pt3> gaifullin_points() {
	double s3=std::sqrt(3.0);
	cx one(1,0), zero(0,0);
	cx omega(-0.5, s3/2.0), omega2(-0.5,-s3/2.0);
	vector<pt3> P(15);
	P[0]=mkpt(zero,zero,one);           // (12)(34)
	P[1]=mkpt(zero,one,zero);           // (13)(24)
	P[2]=mkpt(one,zero,zero);           // (14)(23)
	P[3]=mkpt(-one,omega,omega2);       // (1,1)
	P[4]=mkpt(-one,omega2,omega);       // (1,2)
	P[5]=mkpt(-one,one,one);            // (1,3)
	P[6]=mkpt(one,-omega,omega2);       // (2,1)
	P[7]=mkpt(one,-omega2,omega);       // (2,2)
	P[8]=mkpt(one,-one,one);            // (2,3)
	P[9]=mkpt(one,omega,-omega2);       // (3,1)
	P[10]=mkpt(one,omega2,-omega);      // (3,2)
	P[11]=mkpt(one,one,-one);           // (3,3)
	P[12]=mkpt(one,omega,omega2);       // (4,1)
	P[13]=mkpt(one,omega2,omega);       // (4,2)
	P[14]=mkpt(one,one,one);            // (4,3)
	for(int i=0;i<15;++i) P[i]=normalize3(P[i]);
	if(g_generic) {
		cx M[3][3]; random_unitary(M,g_generic_seed);
		for(int i=0;i<15;++i) P[i]=apply_unitary(M,P[i]);
	}
	return P;
}

// --- The 108 four-simplices of Gaifullin's triangulation, ALREADY in
// Maubach-compatible per-cell vertex order (local position = the color
// a backtracking 5-edge-coloring search assigned to the facet opposite
// that vertex) -- derived and verified offline (facet-owner histogram
// {2:270}, 0 same-index mismatches, 0 orientation-coherence violations
// under the corresponding bipartite 2-coloring). Orientation itself is
// NOT baked in here: main() just seeds it arbitrarily and calls the
// same reorient_via_bfs() Kuhnel's seed already used, which re-derives
// it from this vertex order and (per the above) is guaranteed to find
// it coherent.
typedef vgtl::array<int,5> cell5;
static const int gaifullin_cells_raw[108][5] = {
	{0,3,7,9,13}, {0,3,7,9,14}, {0,3,7,10,12}, {0,3,7,10,14}, {0,3,7,11,12}, {0,3,7,11,13},
	{0,3,8,9,13}, {0,3,8,9,14}, {0,3,8,10,12}, {0,3,8,10,14}, {0,3,8,11,12}, {0,3,8,11,13},
	{0,4,6,9,13}, {0,4,6,9,14}, {0,4,6,10,12}, {0,4,6,10,14}, {0,4,6,11,12}, {0,4,6,11,13},
	{0,4,8,9,13}, {0,4,8,9,14}, {0,4,8,10,12}, {0,4,8,10,14}, {0,4,8,11,12}, {0,4,8,11,13},
	{0,5,6,9,13}, {0,5,6,9,14}, {0,5,6,10,12}, {0,5,6,10,14}, {0,5,6,11,12}, {0,5,6,11,13},
	{0,5,7,9,13}, {0,5,7,9,14}, {0,5,7,10,12}, {0,5,7,10,14}, {0,5,7,11,12}, {0,5,7,11,13},
	{1,3,6,10,13}, {1,3,6,10,14}, {1,3,6,11,13}, {1,3,6,11,14}, {1,3,7,10,12}, {1,3,7,10,14},
	{1,3,7,11,12}, {1,3,7,11,14}, {1,3,8,10,12}, {1,3,8,10,13}, {1,3,8,11,12}, {1,3,8,11,13},
	{1,4,6,9,13}, {1,4,6,9,14}, {1,4,6,11,13}, {1,4,6,11,14}, {1,4,7,9,12}, {1,4,7,9,14},
	{1,4,7,11,12}, {1,4,7,11,14}, {1,4,8,9,12}, {1,4,8,9,13}, {1,4,8,11,12}, {1,4,8,11,13},
	{1,5,6,9,13}, {1,5,6,9,14}, {1,5,6,10,13}, {1,5,6,10,14}, {1,5,7,9,12}, {1,5,7,9,14},
	{1,5,7,10,12}, {1,5,7,10,14}, {1,5,8,9,12}, {1,5,8,9,13}, {1,5,8,10,12}, {1,5,8,10,13},
	{2,3,6,10,13}, {2,3,6,10,14}, {2,3,6,11,13}, {2,3,6,11,14}, {2,3,7,9,13}, {2,3,7,9,14},
	{2,3,7,11,13}, {2,3,7,11,14}, {2,3,8,9,13}, {2,3,8,9,14}, {2,3,8,10,13}, {2,3,8,10,14},
	{2,4,6,10,12}, {2,4,6,10,14}, {2,4,6,11,12}, {2,4,6,11,14}, {2,4,7,9,12}, {2,4,7,9,14},
	{2,4,7,11,12}, {2,4,7,11,14}, {2,4,8,9,12}, {2,4,8,9,14}, {2,4,8,10,12}, {2,4,8,10,14},
	{2,5,6,10,12}, {2,5,6,10,13}, {2,5,6,11,12}, {2,5,6,11,13}, {2,5,7,9,12}, {2,5,7,9,13},
	{2,5,7,11,12}, {2,5,7,11,13}, {2,5,8,9,12}, {2,5,8,9,13}, {2,5,8,10,12}, {2,5,8,10,13},
};

vector<cell5> gaifullin_cells() {
	vector<cell5> cells(108);
	for(int i=0;i<108;++i) for(int k=0;k<5;++k) cells[i][k]=gaifullin_cells_raw[i][k];
	return cells;
}

// --- Orientation coherence check + BFS re-derivation, byte-for-byte
// the same algorithm as examples/top/riemann/riemann.cpp's
// count_incoherent()/orientation-BFS pass -- entirely topological, no
// dependency on CP^2 vs the old C_infty^2 specifics.
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

void reorient_via_bfs(T& t, int ncells) {
	// Byte-for-byte the same BFS as examples/top/riemann/riemann.cpp's,
	// run here right after seeding the 108 Gaifullin cells (all current
	// at this point, unlike the earlier Kuhnel-seed pipeline which ran
	// this AFTER barycentric subdivision, with stale not-current cells
	// also sitting in the container -- the is_current() filtering below
	// is kept anyway, harmless and cheap, in case this is ever called
	// somewhere less pristine).
	set<Cell(T)> seen;
	queue<Cell(T)> q;
	Cell_it(T) ci,cend;
	simplices(t,ci,cend);
	while(ci!=cend && !is_current(t,*ci)) ++ci;
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
			if(!is_current(t,c2)) continue;
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
			<<nconflict<<" conflicts (expected 0 -- Gaifullin's dual graph is "
			<<"bipartite, verified offline before this reordering was written)"<<endl;
}

// --- Maubach ApplyNew: geodesic (Fubini-Study) midpoint for new edge
// vertices, exactly as slerp_midpoint did on the real sphere before.
struct refine_app : do_nothing {
	using do_nothing::apply;
	vector<Cell(T)> new_cells;

	Vertex(T) add_edge_vertex(T& t, Edge(T) e) {
		vgtl::array<Vertex(T),2> vs;
		vertices(t,e,vs);
		Vertex(T) v=add(t);
		vector<pt3> pts(2);
		pts[0]=point(t,vs[0]); pts[1]=point(t,vs[1]);
		point_set(t,v,fs_barycenter(pts));
		return v;
	}

	void apply(T& t, Cell(T) s) {
		new_cells.push_back(s);
	}
};

// --- Refinement priority: gradient dispersion across a cell's own 5
// vertices. Replaces the earlier Hessian/flex-proximity criterion --
// that measured ramification of the curve's Gauss map (a genuinely
// different, projectively meaningful thing), not curvature, and a flex
// is a perfectly smooth point with no elevated risk of a triangle probe
// missing or double-counting a crossing. What actually creates that
// risk is the curve's tangent plane bending sharply (or, worse, nearing
// tangency with a probe triangle's own plane) across a cell -- i.e.
// curvature, which is exactly what varying gradient direction measures.
//
// This is Hermitian "PCA" on the cell's own (unit-normalized) gradients
// ghat_0..ghat_4, without ever building or diagonalizing the 3x3
// covariance-like matrix M = sum_i ghat_i ghat_i^dagger: trace(M) = n
// trivially (each ghat_i is unit), and trace(M^2), which is all that's
// needed, is exactly sum_{i,j} |<ghat_i,ghat_j>|^2 -- a plain double sum
// of Hermitian inner products already available via hdot(). kappa =
// trace(M^2)/trace(M)^2 lies in [1/3,1]: 1 means every vertex's gradient
// points the same way (flat, up to an overall phase), 1/3 is the most
// spread 5 unit vectors in C^3 can be. priority = 1-kappa is then 0 when
// flat and up to 2/3 when the tangent plane twists sharply -- and it is
// automatically scaled by the cell's own size (a bigger cell samples
// gradient direction over more of the curve's bend, hence more
// dispersion, even at constant intrinsic curvature), so no separate
// diameter term is needed the way the old branch_gap-based ratio did.
//
// Optional (--proximity) extra factor: gradient dispersion alone is a
// GLOBAL signal, evaluated identically whether or not a cell is
// anywhere near F=0 -- it can and does chase curvature-irrelevant
// gradient-direction variation elsewhere in CP^2 (see
// riemann_pc2_gradient_dispersion memory: the parabola's capped-
// refinement fraction not shrinking with depth). To bias refinement
// toward the curve itself without breaking convergence, weight by
// proximity = diam/(mind+diam), where mind = min over the cell's own
// vertices of a first-order distance-to-curve estimate |F(v)|/|gradF(v)|
// (same idea as newton_on_triangle's out_geom_resid: dividing by the
// gradient magnitude cancels F's arbitrary scale, giving units
// comparable to the cell's own Fubini-Study diameter) and diam is the
// cell's own diameter. This is NOT plain 1/|F|: multiplying by 1/|F|
// directly would blow up priority exactly as a cell shrinks toward the
// curve (mind->0 together with the cell itself), which is the opposite
// of convergence. Comparing mind against diam instead saturates at 1
// (even when mind=0 exactly) once the curve is within about one cell
// diameter, and decays smoothly toward 0 once the curve is many cell
// diameters away.
bool g_proximity=false;

double cell_priority(const T& t, Cell(T) cv) {
	vgtl::array<Vertex(T),DIM+1> vs;
	vertices(t,cv,vs);
	pt3 ghat[DIM+1], pts[DIM+1];
	int n=0;
	double mind=1e300;
	for(int i=0;i<=DIM;++i) {
		pts[i]=point(t,vs[i]);
		pt3 g=vertex_gradient(t,vs[i]);
		double gn=hnorm(g);
		if(gn<1e-12) continue; // near a singular point (F, gradient both ~0): can't normalize
		if(g_proximity) {
			double dv=std::abs(vertex_Fval(t,vs[i]))/gn;
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
	double kappa=S/(double(n)*double(n));
	double priority=1.0-kappa;
	if(g_proximity) {
		double diam=0;
		for(int i=0;i<=DIM;++i)
			for(int j=i+1;j<=DIM;++j)
				diam=std::max(diam,fs_dist(pts[i],pts[j]));
		priority *= diam/(mind+diam);
	}
	return priority;
}

// --- Bernstein-certified refinement/pruning (--bernstein), a second,
// independent refinement criterion alongside gradient dispersion above
// -- an experiment in exploiting that F is a genuine polynomial (not
// just "some function"), per the user's proposal. Two phases:
//
// (1) For a 2-simplex, restricting F to it (X,Y,Z = l0 P0+l1 P1+l2 P2,
// real barycentric l0,l1,l2) gives a polynomial of the SAME degree d in
// (l0,l1,l2), homogeneous since F is -- i.e. exactly a degree-d
// triangular Bezier/Bernstein patch. Its Bernstein-Bezier coefficients
// enclose its range on the triangle (they're a partition-of-unity
// convex combination), so [min,max] of Re and Im of those coefficients
// is a CERTIFIED enclosure of Re(F) and Im(F) on that triangle -- not a
// heuristic. A cell can safely be dropped from refinement forever, not
// just deprioritized, once NONE of its 10 triangular faces has both
// ranges straddling 0 (dim(triangle)+dim(zero locus)-dim(ambient)=
// 2+2-4=0, same reasoning triangle_intersection already relies on: if
// the zero locus can't cross ANY of a cell's 2-faces, it can't be
// passing through the cell's interior either, by the same argument that
// makes the reconstructed polygon's vertices always triangle-crossings).
//
// (2) Extraction (triangle_intersection) reuses these same cached
// bounds to skip Newton's seed search entirely on a triangle already
// certified to have no root -- kept as Newton otherwise (per the user's
// choice), not replaced by a Bernstein-subdivision root isolator.
bool g_bernstein=false;
int g_bernstein_level=0; // --bernstein-level: see bernstein_bounds_recursive

// Forward declarations -- real definitions are with the rest of the
// chart machinery, further down (F_chart0/1/2 etc.), but pick_chart
// and compute_bernstein_bounds need to match triangle_intersection's
// own chart choice exactly (see the domain-mismatch note below), so
// they're declared here and defined the same as always down there.
void dehomogenize(int chart, const pt3& p, cx& a, cx& b);
pt3 rehomogenize(int chart, cx a, cx b);

// Picks the same best-conditioned chart triangle_intersection does
// (max over charts of the min |coordinate| across the 3 vertices) --
// factored out so both call sites are provably using the identical
// rule, not two copies that could drift apart.
int pick_chart(const pt3 p[3]) {
	int chart=0; double bestscore=-1;
	for(int c=0;c<3;++c) {
		double m=std::min(std::abs(p[0][c]), std::min(std::abs(p[1][c]),std::abs(p[2][c])));
		if(m>bestscore) { bestscore=m; chart=c; }
	}
	return chart;
}

// A degree-d homogeneous polynomial in (l0,l1,l2), stored as its
// monomial coefficients c[a][b] for l0^a l1^b l2^(d-a-b) (a+b<=d<=4 for
// every catalog curve, so a plain 5x5 grid suffices -- no map/tree
// needed for exponents this small).
struct bary_poly {
	int d;
	cx c[5][5];
	bary_poly() : d(0) { for(int a=0;a<5;++a) for(int b=0;b<5;++b) c[a][b]=cx(0,0); }
};

double fact3(int n) { double r=1; for(int i=2;i<=n;++i) r*=i; return r; }
double multinom3(int d,int a,int b,int cc) { return fact3(d)/(fact3(a)*fact3(b)*fact3(cc)); }

// (l0*X0+l1*X1+l2*X2)^n via the trinomial theorem.
bary_poly pow_linear3(cx X0,cx X1,cx X2,int n) {
	bary_poly r; r.d=n;
	for(int p=0;p<=n;++p)
		for(int q=0;q<=n-p;++q) {
			int rr=n-p-q;
			cx coeff=multinom3(n,p,q,rr)*ipow(X0,p)*ipow(X1,q)*ipow(X2,rr);
			r.c[p][q]+=coeff;
		}
	return r;
}

bary_poly mul_bary(const bary_poly& A, const bary_poly& B) {
	bary_poly r; r.d=A.d+B.d;
	for(int a1=0;a1<=A.d;++a1)
		for(int b1=0;a1+b1<=A.d;++b1) {
			if(A.c[a1][b1]==cx(0,0)) continue;
			for(int a2=0;a2<=B.d;++a2)
				for(int b2=0;a2+b2<=B.d;++b2)
					r.c[a1+a2][b1+b2]+=A.c[a1][b1]*B.c[a2][b2];
		}
	return r;
}

// F restricted to the triangle's own best-conditioned CHART (the exact
// same one triangle_intersection picks via pick_chart), as a uniform
// degree-F.d homogeneous polynomial in barycentric (l0,l1,l2) of the
// chart's 2 free coordinates (a,b) -- NOT yet Bernstein coefficients
// (see caller). This replaces an earlier version that composed F in
// the RAW ambient (X,Y,Z) affine barycentric combination instead: that
// bounded a genuinely different 2D patch than the one
// triangle_intersection actually searches (dehomogenize each vertex to
// the chart, THEN interpolate) -- the two coincide at the 3 corners but
// can diverge for larger triangles, and a --bernstein-selftest-style
// check plus real --bernstein-level>0 runs (see riemann_pc2_bernstein_
// refinement memory) showed the mismatch is real, not theoretical:
// tightening the OLD (wrong-domain) enclosure made extraction WORSE,
// silently dropping genuine crossings. Dehomogenizing first and
// building the SAME (a,b)-chart barycentric patch triangle_intersection
// implicitly searches removes that mismatch by construction.
//
// A term c*X^i*Y^j*Z^k, dehomogenized (one of X,Y,Z set to 1 per
// chart), contributes c*a^ea*b^eb of degree ea+eb<=F.d -- generally
// LESS than F.d now (whichever original exponent got dropped), so each
// term needs "degree elevation" (multiply by (l0+l1+l2)^(F.d-ea-eb),
// exactly 1 on the simplex l0+l1+l2=1, so this doesn't change the
// function's values there) before it can be added into a single
// uniform-degree grid -- standard CAGD technique, needed here because
// the dehomogenized F is no longer homogeneous the way the raw ambient
// one was.
bary_poly compose_to_chart(const poly_F3& F, int chart, cx a0,cx a1,cx a2, cx b0,cx b1,cx b2) {
	bary_poly acc; acc.d=F.d;
	for(size_t k=0;k<F.t.size();++k) {
		const term3& tm=F.t[k];
		int ea,eb;
		switch(chart) {
			case 0: ea=tm.e[1]; eb=tm.e[2]; break; // X=1: a=Y,b=Z
			case 1: ea=tm.e[0]; eb=tm.e[2]; break; // Y=1: a=X,b=Z
			default: ea=tm.e[0]; eb=tm.e[1]; break; // Z=1: a=X,b=Y
		}
		bary_poly pa=pow_linear3(a0,a1,a2,ea);
		bary_poly pb=pow_linear3(b0,b1,b2,eb);
		bary_poly term=mul_bary(pa,pb);
		int termdeg=ea+eb;
		if(termdeg<acc.d) {
			bary_poly elev=pow_linear3(cx(1,0),cx(1,0),cx(1,0),acc.d-termdeg); // (l0+l1+l2)^(F.d-termdeg)
			term=mul_bary(term,elev);
		}
		for(int a=0;a<=acc.d;++a)
			for(int b=0;a+b<=acc.d;++b)
				acc.c[a][b]+=tm.c*term.c[a][b];
	}
	return acc;
}

// A point in the triangle's chart (the 2 free affine coordinates, e.g.
// a=Y/X,b=Z/X for chart 0) -- the recursion below works entirely in
// this 2D space, matching triangle_intersection's own domain, rather
// than in the raw 3D ambient one.
struct cpt { cx a,b; };
cpt lerp_cpt(const cpt& p, const cpt& q) { cpt r; r.a=0.5*(p.a+q.a); r.b=0.5*(p.b+q.b); return r; }

void bernstein_leaf_bounds(const poly_F3& F, int chart, const cpt& P0, const cpt& P1, const cpt& P2,
														double& reLo, double& reHi, double& imLo, double& imHi) {
	bary_poly m=compose_to_chart(F,chart,P0.a,P1.a,P2.a,P0.b,P1.b,P2.b);
	reLo=1e300; reHi=-1e300; imLo=1e300; imHi=-1e300;
	for(int a=0;a<=m.d;++a) {
		for(int b=0;a+b<=m.d;++b) {
			int cc=m.d-a-b;
			cx beta=m.c[a][b]/multinom3(m.d,a,b,cc); // Bernstein coefficient
			double re=beta.real(), im=beta.imag();
			if(re<reLo) reLo=re; if(re>reHi) reHi=re;
			if(im<imLo) imLo=im; if(im>imHi) imHi=im;
		}
	}
}

void get_subtriangle(int s, const cpt& P0,const cpt& P1,const cpt& P2,
											const cpt& M01,const cpt& M12,const cpt& M20,
											cpt& A, cpt& B, cpt& C) {
	switch(s) {
		case 0: A=P0;  B=M01; C=M20; break;
		case 1: A=M01; B=P1;  C=M12; break;
		case 2: A=M20; B=M12; C=P2;  break;
		default: A=M01; B=M12; C=M20; break; // the middle, "upside-down" piece
	}
}

// Bernstein bounds for F restricted to the chart-triangle (P0,P1,P2),
// optionally tightened by `level` rounds of standard 1-to-4 triangular
// subdivision (connect edge midpoints, recurse on each of the 4
// children, union their bounds) before taking the enclosure -- level 0
// is exactly the single top-level patch. This is de Casteljau-style
// subdivision: it tightens the CERTIFIED enclosure (the union of 4
// smaller boxes is always a subset of, never bigger than, the parent
// box, and generically much smaller once each child sees less of the
// curve's own bend) WITHOUT touching the real simplicial mesh at all --
// P0,P1,P2 are always the triangle's own fixed (dehomogenized) chart
// coordinates; the recursion's extra points are purely virtual samples
// for this computation.
void bernstein_bounds_recursive(const poly_F3& F, int chart, const cpt& P0, const cpt& P1, const cpt& P2,
																 int level, double& reLo, double& reHi, double& imLo, double& imHi) {
	if(level<=0) { bernstein_leaf_bounds(F,chart,P0,P1,P2,reLo,reHi,imLo,imHi); return; }
	cpt M01=lerp_cpt(P0,P1), M12=lerp_cpt(P1,P2), M20=lerp_cpt(P2,P0);
	reLo=1e300; reHi=-1e300; imLo=1e300; imHi=-1e300;
	for(int s=0;s<4;++s) {
		cpt A,B,C; get_subtriangle(s,P0,P1,P2,M01,M12,M20,A,B,C);
		double srl,srh,sil,sih;
		bernstein_bounds_recursive(F,chart,A,B,C,level-1,srl,srh,sil,sih);
		if(srl<reLo) reLo=srl; if(srh>reHi) reHi=srh;
		if(sil<imLo) imLo=sil; if(sih>imHi) imHi=sih;
	}
}

// --bernstein-selftest: direct numerical verification, independent of
// the mesh, that (a) level-0 and level-1 enclosures both actually
// contain the TRUE range of Re(F)/Im(F) over the triangle's own CHART
// domain (the same one triangle_intersection searches -- see
// compose_to_chart's comment), sampled by brute force, and (b) level-1's
// box is a subset of level-0's (mathematically required: subdivision
// only tightens). Kept as a real regression check, not just a one-off:
// this is what caught the raw-ambient-vs-chart domain mismatch in the
// first place.
bool g_bernstein_selftest=false;

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
	std::mt19937 rng(1);
	std::uniform_real_distribution<double> ud(0.0,1.0);
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

void compute_bernstein_bounds(T& t, Simplex(T,2) tri) {
	extra_data<2>* d=attr(t,tri);
	if(d->bernstein_computed) return;
	d->bernstein_computed=true;

	vgtl::array<Vertex(T),3> vs;
	vertices(t,tri,vs);
	pt3 p[3]; for(int k=0;k<3;++k) p[k]=point(t,vs[k]);

	int chart=pick_chart(p);
	cpt Q0,Q1,Q2;
	dehomogenize(chart,p[0],Q0.a,Q0.b);
	dehomogenize(chart,p[1],Q1.a,Q1.b);
	dehomogenize(chart,p[2],Q2.a,Q2.b);

	// Recursion stays double for numerical stability; only the final,
	// immutable, cached result is narrowed to float. Round outward (one
	// float ULP via nextafterf) so the stored box is still a superset of
	// the double-precision enclosure -- the certified-enclosure property
	// (used by triangle_maybe_zero's pruning test) is preserved exactly,
	// not just approximately.
	double reLo,reHi,imLo,imHi;
	bernstein_bounds_recursive(g_F,chart,Q0,Q1,Q2,g_bernstein_level,
															 reLo,reHi,imLo,imHi);
	d->reLo=nextafterf((float)reLo,-numeric_limits<float>::infinity());
	d->reHi=nextafterf((float)reHi, numeric_limits<float>::infinity());
	d->imLo=nextafterf((float)imLo,-numeric_limits<float>::infinity());
	d->imHi=nextafterf((float)imHi, numeric_limits<float>::infinity());
}

bool triangle_maybe_zero(T& t, Simplex(T,2) tri) {
	compute_bernstein_bounds(t,tri);
	const extra_data<2>* d=attr(t,tri);
	return d->reLo<=0 && d->reHi>=0 && d->imLo<=0 && d->imHi>=0;
}

// Priority for --bernstein mode: the widest candidate-face box (the
// face least converged so far), or -1 (a sentinel meaning "certified
// empty -- drop this cell, never refine it again") if none of the
// cell's 10 triangular faces can contain a zero.
long g_bernstein_pruned_faces=0; // diagnostic: triangles Newton never had to touch

// Shared with cell_priority above (--proximity): diam/(mind+diam),
// mind=min over the cell's own vertices of |F(v)|/|gradF(v)|, diam=the
// cell's own Fubini-Study diameter. For --bernstein the hard hard prune
// above already rules out cells that can't contain a zero at all --
// what's missing is a way to rank the survivors, since box-width alone
// doesn't distinguish "wide because genuinely close to a fast-varying
// part of the curve" from "wide because this triangle is still large
// and the enclosure hasn't tightened yet". Reapplying the same
// proximity signal used for gradient dispersion is a second, cheap
// (already-cached per-vertex Fx/Fy/Fz/Fval) estimate of the latter, to
// bias the *ordering* among certified candidates toward the ones
// nearest the curve, on top of the hard prune (which is unaffected --
// this only ever scales an already-positive box-width priority).
double cell_diam(const T& t, const vgtl::array<Vertex(T),DIM+1>& vs) {
	pt3 pts[DIM+1];
	for(int i=0;i<=DIM;++i) pts[i]=attr(t,vs[i])->p;
	double diam=0;
	for(int i=0;i<=DIM;++i)
		for(int j=i+1;j<=DIM;++j)
			diam=std::max(diam,fs_dist(pts[i],pts[j]));
	return diam;
}

double proximity_factor(const T& t, const vgtl::array<Vertex(T),DIM+1>& vs) {
	double mind=1e300;
	for(int i=0;i<=DIM;++i) {
		pt3 g=vertex_gradient(t,vs[i]);
		double gn=hnorm(g);
		if(gn<1e-12) continue;
		double dv=std::abs(vertex_Fval(t,vs[i]))/gn;
		if(dv<mind) mind=dv;
	}
	double diam=cell_diam(t,vs);
	return diam/(mind+diam);
}

// Raw box width isn't comparable across levels: for a smooth curve it
// shrinks roughly LINEARLY with the cell's own diameter (first-order
// Taylor behavior), so comparing it against one fixed --threshold
// conflates "hasn't shrunk yet because the cell is still large" with
// "genuinely needs more refinement". Dividing by the cell's own
// diameter removes that size-dependent part, leaving a quantity on a
// consistent scale (units of Re/Im per unit Fubini-Study distance,
// i.e. a rate) regardless of level -- the same reasoning --proximity
// already used (comparing a length to the cell's own diameter) applied
// here to the box itself rather than to a separate distance estimate.
double cell_priority_bernstein(T& t, Cell(T) cv) {
	double best=-1.0;
	for(int j=0;j<=DIM;++j) {
		Simplex(T,3) tet=face_op(t,cv,j);
		for(int k=0;k<=3;++k) {
			Simplex(T,2) tri=face_op(t,tet,k);
			if(triangle_maybe_zero(t,tri)) {
				const extra_data<2>* d=attr(t,tri);
				double w=(d->reHi-d->reLo)+(d->imHi-d->imLo);
				if(w>best) best=w;
			}
		}
	}
	if(best>=0) {
		vgtl::array<Vertex(T),DIM+1> vs;
		vertices(t,cv,vs);
		double diam=cell_diam(t,vs);
		if(diam>1e-12) best/=diam;
		if(g_proximity) best*=proximity_factor(t,vs);
	}
	return best;
}

// --- Per-triangle local frame + Newton solve, unchanged from
// examples/top/riemann/riemann.cpp (see there for the full derivation
// comments) -- entirely generic in terms of "2 complex ambient
// coordinates", which here are a chart's 2 free affine coordinates
// (a,b) instead of the old (w,z).
double real_dot(cx a, cx b) { return a.real()*b.real()+a.imag()*b.imag(); }

struct tri_frame {
	int io, ip, iq;
	cx Pow, Poz;
	cx Xw, Xz, Yw, Yz;
	double xp;
};

void build_tri_frame(const cx w[3], const cx z[3], tri_frame& fr) {
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

void bary_to_xy(const tri_frame& fr, double l1, double l2, double& x, double& y) {
	double l[3]={1.0-l1-l2,l1,l2};
	y=l[fr.io];
	x=l[fr.ip]*fr.xp+l[fr.iq];
}

typedef cx (*cxfun2)(cx,cx);

// --- Error statistics for the Newton-solved triangle crossings: every
// accepted root's raw residual |F| (chart-scale-dependent, kept mainly
// for sanity) and a scale-invariant geometric estimate |F|/|dF/d(x,y)|
// -- a first-order estimate, in the SAME (x,y) local-frame units as
// build_tri_frame/bary_to_xy, of how far the converged (x,y) still sits
// from the true zero. Dividing by the local gradient cancels any
// overall scale of F (e.g. switching chart, or an arbitrary constant
// multiple in the catalog polynomial), which raw |F| does not.
vector<double> g_resid_F, g_resid_geom;
vector<int> g_resid_iters;

bool newton_on_triangle(const tri_frame& fr, double x, double y, double& out_l1, double& out_l2,
												 cxfun2 Ffun, cxfun2 Fwfun, cxfun2 Fzfun,
												 double& out_resid, double& out_geom_resid, int& out_iters) {
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

// --- The 3 symmetric affine charts of CP^2 (X!=0, Y!=0, Z!=0). Unlike
// the old P^1 x P^1 setup, there is no distinguished/pole chart: every
// point of CP^2 has at least one homogeneous coordinate with
// |coordinate| >= 1/sqrt(3) (since the representative is unit-norm), so
// picking whichever chart is best-conditioned for the 3 vertices of a
// given triangle is always possible and needs no special-casing.
cx F_chart0(cx a, cx b) { return eval_poly3(g_F, cx(1,0), a, b); }   // a=Y/X, b=Z/X
cx Fa_chart0(cx a, cx b){ return eval_poly3(g_Fy, cx(1,0), a, b); }
cx Fb_chart0(cx a, cx b){ return eval_poly3(g_Fz, cx(1,0), a, b); }
cx F_chart1(cx a, cx b) { return eval_poly3(g_F, a, cx(1,0), b); }   // a=X/Y, b=Z/Y
cx Fa_chart1(cx a, cx b){ return eval_poly3(g_Fx, a, cx(1,0), b); }
cx Fb_chart1(cx a, cx b){ return eval_poly3(g_Fz, a, cx(1,0), b); }
cx F_chart2(cx a, cx b) { return eval_poly3(g_F, a, b, cx(1,0)); }   // a=X/Z, b=Y/Z
cx Fa_chart2(cx a, cx b){ return eval_poly3(g_Fx, a, b, cx(1,0)); }
cx Fb_chart2(cx a, cx b){ return eval_poly3(g_Fy, a, b, cx(1,0)); }

void dehomogenize(int chart, const pt3& p, cx& a, cx& b) {
	switch(chart) {
		case 0: a=p[1]/p[0]; b=p[2]/p[0]; break;
		case 1: a=p[0]/p[1]; b=p[2]/p[1]; break;
		default: a=p[0]/p[2]; b=p[1]/p[2]; break;
	}
}
pt3 rehomogenize(int chart, cx a, cx b) {
	pt3 p;
	switch(chart) {
		case 0: p[0]=cx(1,0); p[1]=a; p[2]=b; break;
		case 1: p[0]=a; p[1]=cx(1,0); p[2]=b; break;
		default: p[0]=a; p[1]=b; p[2]=cx(1,0); break;
	}
	return normalize3(p);
}

void
triangle_intersection(T& t, Simplex(T,2) tri) {
	extra_data<2>* d=attr(t,tri);
	if(d->computed) return;
	d->computed=true;

	if(g_bernstein && !triangle_maybe_zero(t,tri)) {
		++g_bernstein_pruned_faces; // certified empty: Newton never had to run here
		return;
	}

	vgtl::array<Vertex(T),3> vs;
	vertices(t,tri,vs);
	pt3 p[3]; for(int k=0;k<3;++k) p[k]=point(t,vs[k]);

	int chart=pick_chart(p); // same rule compute_bernstein_bounds uses -- see pick_chart

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

	// Found roots stay in local temporaries -- exactly the old d->l1/
	// d->l2/d->pt fields, just not struct members any more -- until the
	// final count is known, then get pushed as one contiguous run into
	// g_roots below (see extra_data<2>::root_idx).
	double tmpl1[2],tmpl2[2]; pt3 tmppt[2]; int nfound=0;
	for(int s=0; s<nseeds && nfound<2; ++s) {
		double sx,sy; bary_to_xy(fr,seeds[s][0],seeds[s][1],sx,sy);
		double ol1,ol2,rF,rG; int riter;
		if(!newton_on_triangle(fr,sx,sy,ol1,ol2,Ffun,Fwfun,Fzfun,rF,rG,riter)) continue;
		bool dup=false;
		for(int r=0;r<nfound;++r)
			if(fabs(tmpl1[r]-ol1)<1e-7 && fabs(tmpl2[r]-ol2)<1e-7) dup=true;
		if(dup) continue;
		double ol0=1.0-ol1-ol2;
		int r=nfound;
		tmpl1[r]=ol1; tmpl2[r]=ol2;
		cx aroot=ol0*w0+ol1*w1+ol2*w2, broot=ol0*z0+ol1*z1+ol2*z2;
		tmppt[r]=rehomogenize(chart,aroot,broot);
		++nfound;
		g_resid_F.push_back(rF); g_resid_geom.push_back(rG); g_resid_iters.push_back(riter);
	}
	d->nroots=(unsigned char)nfound;
	if(nfound>0) {
		d->root_idx=(int)g_roots.size();
		for(int r=0;r<nfound;++r) {
			root_data rd; rd.pt=narrow_pt(tmppt[r]); rd.l1=tmpl1[r]; rd.l2=tmpl2[r];
			g_roots.push_back(rd);
		}
	}
}

int triangle_nroots(const T& t, Simplex(T,2) tri) { return attr(t,tri)->nroots; }
void triangle_point(const T& t, Simplex(T,2) tri, int r, pt3& p) {
	p=widen_pt(g_roots[attr(t,tri)->root_idx+r].pt);
}
void triangle_bary(const T& t, Simplex(T,2) tri, int r, double& l0, double& l1, double& l2) {
	const root_data& rd=g_roots[attr(t,tri)->root_idx+r];
	l1=rd.l1; l2=rd.l2; l0=1.0-l1-l2;
}

// --- Crossing-node identity and cycle extraction: unchanged in spirit
// from examples/top/riemann/riemann.cpp (see there for the full
// rationale), just carrying a single homogeneous point instead of (w,z).
struct crossing_node {
	int dim, desc, sub;
	double param;
	pt3 p;
};

void
compute_crossing_node(const T& t, Simplex(T,2) tri, int r, crossing_node& nd) {
	double l0,l1,l2;
	triangle_bary(t,tri,r,l0,l1,l2);
	double l[3]={l0,l1,l2};
	const double tol=1e-6;
	int zeros=0, zi[3];
	for(int b=0;b<3;++b) if(fabs(l[b])<tol) zi[zeros++]=b;
	triangle_point(t,tri,r,nd.p);
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

// --- Visualization projections -----------------------------------------
// Two modes, selected by --onion:
//
// - flat (default): pick, PER POINT (there's no single global chart
//   anymore -- every point independently uses whichever of its 3
//   homogeneous coordinates is largest), the best-conditioned affine
//   chart, and output (Re a, Im a, Re b). A projection, not an
//   embedding, exactly like examples/top/riemann's flat mode -- Im(b)
//   is dropped.
//
// - onion: a direct CP^2 analogue of examples/top/riemann's onion mode,
//   but simpler -- there, the "companion root" trick was needed because
//   w and z played structurally different roles (w the fiber variable
//   of an explicit branched cover, z the base) and had to be told apart
//   from F's own coefficients. Here every catalog curve is just a single
//   homogeneous F(X,Y,Z)=0 with no distinguished variable, so there's no
//   "other root" to compute -- direction and bump both come straight off
//   the point's own 3 homogeneous coordinates: Z (an arbitrary but fixed
//   choice, uniform across every curve in the catalog, unlike picking
//   "z" would have been before this project's whole point of dropping
//   distinguished directions) is stereographically projected to a
//   direction on S^2 via the same to_sphere() formula as before, and the
//   bump is |X|+|Y| -- the sum of the moduli of the OTHER two
//   homogeneous coordinates. Because every point is kept unit-normalized
//   (|X|^2+|Y|^2+|Z|^2=1), |X|+|Y| is automatically bounded (<=sqrt(2)),
//   so the bump stays well-behaved with no extra clamping.
double g_cutoff=1e300;
bool g_onion=false;
double g_onion_scale=0.5;

vec<3,double> to_sphere(cx w) {
	double n=std::norm(w);
	double s=n+1.0;
	vec<3,double> p;
	p[0]=2.0*w.real()/s; p[1]=2.0*w.imag()/s; p[2]=(n-1.0)/s;
	return p;
}

// Picks ONE chart for an entire output polygon (worst-case-best rule,
// same shape as pick_chart/triangle_intersection's own choice, but
// maximizing the minimum |coordinate| across ALL of the polygon's
// nodes, not just one triangle's 3 vertices). Fixes a real bug: the
// previous project_for_viz picked a chart independently PER POINT, so
// a polygon whose nodes came from different triangles/tetrahedra could
// get its vertices dehomogenized in DIFFERENT (a,b) frames -- each
// point individually correct, but stitched together into a polygon
// spanning two incompatible coordinate systems. Confirmed empirically
// before this fix (not just suspected): 1206/173725 polygons (0.69%)
// had a per-point chart mismatch, matching almost exactly the 3106
// edges (0.52%, ~2.5 per mismatched polygon) sitting 10x+ longer than
// the median in the exported mesh -- the source of the "ribbon spray"
// artifact seen in renders despite a low bad-cell rate.
int pick_chart_polygon(const vector<crossing_node>& nodes, const vector<int>& cyc) {
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

vec<3,double> project_for_viz(const pt3& p, int chart) {
	if(g_onion) {
		vec<3,double> dir=to_sphere(p[2]); // Z: direction on S^2
		double bump=std::abs(p[0])+std::abs(p[1]); // |X|+|Y|: radial displacement
		double r=1.0+g_onion_scale*bump;
		vec<3,double> out;
		for(int i=0;i<3;++i) out[i]=r*dir[i];
		return out;
	}
	cx a,b; dehomogenize(chart,p,a,b);
	vec<3,double> r;
	r[0]=a.real(); r[1]=a.imag(); r[2]=b.real();
	return r;
}

struct obj_writer {
	ofstream out;
	int nverts, nfaces, nclipped;
	obj_writer(const char* path) : out(path), nverts(0), nfaces(0), nclipped(0) {
		if(g_onion)
			out<<"# Riemann surface (CP^2 approach) extraction, onion projection: "
					<<"direction=Z on S^2, radius=1+"<<g_onion_scale<<"*(|X|+|Y|)\n";
		else
			out<<"# Riemann surface (CP^2 approach) extraction, projected via "
					<<"per-point best affine chart (Re a, Im a, Re b)\n";
		if(!g_onion && g_cutoff<1e299)
			out<<"# cutoff: polygons with a vertex farther than "<<g_cutoff<<" from the origin dropped\n";
	}
	void write_polygon(const vector<crossing_node>& nodes, const vector<int>& cyc) {
		int chart=g_onion ? 0 : pick_chart_polygon(nodes,cyc); // one chart for the whole polygon
		vector<vec<3,double> > pts(cyc.size());
		for(size_t i=0;i<cyc.size();++i) pts[i]=project_for_viz(nodes[cyc[i]].p,chart);
		if(!g_onion) {
			for(size_t i=0;i<pts.size();++i) {
				double dd=sqrt(pts[i][0]*pts[i][0]+pts[i][1]*pts[i][1]+pts[i][2]*pts[i][2]);
				if(dd>g_cutoff) { ++nclipped; return; }
			}
		}
		for(size_t i=0;i<cyc.size();++i)
			out<<"v "<<pts[i][0]<<" "<<pts[i][1]<<" "<<pts[i][2]<<"\n";
		out<<"f";
		for(size_t i=0;i<cyc.size();++i) out<<" "<<(nverts+(int)i+1);
		out<<"\n";
		nverts+=(int)cyc.size();
		++nfaces;
	}
};

int main(int argc, char* argv[]) {

	int max_depth=10;
	double threshold=0.1;
	int function_index=1; // fermat_cubic by default
	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--depth" && i+1<argc) {
			max_depth=atoi(argv[++i]);
		} else if(arg=="--threshold" && i+1<argc) {
			threshold=atof(argv[++i]);
		} else if(arg=="--function" && i+1<argc) {
			function_index=atoi(argv[++i]);
		} else if(arg=="--cutoff" && i+1<argc) {
			g_cutoff=atof(argv[++i]);
		} else if(arg=="--onion") {
			g_onion=true;
		} else if(arg=="--onion-scale" && i+1<argc) {
			g_onion_scale=atof(argv[++i]);
		} else if(arg=="--proximity") {
			g_proximity=true;
		} else if(arg=="--generic") {
			g_generic=true;
		} else if(arg=="--generic-seed" && i+1<argc) {
			g_generic=true; g_generic_seed=(unsigned)atoi(argv[++i]);
		} else if(arg=="--bernstein") {
			g_bernstein=true;
		} else if(arg=="--bernstein-level" && i+1<argc) {
			g_bernstein=true; g_bernstein_level=atoi(argv[++i]);
		} else if(arg=="--bernstein-selftest") {
			g_bernstein_selftest=true;
		} else if(arg=="--list-functions") {
			cout<<"available functions:"<<endl;
			print_function_catalog_pc2(cout);
			return 0;
		} else {
			cerr<<"unrecognized argument: "<<arg<<endl;
			cerr<<"usage: "<<argv[0]<<" [--function N] [--depth N] [--threshold X] "
					<<"[--cutoff X] [--onion] [--onion-scale X] [--proximity] "
					<<"[--generic] [--generic-seed N] [--bernstein] [--bernstein-level N] "
					<<"[--list-functions]"<<endl;
			return 1;
		}
	}

	vector<catalog_entry_pc2>& catalog=function_catalog_pc2();
	if(function_index<0 || function_index>=(int)catalog.size()) {
		cerr<<"--function "<<function_index<<" out of range; available:"<<endl;
		print_function_catalog_pc2(cerr);
		return 1;
	}
	set_curve(catalog[function_index].F);

	if(g_bernstein_selftest) {
		bernstein_selftest();
		return 0;
	}

	cout<<"function: "<<function_index<<" ("<<catalog[function_index].name<<") -- "
			<<catalog[function_index].description<<endl;
	cout<<"degree: "<<g_d<<"  max_depth="<<max_depth<<" threshold="<<threshold
			<<"  proximity="<<(g_proximity?"on":"off")
			<<"  generic="<<(g_generic?"on":"off");
	if(g_generic) cout<<" (seed="<<g_generic_seed<<")";
	cout<<"  bernstein="<<(g_bernstein?"on":"off");
	if(g_bernstein) cout<<" (level="<<g_bernstein_level<<")";
	cout<<endl;
	cout<<"projection: "<<(g_onion?"onion":"flat");
	if(g_onion) cout<<" (scale="<<g_onion_scale<<")";
	cout<<endl;
	if(!g_onion && g_cutoff<1e299) cout<<"cutoff: "<<g_cutoff<<endl;

	string obj_path_s="riemann_pc2_"+catalog[function_index].name+(g_onion?"_onion":"")+".obj";
	const char* obj_path=obj_path_s.c_str();

	// --- Phase 1: seed mesh -- Gaifullin's 15-vertex/108-cell
	// triangulation of CP^2, already Maubach-compatible and orientable
	// with no subdivision step (see gaifullin_points()/gaifullin_cells()
	// above for why this replaced Kuhnel's 9-vertex/36-cell + barycentric
	// subdivision).
	cout<<endl<<"--- seed mesh: Gaifullin CP^2_15 ---"<<endl;

	vector<pt3> gp=gaifullin_points();
	vector<Vertex(T)> V(15);
	for(int i=0;i<15;++i) {
		V[i]=add(t);
		point_set(t,V[i],gp[i]);
	}
	vector<cell5> gc=gaifullin_cells();
	cout<<"Gaifullin base cells: "<<gc.size()<<" (expected 108)"<<endl;
	{
		complex_builder<T> cb(t);
		for(size_t ci=0; ci<gc.size(); ++ci) {
			vgtl::array<Vertex(T),5> vs;
			for(int k=0;k<5;++k) vs[k]=V[gc[ci][k]];
			Cell(T) cv=add(cb,vs);
			orientation_set(t,cv,+1); // arbitrary seed -- overwritten by the BFS below
			level_set(t,cv,0);
		}
	}

	int nv=0,n4=0;
	{ Vertex_it(T) i,end; for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++nv; }
	{ Cell_it(T) i,end; for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++n4; }
	cout<<"vertices: "<<nv<<" (expected 15)"<<endl;
	cout<<"4-simplices: "<<n4<<" (expected 108)"<<endl;
	cout<<"euler characteristic: "<<euler_characteristic(t)<<" (expected 3, = CP^2's)"<<endl;

	int nbnd=0;
	{ Simplex_it(T,3) i,end; for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i) && boundary(t,*i)) ++nbnd; }
	cout<<"boundary 3-simplices: "<<nbnd<<" (expected 0)"<<endl;

	cout<<"incoherent facets from ord_split's propagated orientation: "<<count_incoherent(t)<<endl;
	reorient_via_bfs(t,n4);
	cout<<"incoherent facets after BFS re-derivation: "<<count_incoherent(t)<<" (expected 0)"<<endl;

	// --- Phase 2: adaptive refinement (gradient-dispersion cell_priority,
	// or --bernstein's certified enclosure test -- see cell_priority vs
	// cell_priority_bernstein above) --
	cout<<endl<<"--- adaptive refinement ---"<<endl;

	priority_queue<pair<double,Cell(T)> > pq;
	int ndropped=0;
	{
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend); ci!=cend; ++ci) {
			if(!is_current(t,*ci)) continue;
			double p=g_bernstein ? cell_priority_bernstein(t,*ci) : cell_priority(t,*ci);
			if(g_bernstein && p<0) { ++ndropped; continue; } // certified empty
			pq.push(make_pair(p,*ci));
		}
	}

	int nsubdivisions=0;
	while(!pq.empty()) {
		pair<double,Cell(T)> top=pq.top(); pq.pop();
		Cell(T) cv=top.second;
		if(!is_current(t,cv)) continue;
		if(top.first<threshold) break;
		if(level(t,cv)>=max_depth) continue;
		refine_app app;
		maubach_subdivide(t,cv,app);
		++nsubdivisions;
		for(size_t i=0;i<app.new_cells.size();++i) {
			Cell(T) nc=app.new_cells[i];
			double p=g_bernstein ? cell_priority_bernstein(t,nc) : cell_priority(t,nc);
			if(g_bernstein && p<0) { ++ndropped; continue; }
			pq.push(make_pair(p,nc));
		}
	}
	cout<<"subdivisions performed: "<<nsubdivisions<<endl;
	if(g_bernstein) cout<<"cells certified empty (dropped, never refined): "<<ndropped<<endl;

	int n42=0;
	{ Cell_it(T) i,end; for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++n42; }
	cout<<"4-simplices after refinement: "<<n42<<endl;
	cout<<"euler characteristic after refinement: "<<euler_characteristic(t)<<" (expected 3)"<<endl;
	{
		map<int,int> hist;
		Cell_it(T) i,end;
		for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++hist[level(t,*i)];
		cout<<"final cell level histogram:";
		for(map<int,int>::iterator hi=hist.begin(); hi!=hist.end(); ++hi)
			cout<<" ["<<hi->first<<"]="<<hi->second;
		cout<<endl;
	}

	// --- Phase 3: extraction --------------------------------------------
	cout<<endl<<"--- surface extraction ---"<<endl;

	int npoly_out[8]={0,0,0,0,0,0,0,0};
	int ntouching_tets=0, nbad_tets=0, ok_cells=0, bad_cells=0;

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
					static const int pairings[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
					int best=0; double bestcost=1e300;
					for(int c=0;c<3;++c) {
						double cost=0;
						for(int e=0;e<2;++e) {
							const crossing_node& A=nodes[tnodes[pairings[c][2*e]]];
							const crossing_node& B=nodes[tnodes[pairings[c][2*e+1]]];
							cost+=fs_dist(A.p,B.p);
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
						int nxt=(adj[cur][0]==prev) ? adj[cur][1] : adj[cur][0];
						prev=cur; cur=nxt;
						if(cyc.size()>(size_t)(DIM+1)) { decompose_ok=false; break; }
					}
					if(!decompose_ok) break;
					if(cyc.size()<3) { decompose_ok=false; break; }
					cycles_idx.push_back(cyc);
				}
			}

			if(!decompose_ok) { ++bad_cells; continue; }
			++ok_cells;

			for(size_t p=0;p<cycles_idx.size();++p) {
				int sz=(int)cycles_idx[p].size();
				if(sz>=3 && sz<8) ++npoly_out[sz];
				obj.write_polygon(nodes,cycles_idx[p]);
			}
		}
	}

	cout<<"extracted cells: ok="<<ok_cells<<" bad="<<bad_cells<<endl;
	cout<<"polygon sizes:";
	for(int sz=3; sz<8; ++sz) if(npoly_out[sz]) cout<<" "<<sz<<"-gon="<<npoly_out[sz];
	cout<<endl;
	cout<<"touching tetrahedra: "<<ntouching_tets<<"  unhandled node count: "<<nbad_tets<<endl;
	cout<<"wrote "<<obj_path<<": "<<obj.nverts<<" vertices, "<<obj.nfaces<<" faces"<<endl;
	if(obj.nclipped) cout<<"polygons dropped by --cutoff: "<<obj.nclipped<<endl;
	if(g_bernstein) cout<<"triangles Bernstein-pruned (Newton skipped, certified no root): "
			<<g_bernstein_pruned_faces<<endl;

	// Diagnostic: how many 2-simplices actually carry 2 roots (the
	// nroots==2 case extra_data<2>::pt[2]/l1[2]/l2[2] exist for), vs 1 or
	// 0 -- triangle_intersection() is idempotent (see d->computed), so
	// this re-walk is free, every triangle here was already resolved
	// during extraction above.
	{
		long hist[3]={0,0,0};
		Simplex_it(T,2) si,send;
		for(simplices(t,si,send); si!=send; ++si) {
			if(!is_current(t,*si)) continue;
			triangle_intersection(t,*si);
			++hist[triangle_nroots(t,*si)];
		}
		cout<<"triangle nroots histogram: 0="<<hist[0]<<" 1="<<hist[1]<<" 2="<<hist[2]<<endl;
	}

	// --- Solution-error statistics: |F| residual and the scale-invariant
	// geometric estimate |F|/|dF/d(x,y)| (see newton_on_triangle) over
	// every accepted triangle-crossing root.
	if(!g_resid_F.empty()) {
		cout<<endl<<"--- solution error statistics ("<<g_resid_F.size()<<" accepted roots) ---"<<endl;
		vector<double> sF=g_resid_F, sG=g_resid_geom;
		std::sort(sF.begin(),sF.end());
		std::sort(sG.begin(),sG.end());
		auto pct=[&](const vector<double>& v, double q) {
			size_t idx=(size_t)(q*(v.size()-1));
			return v[idx];
		};
		double sumF=0,sumG=0; long long sumIter=0;
		for(size_t i=0;i<sF.size();++i) sumF+=g_resid_F[i];
		for(size_t i=0;i<sG.size();++i) sumG+=g_resid_geom[i];
		for(size_t i=0;i<g_resid_iters.size();++i) sumIter+=g_resid_iters[i];
		cout<<"raw |F| residual:      min="<<sF.front()<<"  median="<<pct(sF,0.5)
				<<"  mean="<<(sumF/sF.size())<<"  p95="<<pct(sF,0.95)<<"  max="<<sF.back()<<endl;
		cout<<"geometric |F|/|dF|:    min="<<sG.front()<<"  median="<<pct(sG,0.5)
				<<"  mean="<<(sumG/sG.size())<<"  p95="<<pct(sG,0.95)<<"  max="<<sG.back()<<endl;
		cout<<"newton iterations:     mean="<<(double(sumIter)/g_resid_iters.size())<<endl;
	}

	return 0;
}
