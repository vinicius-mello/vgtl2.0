#ifndef GLPT_POINTS_HPP
#define GLPT_POINTS_HPP

/*! \file
 * \brief Real (Fubini-Study) coordinates for glpt cells' vertices --
 * the first missing piece for a glpt_tree-based CP^2 mesh (see the
 * conversation that led here: glpt.hpp/glpt_vertex_ids.hpp/glpt_tree.hpp
 * in ~/code/lpt only give combinatorics; nothing there yet produces an
 * actual point in CP^2).
 *
 * ============================================================
 * WHY THIS ISN'T JUST glpt::vertex_weights() DOT gaifullin_points()
 * ============================================================
 * vertex_weights() gives a cell's vertices as REAL, LINEAR barycentric
 * weights relative to its seed's 5 original points -- exact and
 * closed-form (Theorem 1). It would be natural to assume a deep cell's
 * actual embedded vertex is just that linear combination of the seed's
 * 5 real (complex-projective) points. It is NOT: riemann_cp2.cpp's own
 * Maubach refinement (refine_app::add_edge_vertex, riemann_cp2.cpp)
 * places every NEW edge-bisection vertex at the Fubini-Study GEODESIC
 * midpoint of its two immediate parents -- align_phase (resolve the
 * U(1) representative-phase ambiguity between the two homogeneous
 * points) then normalize3(sum) -- not a linear average, and NOT
 * expressible as a single closed-form linear combination of the
 * ORIGINAL 5 (or 15) points once you're more than one bisection deep
 * (the geodesic/spherical combination doesn't distribute like a linear
 * one does). So a cell's real coordinates are, like glpt_vertex_ids.hpp's
 * vertex ids, inherently STATEFUL: they must be threaded through the
 * actual bisection PATH (root -> child -> child -> ...), recomputing
 * each new midpoint from its own two immediate parents in turn, not
 * derived from a bare glpt code in isolation.
 *
 * This file mirrors glpt_vertex_ids.hpp's own child-vertex matching
 * approach (compare a child's exact barycentric weights against its
 * parent's to find which DIM positions are inherited unchanged and
 * which ONE is new -- see vertex_weights_exact()'s own comment in
 * glpt.hpp) to identify the new position, then computes ITS point via
 * align_phase+normalize3 on the two parent points at local positions
 * level() and DIM -- exactly what add_edge_vertex does, just walking
 * glpt's own combinatorics instead of vgtl's nmt<4>/Edge(T).
 *
 * Lives here, not in ~/code/lpt, because it depends on CP^2-specific
 * geometry (pt3, the Fubini-Study inner product, gaifullin_points())
 * that has no place in a dimension/curve-agnostic LPT library -- see
 * glpt_tree.hpp's own header comment for the same reasoning applied to
 * why glpt_tree isn't a port of lpt_tree.
 *
 * gaifullin_points()/pt3/mkpt/hdot/hnorm/normalize3/align_phase/
 * fs_barycenter below are copied verbatim from riemann_cp2.cpp (not
 * included from it -- that file is a monolithic .cpp, not a reusable
 * header) -- the same precedent glpt.hpp itself already set by copying
 * gaifullin_cells_raw verbatim rather than depending on that file.
 */

#include <cmath>
#include <complex>
#include <vector>
#include <vgtl/alg/vec.hpp>
#include "glpt.hpp"

typedef std::complex<double> cx;
typedef vgtl::vec<3,cx> pt3;

inline pt3 mkpt(cx a, cx b, cx c) { pt3 p; p[0]=a; p[1]=b; p[2]=c; return p; }

inline cx hdot(const pt3& a, const pt3& b) {
	cx r(0,0);
	for(int i=0;i<3;++i) r += a[i]*std::conj(b[i]);
	return r;
}
inline double hnorm(const pt3& a) { return std::sqrt(std::max(0.0, hdot(a,a).real())); }
inline pt3 normalize3(pt3 a) {
	double n=hnorm(a);
	for(int i=0;i<3;++i) a[i]/=n;
	return a;
}
// Rephase b so <a,b> is real and non-negative -- resolves the U(1)
// ambiguity in b's homogeneous representative relative to a, needed
// before b can be meaningfully averaged with a.
inline pt3 align_phase(const pt3& a, pt3 b) {
	cx ip=hdot(a,b);
	double m=std::abs(ip);
	if(m<1e-14) return b; // near-orthogonal: no sensible phase to lock onto
	cx phase=ip/m;
	for(int i=0;i<3;++i) b[i]*=phase;
	return b;
}
// Geodesic (Fubini-Study) midpoint of two points -- the k=1 case of
// riemann_cp2.cpp's own fs_barycenter, which is all a single edge
// bisection ever needs.
inline pt3 fs_midpoint(const pt3& a, const pt3& b) {
	pt3 ab = align_phase(a,b);
	pt3 acc = a;
	for(int c=0;c<3;++c) acc[c]+=ab[c];
	return normalize3(acc);
}

// Gaifullin's 15 original points (arXiv:0904.4222, Construction 1.1 +
// Sec. 3), verbatim from riemann_cp2.cpp's gaifullin_points() -- fixed
// ids 0..14, matching glpt_vertex_ids.hpp's GLPT_BASE_VERTEX_COUNT and
// glpt_gaifullin_cells' own indexing.
inline std::vector<pt3> gaifullin_points() {
	double s3=std::sqrt(3.0);
	cx one(1,0), zero(0,0);
	cx omega(-0.5, s3/2.0), omega2(-0.5,-s3/2.0);
	std::vector<pt3> P(15);
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
	return P;
}

//! Root cell's vertex points: direct lookup via glpt_gaifullin_cells,
//! into the SAME gaifullin_points() every riemann_cp2.cpp curve uses.
inline void glpt_root_vertex_points(int seed, const std::vector<pt3>& gp, pt3 points[glpt::DIM+1]) {
	for(int k=0;k<=glpt::DIM;++k) points[k] = gp[glpt_gaifullin_cells[seed][k]];
}

//! Child (zo=0 or 1) vertex points, given the PARENT's (already-known)
//! points. See this file's own header comment for why the DIM
//! unchanged positions are found by exact-weight matching (same
//! approach as glpt_vertex_ids.hpp's glpt_child_vertex_ids(), reusing
//! glpt.hpp's own vertex_weights_exact()) rather than a hand-derived
//! permutation rule, and why the new one is a geodesic (not linear)
//! midpoint.
inline void glpt_child_vertex_points(const glpt& parent, const pt3 parent_points[glpt::DIM+1],
		int zo, pt3 child_points[glpt::DIM+1]) {
	const int DIM = glpt::DIM;
	glpt child = parent.child(zo);

	int pw[DIM+1][DIM+1], cw[DIM+1][DIM+1];
	int p_shift, c_shift;
	parent.vertex_weights_exact(pw, p_shift);
	child.vertex_weights_exact(cw, c_shift);
	int scale = 1<<(c_shift-p_shift);

	int n_unmatched=0;
	for(int k=0;k<=DIM;++k) {
		int match_j=-1;
		for(int j=0;j<=DIM;++j) {
			bool eq=true;
			for(int c=0;c<=DIM && eq;++c) if(cw[k][c] != pw[j][c]*scale) eq=false;
			if(eq) { match_j=j; break; }
		}
		if(match_j>=0) {
			child_points[k] = parent_points[match_j];
		} else {
			++n_unmatched;
			child_points[k] = fs_midpoint(parent_points[parent.level()], parent_points[DIM]);
		}
	}
	assert(n_unmatched==1 && "glpt_child_vertex_points: expected exactly one new vertex per bisection");
	(void)n_unmatched;
}

#endif // GLPT_POINTS_HPP
