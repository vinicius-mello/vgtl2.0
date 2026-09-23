// Self-tests for glpt_points.hpp. See that file's own header comment
// for the design rationale (why vertex coordinates must be threaded
// through the actual bisection path, geodesic not linear).
//
//   1. Root sanity: a root's vertex points match gaifullin_points()
//      directly, in glpt_gaifullin_cells' own order.
//   2. THE CRITICAL TEST -- neighbor point consistency: build a cell C
//      via a random bisection path, compute its points; compute its
//      neighbor N = C.neighbor(i); independently RECONSTRUCT N's own
//      points from scratch (its own root + its own path, reconstructed
//      via parent()/is_child0() exactly as
//      glpt_vertex_ids_test.cpp's path_from_root() does); check the
//      shared facet's points agree numerically on both sides. Unlike
//      glpt_vertex_ids.hpp's ids (exact integers), points are floating-
//      point results of an iterative (align_phase+normalize) geodesic
//      construction -- this checks the two independent reconstructions
//      of the SAME physical vertex actually agree, not just that they
//      SHOULD in principle.
//
// Build: g++ -std=c++17 -I ../../../include -I ~/code/lpt -O2
//        glpt_points_test.cpp -o glpt_points_test

#include <cstdio>
#include <cstdlib>
#include "glpt_points.hpp"

static int g_fail = 0;
#define CHECK(cond, msg) do { if(!(cond)) { ++g_fail; printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); } } while(0)

double pt_dist(const pt3& a, const pt3& b) {
	double d=0;
	for(int i=0;i<3;++i) d += std::norm(a[i]-b[i]);
	return std::sqrt(d);
}

void test_root_sanity(const std::vector<pt3>& gp) {
	printf("--- test 1: root vertex points match gaifullin_points() ---\n");
	int n_checked=0;
	for(int s=0;s<GLPT_NCELLS;++s) {
		pt3 pts[glpt::DIM+1];
		glpt_root_vertex_points(s, gp, pts);
		for(int k=0;k<=glpt::DIM;++k) {
			pt3 expected = gp[glpt_gaifullin_cells[s][k]];
			CHECK(pt_dist(pts[k],expected)<1e-12, "root vertex point doesn't match gaifullin_points()");
			++n_checked;
		}
	}
	printf("checked %d root vertex points\n", n_checked);
}

int lstar_of(const glpt& g) {
	int lvl = g.level();
	int lminus = (lvl==0) ? (glpt::DIM-1) : (lvl-1);
	return lminus+1;
}
void path_from_root(const glpt& c, int path[], int& path_len) {
	path_len = 0; glpt cur = c;
	while(!cur.is_root()) { path[path_len++] = cur.is_child0()?0:1; cur = cur.parent(); }
	for(int i=0;i<path_len/2;++i) { int t=path[i]; path[i]=path[path_len-1-i]; path[path_len-1-i]=t; }
}
void points_from_scratch(const glpt& c, const std::vector<pt3>& gp, pt3 points[glpt::DIM+1]) {
	int path[64], path_len;
	path_from_root(c, path, path_len);
	glpt_root_vertex_points(c.seed_index(), gp, points);
	glpt cur(c.seed_index());
	for(int i=0;i<path_len;++i) {
		pt3 child_points[glpt::DIM+1];
		glpt_child_vertex_points(cur, points, path[i], child_points);
		cur = cur.child(path[i]);
		for(int k=0;k<=glpt::DIM;++k) points[k]=child_points[k];
	}
}

void test_neighbor_point_consistency(const std::vector<pt3>& gp) {
	printf("--- test 2: neighbor point consistency (the critical test) ---\n");
	srand(24242);
	long n_checked=0, n_deep_skipped=0;
	double max_err=0.0;
	for(int trial=0; trial<20000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = rand()%20;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);

		pt3 c_pts[glpt::DIM+1];
		points_from_scratch(cur, gp, c_pts);

		bool c0 = cur.is_root() ? false : cur.is_child0();
		int ls = lstar_of(cur);
		for(int i=0;i<=glpt::DIM;++i) {
			bool sibling = !cur.is_root() && ((c0 && i==glpt::DIM) || (!c0 && i==ls));
			glpt nb;
			glpt::result r = cur.neighbor(i, nb);
			if(r!=glpt::OK) { ++n_deep_skipped; continue; }
			pt3 n_pts[glpt::DIM+1];
			points_from_scratch(nb, gp, n_pts);

			// The shared facet's points (all of C's except position i)
			// must appear, up to numerical tolerance, among N's points.
			for(int k=0;k<=glpt::DIM;++k) {
				if(k==i) continue;
				double best=1e300;
				for(int j=0;j<=glpt::DIM;++j) {
					double d=pt_dist(c_pts[k], n_pts[j]);
					if(d<best) best=d;
				}
				if(best>max_err) max_err=best;
				CHECK(best<1e-9, "a shared-facet vertex point doesn't match on the neighbor's side");
				++n_checked;
			}
			(void)sibling;
		}
	}
	printf("checked %ld shared-vertex point pairs, %ld unexpectedly needed deep-crossing, max mismatch=%.3e\n",
		n_checked, n_deep_skipped, max_err);
}

int main() {
	std::vector<pt3> gp = gaifullin_points();
	test_root_sanity(gp);
	test_neighbor_point_consistency(gp);

	printf("\n%s (%d failures)\n", g_fail==0 ? "ALL CHECKS PASSED" : "CHECKS FAILED", g_fail);
	return g_fail==0 ? 0 : 1;
}
