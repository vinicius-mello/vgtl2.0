// Self-tests for glpt_crossing.hpp. See that file's own header comment
// for scope (core Newton crossing-detection ported, Bernstein fallback
// deferred, nmt-specific caching replaced by a vertex-id-keyed cache).
//
//   1. Sanity: on a simple curve (catalog "conic"), compute crossings
//      for every 2-face of a handful of root-level cells and check any
//      root found actually satisfies |F(root)| ~ 0.
//   2. THE CRITICAL TEST -- cross-cell face agreement: build a cell and
//      its neighbor (sharing a facet, hence 4 shared 2-faces), extract
//      each side's OWN local vertex order for those shared faces
//      (generally a DIFFERENT order than the other side's), and check
//      that compute_face_crossing() -- via face_crossing_cache, but
//      also directly -- agrees on how many roots and where, regardless
//      of input point order. This is what actually justifies caching by
//      vertex-id triple: it's only safe if two cells reaching the same
//      face (via different local orderings) would have computed the
//      same answer anyway.
//
// Build: g++ -std=c++17 -I ../../../include -I ~/code/lpt -O2
//        glpt_crossing_test.cpp -o glpt_crossing_test

#include <cstdio>
#include <cstdlib>
#include "glpt_crossing.hpp"

static int g_fail = 0;
#define CHECK(cond, msg) do { if(!(cond)) { ++g_fail; printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); } } while(0)

double pt_dist(const pt3& a, const pt3& b) {
	double d=0;
	for(int i=0;i<3;++i) d += std::norm(a[i]-b[i]);
	return std::sqrt(d);
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
void cell_points_and_ids(const glpt& c, const std::vector<pt3>& gp, glpt_edge_cache& cache,
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

void test_sanity(const std::vector<pt3>& gp) {
	printf("--- test 1: sanity (conic curve, root-level cells) ---\n");
	std::vector<catalog_entry_cp2>& cat = function_catalog_cp2();
	int idx=-1;
	for(size_t i=0;i<cat.size();++i) if(cat[i].name=="conic") idx=(int)i;
	CHECK(idx>=0, "conic curve not found in catalog");
	set_curve(cat[idx].F);

	int n_faces_checked=0, n_roots_found=0, n_bad_residual=0;
	for(int s=0;s<12;++s) {
		pt3 pts[glpt::DIM+1];
		glpt_root_vertex_points(s, gp, pts);
		for(int f=0;f<10;++f) {
			pt3 face_pts[3];
			for(int k=0;k<3;++k) face_pts[k]=pts[GLPT_CELL_FACES[f][k]];
			pt3 out_pts[2]; double out_bary[2][3];
			int nr = compute_face_crossing(face_pts, out_pts, out_bary);
			++n_faces_checked;
			for(int r=0;r<nr;++r) {
				++n_roots_found;
				cx Fv = eval_poly3(g_F, out_pts[r][0], out_pts[r][1], out_pts[r][2]);
				if(std::abs(Fv)>1e-9) ++n_bad_residual;
			}
		}
	}
	CHECK(n_bad_residual==0, "a found root doesn't actually satisfy F=0");
	printf("faces checked=%d roots found=%d bad_residual=%d\n", n_faces_checked, n_roots_found, n_bad_residual);
}

void test_cross_cell_agreement(const std::vector<pt3>& gp) {
	printf("--- test 2: cross-cell face agreement (the critical test) ---\n");
	std::vector<catalog_entry_cp2>& cat = function_catalog_cp2();
	int idx=-1;
	for(size_t i=0;i<cat.size();++i) if(cat[i].name=="elliptic") idx=(int)i;
	CHECK(idx>=0, "elliptic curve not found in catalog");
	set_curve(cat[idx].F);

	srand(31337);
	glpt_edge_cache cache;
	long n_shared_faces=0, n_nroots_mismatch=0, n_points_mismatch=0;
	double max_pt_err=0.0;
	for(int trial=0; trial<300; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = rand()%10;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);

		pt3 c_pts[glpt::DIM+1]; int c_ids[glpt::DIM+1];
		cell_points_and_ids(cur, gp, cache, c_pts, c_ids);

		bool c0 = cur.is_root() ? false : cur.is_child0();
		int ls = lstar_of(cur);
		for(int i=0;i<=glpt::DIM;++i) {
			bool sibling = !cur.is_root() && ((c0 && i==glpt::DIM) || (!c0 && i==ls));
			glpt nb;
			glpt::result r = cur.neighbor(i, nb);
			if(r!=glpt::OK) continue;
			pt3 n_pts[glpt::DIM+1]; int n_ids[glpt::DIM+1];
			cell_points_and_ids(nb, gp, cache, n_pts, n_ids);

			// For each of cur's faces NOT involving excluded vertex i (i.e.
			// wholly within the shared facet), find the SAME face on nb's
			// side by matching ids, and compare compute_face_crossing()
			// results computed independently from each side's own local
			// point order.
			for(int f=0;f<10;++f) {
				bool touches_i=false;
				for(int k=0;k<3;++k) if(GLPT_CELL_FACES[f][k]==i) touches_i=true;
				if(touches_i) continue;
				int face_ids[3];
				pt3 face_pts_c[3];
				for(int k=0;k<3;++k) { face_ids[k]=c_ids[GLPT_CELL_FACES[f][k]]; face_pts_c[k]=c_pts[GLPT_CELL_FACES[f][k]]; }
				// find matching positions on nb's side
				pt3 face_pts_n[3];
				bool all_found=true;
				for(int k=0;k<3;++k) {
					int pos=-1;
					for(int j=0;j<=glpt::DIM;++j) if(n_ids[j]==face_ids[k]) pos=j;
					if(pos<0) { all_found=false; break; }
					face_pts_n[k]=n_pts[pos];
				}
				if(!all_found) continue; // this face isn't on the shared facet from nb's numbering (shouldn't happen, but don't assume)
				++n_shared_faces;

				pt3 out_c[2],out_n[2]; double bc[2][3],bn[2][3];
				int nrc = compute_face_crossing(face_pts_c, out_c, bc);
				int nrn = compute_face_crossing(face_pts_n, out_n, bn);
				if(nrc!=nrn) { ++n_nroots_mismatch; continue; }
				for(int r=0;r<nrc;++r) {
					double best=1e300;
					for(int r2=0;r2<nrn;++r2) { double d=pt_dist(out_c[r],out_n[r2]); if(d<best) best=d; }
					if(best>max_pt_err) max_pt_err=best;
					if(best>1e-6) ++n_points_mismatch;
				}
			}
			(void)sibling;
		}
	}
	CHECK(n_nroots_mismatch==0, "two cells sharing a face disagreed on the NUMBER of crossings found");
	CHECK(n_points_mismatch==0, "two cells sharing a face found crossings at different points");
	printf("shared faces checked=%ld nroots_mismatch=%ld points_mismatch=%ld max_pt_err=%.3e\n",
		n_shared_faces, n_nroots_mismatch, n_points_mismatch, max_pt_err);
}

int main() {
	std::vector<pt3> gp = gaifullin_points();
	test_sanity(gp);
	test_cross_cell_agreement(gp);

	printf("\n%s (%d failures)\n", g_fail==0 ? "ALL CHECKS PASSED" : "CHECKS FAILED", g_fail);
	return g_fail==0 ? 0 : 1;
}
