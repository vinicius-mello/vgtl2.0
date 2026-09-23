// riemann_cp2_glpt.cpp -- glpt_tree-based rebuild of riemann_cp2.cpp's
// pipeline (see glpt_points.hpp/glpt_crossing.hpp's own header comments
// for the two pieces this builds on, and the conversation that led to
// this file for the full rationale: an O(1)-per-current-cell mesh,
// ~/code/lpt/glpt_tree.hpp, replacing riemann_cp2.cpp's append-only
// nmt<4>-backed one).
//
// SCOPE OF THIS FIRST VERSION: seeds all 108 Gaifullin cells, refines by
// the SAME gradient-dispersion cell_priority as riemann_cp2.cpp (ported
// below) through a priority-queue loop driven by glpt_tree::
// compat_bisect()/recent_leaves(), computes curve crossings on every
// final leaf's ten 2-faces (face_crossing_cache, so a face shared by
// many cells is solved once), and writes every distinct crossing point
// found as an OBJ point cloud. Deliberately NOT yet ported:
// riemann_cp2.cpp's crossing_node/cycle-extraction machinery (which
// stitches a cell's own crossing points into actual curve-segment
// connectivity) and its onion/flat projection modes -- both substantial
// pieces in their own right, left for a follow-up once this simpler
// end-to-end pipeline (mesh -> refine -> find crossings -> see them) is
// confirmed working. A point cloud is enough to visually confirm the
// extracted locus is in the right place.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <iostream>
#include "glpt_crossing.hpp"
#include "glpt_tree.hpp"

using namespace std;

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
double cell_priority(const pt3 pts[glpt::DIM+1]) {
	pt3 ghat[glpt::DIM+1];
	int n=0;
	for(int i=0;i<=glpt::DIM;++i) {
		pt3 g = vertex_gradient(pts[i]);
		double gn = hnorm(g);
		if(gn<1e-12) continue; // near a singular point: can't normalize
		for(int c=0;c<3;++c) g[c]/=gn;
		ghat[n++]=g;
	}
	if(n<2) return 1e18; // most of the cell sits right at a singular point: force refinement
	double S=0;
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			S += std::norm(hdot(ghat[i],ghat[j]));
	double kappa = S/(double(n)*double(n));
	return 1.0-kappa;
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

int main(int argc, char* argv[]) {
	int function_idx=0;
	int max_depth=8;
	double threshold=0.05;
	string obj_path="riemann_cp2_glpt.obj";

	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--function" && i+1<argc) { function_idx=atoi(argv[++i]); }
		else if(arg=="--depth" && i+1<argc) { max_depth=atoi(argv[++i]); }
		else if(arg=="--threshold" && i+1<argc) { threshold=atof(argv[++i]); }
		else if(arg=="--obj" && i+1<argc) { obj_path=argv[++i]; }
		else if(arg=="--list") { print_function_catalog_cp2(cout); return 0; }
		else {
			cerr<<"usage: "<<argv[0]<<" [--function N] [--depth N] [--threshold X] [--obj PATH] [--list]"<<endl;
			if(arg!="--help" && arg!="-h") return 1;
			return 0;
		}
	}

	vector<catalog_entry_cp2>& cat = function_catalog_cp2();
	if(function_idx<0 || function_idx>=(int)cat.size()) {
		cerr<<"bad --function index "<<function_idx<<" (--list to see the catalog)"<<endl;
		return 1;
	}
	cout<<"curve: "<<cat[function_idx].name<<" -- "<<cat[function_idx].description<<endl;
	set_curve(cat[function_idx].F);

	vector<pt3> gp = gaifullin_points();
	glpt_edge_cache id_cache;
	glpt_tree tree;
	tree.seed_all_roots();
	cout<<"seeded "<<tree.leaf_count()<<" root cells"<<endl;

	// --- Phase 1: adaptive refinement (priority-queue, gradient dispersion) --
	cout<<endl<<"--- adaptive refinement ---"<<endl;
	priority_queue<pair<double,glpt> > pq;
	for(int s=0;s<GLPT_NCELLS;++s) {
		pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
		glpt c(s);
		cell_points_and_ids(c, gp, id_cache, pts, ids);
		pq.push(make_pair(cell_priority(pts), c));
	}
	int nsubdivisions=0;
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
			pq.push(make_pair(cell_priority(pts), recent[i]));
		}
	}
	cout<<"subdivisions performed: "<<nsubdivisions<<endl;
	cout<<"final leaf count: "<<tree.leaf_count()<<endl;
	cout<<"distinct vertices minted: "<<id_cache.next_id()<<" ("<<GLPT_BASE_VERTEX_COUNT<<" base + "
		<<(id_cache.next_id()-GLPT_BASE_VERTEX_COUNT)<<" from bisection)"<<endl;

	// --- Phase 2: extraction (point cloud only -- see this file's own
	// header comment on scope) --
	cout<<endl<<"--- crossing extraction ---"<<endl;
	// One traversal: query every final leaf's ten faces (face_crossing_
	// cache solves each DISTINCT face -- by vertex-id triple -- only
	// once, however many cells touch it) and collect the crossing points
	// found. A point can still be collected more than once here (once
	// per cell whose face enumeration reaches it), so a final distance-
	// based dedup pass removes repeats before writing.
	face_crossing_cache fcache;
	vector<pt3> emitted;
	long n_cells_visited=0, n_face_queries=0;
	struct Collector {
		const vector<pt3>* gp; glpt_edge_cache* idc; face_crossing_cache* fc;
		vector<pt3>* out; long* n_cells; long* n_faces;
		void operator()(const glpt& c) const {
			++(*n_cells);
			pt3 pts[glpt::DIM+1]; int ids[glpt::DIM+1];
			cell_points_and_ids(c, *gp, *idc, pts, ids);
			for(int f=0;f<10;++f) {
				pt3 face_pts[3]; int face_ids[3];
				for(int k=0;k<3;++k) { face_pts[k]=pts[GLPT_CELL_FACES[f][k]]; face_ids[k]=ids[GLPT_CELL_FACES[f][k]]; }
				const face_result& fr = fc->get(face_pts, face_ids);
				++(*n_faces);
				for(int r=0;r<fr.nroots;++r) out->push_back(fc->root_point(fr,r));
			}
		}
	};
	Collector col; col.gp=&gp; col.idc=&id_cache; col.fc=&fcache; col.out=&emitted;
	col.n_cells=&n_cells_visited; col.n_faces=&n_face_queries;
	tree.for_each_leaf(col);
	cout<<"cells visited: "<<n_cells_visited<<", face queries: "<<n_face_queries
		<<", distinct faces solved: "<<fcache.size()<<", crossing points before dedup: "<<emitted.size()<<endl;

	// --- Phase 3: output (point cloud OBJ) --
	ofstream out(obj_path.c_str());
	if(!out) { cerr<<"couldn't open "<<obj_path<<" for writing"<<endl; return 1; }
	out<<"# riemann_cp2_glpt point cloud: curve="<<cat[function_idx].name
		<<" depth="<<max_depth<<" threshold="<<threshold<<"\n";
	long npts=0;
	for(size_t i=0;i<emitted.size();++i) {
		bool dup=false;
		for(size_t j=0;j<i && !dup;++j) {
			double d=0; for(int k=0;k<3;++k) d+=std::norm(emitted[i][k]-emitted[j][k]);
			if(d<1e-16) dup=true;
		}
		if(dup) continue;
		// Real, chart-independent embedding: project onto a fixed affine
		// chart (largest-magnitude homogeneous coordinate) just for a
		// viewable point cloud -- not the onion/flat projection modes
		// riemann_cp2.cpp offers (deferred, see this file's own header
		// comment).
		const pt3& p = emitted[i];
		int c=0; double best=-1;
		for(int k=0;k<3;++k) if(std::abs(p[k])>best) { best=std::abs(p[k]); c=k; }
		cx a,b; dehomogenize(c,p,a,b);
		out<<"v "<<a.real()<<" "<<a.imag()<<" "<<b.real()<<"\n";
		++npts;
	}
	cout<<"points written: "<<npts<<" -> "<<obj_path<<endl;
	return 0;
}
