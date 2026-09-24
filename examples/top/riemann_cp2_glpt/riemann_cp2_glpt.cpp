// riemann_cp2_glpt.cpp -- glpt_tree-based rebuild of riemann_cp2.cpp's
// pipeline (see glpt_points.hpp/glpt_crossing.hpp/glpt_extraction.hpp's
// own header comments for the three pieces this builds on, and the
// conversation that led to this file for the full rationale: an
// O(1)-per-current-cell mesh, ~/code/lpt/glpt_tree.hpp, replacing
// riemann_cp2.cpp's append-only nmt<4>-backed one).
//
// Seeds all 108 Gaifullin cells, refines by the SAME gradient-dispersion
// cell_priority as riemann_cp2.cpp (ported below) through a priority-
// queue loop driven by glpt_tree::compat_bisect()/recent_leaves(), then
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
// -- see random_unitary()'s own comment below). Deliberately NOT
// ported: the certified Bernstein-Bezier fallback (only needed when
// every Newton seed misses a root -- matches riemann_cp2.cpp's own
// behavior before that fallback existed, not a new gap) and
// riemann_cp2.cpp's --onion mode (flat+alpha are enough to see the
// extracted surface; onion is deferred, not required).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <fstream>
#include <iostream>
#include <random>
#include "glpt_extraction.hpp"
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
int pick_chart_polygon(const std::vector<crossing_node>& nodes, const std::vector<int>& cyc) {
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
	else { k.a=int(nd.fkey&0xFFFFFFFFu); k.b=int((nd.fkey>>32)&0xFFFFFFFFu); k.qparam=nd.sub; }
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
	// projection: same default/flags as riemann_cp2.cpp's own --flat
	// (implicit default)/--alpha-projection/--cutoff (onion deliberately
	// not ported here -- see this file's own header comment).
	bool alpha_mode=false;
	double cutoff=1e300;
	bool generic=false;
	unsigned generic_seed=12345;

	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--function" && i+1<argc) { function_idx=atoi(argv[++i]); }
		else if(arg=="--depth" && i+1<argc) { max_depth=atoi(argv[++i]); }
		else if(arg=="--threshold" && i+1<argc) { threshold=atof(argv[++i]); }
		else if(arg=="--obj" && i+1<argc) { obj_path=argv[++i]; }
		else if(arg=="--flat") { alpha_mode=false; }
		else if(arg=="--alpha-projection") { alpha_mode=true; }
		else if(arg=="--cutoff" && i+1<argc) { cutoff=atof(argv[++i]); }
		else if(arg=="--generic") { generic=true; }
		else if(arg=="--generic-seed" && i+1<argc) { generic=true; generic_seed=(unsigned)atoi(argv[++i]); }
		else if(arg=="--list") { print_function_catalog_cp2(cout); return 0; }
		else {
			cerr<<"usage: "<<argv[0]<<" [--function N] [--depth N] [--threshold X] [--obj PATH] "
				<<"[--flat] [--alpha-projection] [--cutoff X] [--generic] [--generic-seed N] [--list]"<<endl;
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
	cout<<"generic="<<(generic?"on":"off");
	if(generic) cout<<" (seed="<<generic_seed<<")";
	cout<<endl;
	set_curve(cat[function_idx].F);

	vector<pt3> gp = gaifullin_points();
	if(generic) {
		cx M[3][3]; random_unitary(M,generic_seed);
		for(size_t i=0;i<gp.size();++i) gp[i]=apply_unitary(M,gp[i]);
	}
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
	cout<<endl<<"--- surface extraction ---"<<endl;
	face_crossing_cache fcache;
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

	struct Extractor {
		const vector<pt3>* gp; glpt_edge_cache* idc; face_crossing_cache* fc;
		long* n_cells; int *ok, *bad, *ntouch, *nbad_t, *npoly, *nclip;
		bool alpha_mode; double cutoff;
		vector<pt3>* vpts; map<node_key,int>* vidx; vector<vector<int> >* faces;
		vector<vec3d>* fverts; vector<vector<int> >* ffaces;

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
			if(!res.any_edge) return; // curve doesn't cross this cell at all -- not a failure, nothing to count
			*ntouch += res.ntouching_tets;
			*nbad_t += res.nbad_tets;
			if(!res.decompose_ok) { ++(*bad); return; }
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
	Extractor ext;
	ext.gp=&gp; ext.idc=&id_cache; ext.fc=&fcache; ext.n_cells=&n_cells_visited;
	ext.ok=&ok_cells; ext.bad=&bad_cells; ext.ntouch=&ntouching_tets; ext.nbad_t=&nbad_tets; ext.npoly=npoly_out;
	ext.nclip=&nclipped; ext.alpha_mode=alpha_mode; ext.cutoff=cutoff;
	ext.vpts=&vert_pts; ext.vidx=&vert_index; ext.faces=&faces_out;
	ext.fverts=&flat_verts; ext.ffaces=&flat_faces;
	tree.for_each_leaf(ext);

	cout<<"cells visited: "<<n_cells_visited<<", distinct faces solved: "<<fcache.size()<<endl;
	cout<<"extracted cells: ok="<<ok_cells<<" bad="<<bad_cells<<endl;
	cout<<"polygon sizes:";
	for(int sz=3;sz<8;++sz) if(npoly_out[sz]) cout<<" "<<sz<<"-gon="<<npoly_out[sz];
	cout<<endl;
	cout<<"touching tetrahedra: "<<ntouching_tets<<"  unhandled node count: "<<nbad_tets<<endl;
	if(!alpha_mode && cutoff<1e299) cout<<"polygons dropped by --cutoff: "<<nclipped<<endl;

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
