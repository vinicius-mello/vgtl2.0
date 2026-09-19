#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <queue>
#include <list>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/pm.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/top/cube.hpp>
#include <vgtl/top/do_nothing.hpp>
#include <vgtl/top/iso/iso.hpp>
#include <vgtl/alg/vec.hpp>
#include "functions_implicit3d.hpp"

// Plain codimension-1 isosurface triangulation in R^3: f(x,y,z)=0 for a
// single real scalar field, extracted via vgtl's own marching-simplex
// primitive (vgtl::iso::triangulate, include/vgtl/top/iso/iso.hpp) on
// an adaptively Maubach-refined tetrahedral mesh. Much simpler than
// examples/top/riemann(_PC2)'s codimension-2 curves in C^2/CP^2: one
// real equation instead of two, so a crossing point along a mesh EDGE
// is already isolated (dim(edge)+dim(surface)-dim(ambient)=1+2-3=0) --
// no per-triangle Newton solve needed, just linear interpolation of f
// along each sign-changing edge, refined by the same curvature-driven
// adaptive criterion validated in riemann_pc2.cpp (see cell_priority
// below: gradient dispersion + optional --proximity weighting, ported
// from C^3/CP^2's Hermitian inner product to R^3's plain dot()).

#define DIM 3

using namespace std;
using namespace vgtl;

typedef vgtl::nmt<DIM> T;

T t;

namespace vgtl {

	template <>
	struct extra_data<0> {
		vec<3,double> p;       // position in R^3
		double fx,fy,fz,fval;  // cached gradient and value of f at p
	};

	// Required by vgtl::iso (include/vgtl/top/iso/iso.hpp) via
	// argument-dependent lookup from within namespace vgtl::iso -- must
	// live in namespace vgtl itself, same convention as
	// examples/top/iso/iso.cpp's own signal()/signal_set().
	int g_signal_of(const extra_data<0>* d) { return d->fval>=0.0 ? +1 : -1; }
	int signal(const T& t, Vertex(T) v) { return g_signal_of(attr(t,v)); }

}

vec<3,double> point(const T& t, Vertex(T) v) { return attr(t,v)->p; }
void point_set(T& t, Vertex(T) v, const vec<3,double>& p) { attr(t,v)->p=p; }

// --- The surface and its gradient (see functions_implicit3d.hpp) ---
scalar_fun3 g_f, g_fx, g_fy, g_fz;

void compute_vertex_data(T& t, Vertex(T) v) {
	vec<3,double> p=point(t,v);
	double x=p[0],y=p[1],z=p[2];
	extra_data<0>* d=attr(t,v);
	d->fval=g_f(x,y,z);
	d->fx=g_fx(x,y,z);
	d->fy=g_fy(x,y,z);
	d->fz=g_fz(x,y,z);
}

// --- Maubach ApplyNew: plain Euclidean midpoint (flat R^3, unlike the
// Fubini-Study/spherical barycenters the projective examples needed). ---
struct refine_app : do_nothing {
	using do_nothing::apply;
	vector<Cell(T)> new_cells;

	Vertex(T) add_edge_vertex(T& t, Edge(T) e) {
		vgtl::array<Vertex(T),2> vs;
		vertices(t,e,vs);
		Vertex(T) v=add(t);
		vec<3,double> p=(point(t,vs[0])+point(t,vs[1]))*0.5;
		point_set(t,v,p);
		compute_vertex_data(t,v);
		return v;
	}
	void apply(T& t, Cell(T) s) { new_cells.push_back(s); }
};

// --- Refinement priority: same curvature (gradient dispersion) +
// optional proximity-to-surface weighting validated in
// riemann_pc2.cpp's cell_priority, ported from CP^2's Hermitian inner
// product (hdot, with its U(1) phase-gauge subtlety) to R^3's plain
// real dot() -- no gauge ambiguity here at all, gradient direction is
// unambiguous up to nothing (not even a phase), so this is if anything
// simpler than the projective case. kappa=trace(M^2)/trace(M)^2 for
// M=sum ghat_i ghat_i^T over the cell's DIM+1=4 vertices lies in
// [1/3,1] (same bound as the CP^2 case, since both have ambient
// dimension 3 -- real here, complex there); priority=1-kappa.
//
// --proximity multiplies by diam/(mind+diam), mind=min_v |f(v)|/|grad
// f(v)| (first-order distance-to-surface estimate), diam=the cell's
// own Euclidean diameter -- saturates at 1 near the surface (no blowup
// even at mind=0), decays toward 0 far from it. See
// examples/top/riemann_PC2/riemann_pc2.cpp's cell_priority for the
// full derivation/rationale (identical idea, real instead of complex).
bool g_proximity=false;

double cell_priority(const T& t, Cell(T) cv) {
	vgtl::array<Vertex(T),DIM+1> vs;
	vertices(t,cv,vs);
	vec<3,double> ghat[DIM+1], pts[DIM+1];
	int n=0;
	double mind=1e300;
	for(int i=0;i<=DIM;++i) {
		const extra_data<0>* d=attr(t,vs[i]);
		pts[i]=d->p;
		vec<3,double> g; g[0]=d->fx; g[1]=d->fy; g[2]=d->fz;
		double gn=std::sqrt(dot(g,g));
		if(gn<1e-12) continue; // near a critical point of f: can't normalize
		if(g_proximity) {
			double dv=std::fabs(d->fval)/gn;
			if(dv<mind) mind=dv;
		}
		ghat[n++]=g*(1.0/gn);
	}
	if(n<2) return 1e18; // most of the cell sits right at a critical point: force refinement
	double S=0;
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			S += dot(ghat[i],ghat[j])*dot(ghat[i],ghat[j]);
	double kappa=S/(double(n)*double(n));
	double priority=1.0-kappa;
	if(g_proximity) {
		double diam=0;
		for(int i=0;i<=DIM;++i)
			for(int j=i+1;j<=DIM;++j) {
				vec<3,double> dv=pts[i]-pts[j];
				double dd=std::sqrt(dot(dv,dv));
				if(dd>diam) diam=dd;
			}
		priority *= diam/(mind+diam);
	}
	return priority;
}

// --- OBJ export: crossing points from linear interpolation of f along
// each sign-changing edge vgtl::iso::triangulate reports, oriented via
// its own orientation() so winding stays globally consistent. ---
struct obj_writer3 {
	ofstream out;
	int nverts, nfaces;
	obj_writer3(const char* path) : out(path), nverts(0), nfaces(0) {
		out<<"# Implicit surface in R^3, extracted via vgtl::iso (marching tetrahedra)\n";
	}
	void write_triangle(const vec<3,double> pts[3]) {
		for(int i=0;i<3;++i) out<<"v "<<pts[i][0]<<" "<<pts[i][1]<<" "<<pts[i][2]<<"\n";
		out<<"f "<<(nverts+1)<<" "<<(nverts+2)<<" "<<(nverts+3)<<"\n";
		nverts+=3; ++nfaces;
	}
};

int main(int argc, char* argv[]) {

	int max_depth=8;
	double threshold=0.05;
	int function_index=0;
	for(int i=1;i<argc;++i) {
		string arg=argv[i];
		if(arg=="--depth" && i+1<argc) {
			max_depth=atoi(argv[++i]);
		} else if(arg=="--threshold" && i+1<argc) {
			threshold=atof(argv[++i]);
		} else if(arg=="--function" && i+1<argc) {
			function_index=atoi(argv[++i]);
		} else if(arg=="--proximity") {
			g_proximity=true;
		} else if(arg=="--list-functions") {
			cout<<"available functions:"<<endl;
			print_function_catalog_implicit3d(cout);
			return 0;
		} else {
			cerr<<"unrecognized argument: "<<arg<<endl;
			cerr<<"usage: "<<argv[0]<<" [--function N] [--depth N] [--threshold X] "
					<<"[--proximity] [--list-functions]"<<endl;
			return 1;
		}
	}

	vector<catalog_entry_implicit3d>& catalog=function_catalog_implicit3d();
	if(function_index<0 || function_index>=(int)catalog.size()) {
		cerr<<"--function "<<function_index<<" out of range; available:"<<endl;
		print_function_catalog_implicit3d(cerr);
		return 1;
	}
	catalog_entry_implicit3d& ce=catalog[function_index];
	g_f=ce.f; g_fx=ce.fx; g_fy=ce.fy; g_fz=ce.fz;
	cout<<"function: "<<function_index<<" ("<<ce.name<<") -- "<<ce.description<<endl;
	cout<<"max_depth="<<max_depth<<" threshold="<<threshold
			<<"  proximity="<<(g_proximity?"on":"off")<<endl;

	string obj_path_s="implicit3d_"+ce.name+".obj";
	const char* obj_path=obj_path_s.c_str();

	// --- Phase 1: seed mesh -- an axis-aligned box, Kuhn-triangulated
	// into DIM!=6 tetrahedra by vgtl's own add_cube (orientation is
	// already coherent by construction: no BFS re-derivation needed,
	// unlike the Kuhnel CP^2 seed in riemann_pc2.cpp). ---
	cout<<endl<<"--- seed mesh: axis-aligned box, Kuhn triangulation ---"<<endl;

	Vertex(T) V[8];
	for(int i=0;i<8;++i) {
		V[i]=add(t);
		vec<3,double> p;
		p[0]=(i&1)?ce.Lx:-ce.Lx;
		p[1]=(i&2)?ce.Ly:-ce.Ly;
		p[2]=(i&4)?ce.Lz:-ce.Lz;
		point_set(t,V[i],p);
		compute_vertex_data(t,V[i]);
	}
	{
		complex_builder<T> cb(t);
		add_cube(cb,V);
	}
	{
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend); ci!=cend; ++ci) level_set(t,*ci,0);
	}

	int n4=0;
	{ Cell_it(T) i,end; for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++n4; }
	cout<<"tetrahedra: "<<n4<<" (expected "<<1<<"*2*3="<<6<<")"<<endl;

	// --- Phase 2: adaptive refinement ---
	cout<<endl<<"--- adaptive refinement ---"<<endl;

	priority_queue<pair<double,Cell(T)> > pq;
	{
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend); ci!=cend; ++ci)
			if(is_current(t,*ci)) pq.push(make_pair(cell_priority(t,*ci),*ci));
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
		for(size_t i=0;i<app.new_cells.size();++i)
			pq.push(make_pair(cell_priority(t,app.new_cells[i]),app.new_cells[i]));
	}
	cout<<"subdivisions performed: "<<nsubdivisions<<endl;

	int n42=0;
	{ Cell_it(T) i,end; for(simplices(t,i,end); i!=end; ++i) if(is_current(t,*i)) ++n42; }
	cout<<"tetrahedra after refinement: "<<n42<<endl;
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

	obj_writer3 obj(obj_path);
	vector<double> resid; // |f| at each linearly-interpolated crossing point

	{
		Cell_it(T) ci,cend;
		for(simplices(t,ci,cend); ci!=cend; ++ci) {
			Cell(T) cv=*ci;
			if(!is_current(t,cv)) continue;

			vgtl::array<Vertex(T),DIM+1> vs;
			vertices(t,cv,vs);
			bool haspos=false,hasneg=false;
			for(int i=0;i<=DIM;++i) { if(signal(t,vs[i])>0) haspos=true; else hasneg=true; }
			if(!(haspos && hasneg)) continue;

			list<vgtl::array<Edge(T),DIM> > ls;
			list<int> ori;
			iso::triangulate(t,cv,ls,ori);

			list<vgtl::array<Edge(T),DIM> >::iterator li=ls.begin();
			list<int>::iterator oi=ori.begin();
			for(; li!=ls.end(); ++li,++oi) {
				vec<3,double> pts[3];
				for(int k=0;k<3;++k) {
					Edge(T) e=(*li)[k];
					vgtl::array<Vertex(T),2> ev;
					vertices(t,e,ev);
					const extra_data<0>* d0=attr(t,ev[0]);
					const extra_data<0>* d1=attr(t,ev[1]);
					double f0=d0->fval, f1=d1->fval;
					double tpar=f0/(f0-f1);
					pts[k]=d0->p*(1.0-tpar)+d1->p*tpar;
				}
				if(*oi<0) std::swap(pts[1],pts[2]);
				obj.write_triangle(pts);
				for(int k=0;k<3;++k)
					resid.push_back(std::fabs(g_f(pts[k][0],pts[k][1],pts[k][2])));
			}
		}
	}

	cout<<"wrote "<<obj_path<<": "<<obj.nverts<<" vertices, "<<obj.nfaces<<" faces"<<endl;

	if(!resid.empty()) {
		vector<double> s=resid;
		std::sort(s.begin(),s.end());
		double sum=0; for(size_t i=0;i<s.size();++i) sum+=s[i];
		cout<<endl<<"--- linear-interpolation error (|f| at each extracted crossing point, "
				<<s.size()<<" points) ---"<<endl;
		cout<<"min="<<s.front()<<"  median="<<s[s.size()/2]
				<<"  mean="<<(sum/s.size())
				<<"  p95="<<s[(size_t)(0.95*(s.size()-1))]
				<<"  max="<<s.back()<<endl;
	}

	return 0;
}
