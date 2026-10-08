// riemann_s2s2_glpt.cpp -- the riemann_cp2_glpt pipeline on CP^1 x CP^1.
//
// Same method as examples/top/riemann_cp2_glpt (continuation through a
// pointerless glpt mesh, per-cell extraction from crossings on geodesic-
// cone 2-faces, repair rounds, tangency threshold, a-posteriori topology),
// with only the seed and the geometry changed (s2s2_geometry.hpp): the
// seed is the product of two octahedra (balanced, 384 cells), new vertices
// are Fubini-Study midpoints in each factor, faces are geodesic cones built
// factor by factor, and the curve F(w,z) of examples/top/riemann's catalog
// is evaluated bihomogeneously -- no charts, no special case at infinity.
//
// Written to compare the two compactifications on the same code: the
// closure of w^2=z^3-z has a cusp at (inf,inf) in CP^1 x CP^1 but is
// smooth in CP^2. --cusp-diag bins cells and area by their distance to
// that corner.
//
// Build (from this directory):
//   g++ -std=c++17 -O2 -I../../../include -I$HOME/code/lpt riemann_s2s2_glpt.cpp
//       $HOME/code/lpt/lpt.cpp -o riemann_s2s2_glpt

#define GLPT_SEED_BITS 9 // 384 seed cells
#define GLPT_ORTH_BITS 9 // depth up to 39
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include "s2s2_geometry.hpp"
#include "glpt_tree.hpp"

using namespace std;

inline double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
inline long peak_rss_kb() {
	std::ifstream f("/proc/self/status"); std::string line;
	while(std::getline(f,line)) if(line.rfind("VmHWM:",0)==0) return std::atol(line.c_str()+6);
	return -1;
}

// Random unitary 2x2 (Gram-Schmidt of a complex Gaussian matrix).
void random_unitary2(cx U[2][2], std::mt19937& rng) {
	std::normal_distribution<double> nd(0.0,1.0);
	c2 a=normalize2(mk2(cx(nd(rng),nd(rng)),cx(nd(rng),nd(rng))));
	c2 b=mk2(cx(nd(rng),nd(rng)),cx(nd(rng),nd(rng)));
	cx p=dot2(b,a); b.v[0]-=p*a.v[0]; b.v[1]-=p*a.v[1]; b=normalize2(b);
	for(int i=0;i<2;++i) { U[i][0]=a.v[i]; U[i][1]=b.v[i]; }
}

void path_from_root(const glpt& c, int path[], int& n) {
	n=0; glpt cur=c;
	while(!cur.is_root()) { path[n++]=cur.is_child0()?0:1; cur=cur.parent(); }
	for(int i=0;i<n/2;++i) std::swap(path[i],path[n-1-i]);
}
void cell_points_and_ids(const glpt& c, const vector<pp>& gp, glpt_edge_cache& cache, pp pts[5], int ids[5]) {
	int path[64], n; path_from_root(c,path,n);
	s2s2_root_points(c.seed_index(),gp,pts);
	glpt_root_vertex_ids(c.seed_index(),ids);
	glpt cur(c.seed_index());
	for(int i=0;i<n;++i) {
		pp cp[5]; int ci[5];
		s2s2_child_points(cur,pts,path[i],cp);
		glpt_child_vertex_ids(cur,ids,path[i],cache,ci);
		cur=cur.child(path[i]);
		for(int k=0;k<5;++k) { pts[k]=cp[k]; ids[k]=ci[k]; }
	}
}
double cell_diam(const pp pts[5]) {
	double d=0; for(int a=0;a<5;++a) for(int b=a+1;b<5;++b) d=std::max(d,pdist(pts[a],pts[b]));
	return d;
}

// --- output-vertex identity (riemann_cp2_glpt's node_key) ------------------
struct node_key {
	int dim, a, b; long long q;
	bool operator<(const node_key& o) const {
		if(dim!=o.dim) return dim<o.dim;
		if(a!=o.a) return a<o.a;
		if(b!=o.b) return b<o.b;
		return q<o.q;
	}
};
node_key key_of(const crossing_node& nd) {
	node_key k; k.dim=nd.dim; k.a=0; k.b=0; k.q=0;
	if(nd.dim==0) k.a=nd.a;
	else if(nd.dim==1) { k.a=nd.a; k.b=nd.b; k.q=(long long)(nd.param/CROSSING_MERGE_TOL+0.5); }
	else { k.a=int(nd.fkey.hi>>32); k.b=int(nd.fkey.hi&0xFFFFFFFFu); k.q=(long long)nd.fkey.lo*2+nd.sub; }
	return k;
}

// --- topology and area (riemann_cp2_glpt's genus_check/close_surface) ------
static int dsu_find(std::vector<int>& p, int x) { while(p[x]!=x) { p[x]=p[p[x]]; x=p[x]; } return x; }
//! Product-metric area: each polygon fanned into triangles measured in the
//! horizontal tangent space (per factor) at the fan's first vertex.
double surface_area(const vector<pp>& P, const vector<vector<int> >& faces) {
	double area=0;
	auto tang=[](const c2& a, const c2& b){ c2 e=align2(a,b); e.v[0]-=a.v[0]; e.v[1]-=a.v[1];
		cx c=dot2(e,a); e.v[0]-=c*a.v[0]; e.v[1]-=c*a.v[1]; return e; };
	for(const vector<int>& f: faces) {
		const pp& a=P[f[0]];
		for(size_t i=1;i+1<f.size();++i) {
			c2 e1w=tang(a.W,P[f[i]].W), e1z=tang(a.Z,P[f[i]].Z), e2w=tang(a.W,P[f[i+1]].W), e2z=tang(a.Z,P[f[i+1]].Z);
			double n1=dot2(e1w,e1w).real()+dot2(e1z,e1z).real(), n2=dot2(e2w,e2w).real()+dot2(e2z,e2z).real();
			double r=dot2(e1w,e2w).real()+dot2(e1z,e2z).real();
			area+=0.5*std::sqrt(std::max(0.0,n1*n2-r*r));
		}
	}
	return area;
}
void report_area(const vector<pp>& P, const vector<vector<int> >& faces) {
	double A=surface_area(P,faces), W=M_PI*(g_P.nw+g_P.nz);
	vector<double> res; for(const pp& p: P) res.push_back(presid(p));
	sort(res.begin(),res.end());
	cout<<"area = "<<A<<"  (Wirtinger pi*(n_w+n_z) = "<<W<<", ratio "<<A/W<<")"<<endl;
	if(!res.empty()) cout<<"|P| at "<<res.size()<<" vertices: median="<<res[res.size()/2]<<" max="<<res.back()<<endl;
}
pp corner_point();
static const double CUSP_EDGES_G[]={0.01,0.02,0.05,0.1,0.2,0.5,1e300};
static const int NCUSP_G=7;
inline int cusp_bin_g(double d) { int b=0; while(b<NCUSP_G-1 && d>=CUSP_EDGES_G[b]) ++b; return b; }
void genus_check(const vector<vector<int> >& faces_in, int nverts, const vector<pp>* P=0) {
	vector<vector<pair<int,int> > > corners(nverts);
	for(size_t f=0;f<faces_in.size();++f) for(size_t i=0;i<faces_in[f].size();++i) corners[faces_in[f][i]].push_back(make_pair((int)f,(int)i));
	vector<vector<int> > faces(faces_in);
	int nid=0, pinched=0;
	for(int v=0;v<nverts;++v) {
		const vector<pair<int,int> >& cs=corners[v];
		if(cs.empty()) continue;
		map<int,int> li; vector<int> par;
		auto idx=[&](int x){ auto it=li.find(x); if(it!=li.end()) return it->second; int k=(int)par.size(); li[x]=k; par.push_back(k); return k; };
		vector<int> ca(cs.size());
		for(size_t j=0;j<cs.size();++j) {
			const vector<int>& F=faces_in[cs[j].first]; int n=(int)F.size(), i=cs[j].second;
			int a=idx(F[(i+n-1)%n]), b=idx(F[(i+1)%n]);
			int ra=dsu_find(par,a), rb=dsu_find(par,b); if(ra!=rb) par[ra]=rb;
			ca[j]=a;
		}
		map<int,int> comp;
		for(size_t j=0;j<cs.size();++j) { int r=dsu_find(par,ca[j]); if(!comp.count(r)) { int k=(int)comp.size(); comp[r]=k; } }
		if(comp.size()>1) ++pinched;
		for(size_t j=0;j<cs.size();++j) faces[cs[j].first][cs[j].second]=nid+comp[dsu_find(par,ca[j])];
		nid+=(int)comp.size();
	}
	map<pair<int,int>,int> edges;
	for(auto& f: faces) { int n=(int)f.size(); for(int i=0;i<n;++i) { int a=f[i], b=f[(i+1)%n]; edges[make_pair(min(a,b),max(a,b))]++; } }
	long V=nid, E=(long)edges.size(), F=(long)faces.size(), chi=V-E+F;
	vector<int> pb(nid), ps(nid); for(int i=0;i<nid;++i) pb[i]=ps[i]=i;
	vector<char> onb(nid,0); long nonman=0;
	for(auto& e: edges) {
		if(e.second==1) { onb[e.first.first]=onb[e.first.second]=1; int r1=dsu_find(pb,e.first.first), r2=dsu_find(pb,e.first.second); if(r1!=r2) pb[r1]=r2; }
		if(e.second>2) ++nonman;
		int s1=dsu_find(ps,e.first.first), s2=dsu_find(ps,e.first.second); if(s1!=s2) ps[s1]=s2;
	}
	set<int> loops, comps;
	for(int i=0;i<nid;++i) { if(onb[i]) loops.insert(dsu_find(pb,i)); comps.insert(dsu_find(ps,i)); }
	long b=(long)loops.size(), c=(long)comps.size();
	cout<<endl<<"--- genus check (complex in CP^1 x CP^1) ---"<<endl;
	cout<<"pinched vertices split: "<<pinched<<", non-manifold edges: "<<nonman<<", components c="<<c<<endl;
	cout<<"V="<<V<<" E="<<E<<" F="<<F<<" chi="<<chi<<" boundary loops b="<<b<<endl;
	cout<<"genus estimate (2c-b-chi)/2 = "<<(2.0*c-b-chi)/2.0
		<<"   (smooth curve of bidegree ("<<g_P.nw<<","<<g_P.nz<<"): "<<(g_P.nw-1)*(g_P.nz-1)<<")"<<endl;
	// Orientability: polygons come with arbitrary cycle direction; orient
	// them by BFS across manifold edges (an edge shared by two faces must be
	// traversed in opposite directions) and count the edges where this fails.
	{
		map<pair<int,int>,vector<pair<int,int> > > ef; // undirected edge -> (face, +1 if traversed a->b with a<b)
		for(size_t f=0;f<faces.size();++f) { int n=(int)faces[f].size();
			for(int i=0;i<n;++i) { int a=faces[f][i], bb=faces[f][(i+1)%n]; ef[make_pair(min(a,bb),max(a,bb))].push_back(make_pair((int)f, a<bb?1:-1)); } }
		vector<vector<pair<int,int> > > nbr(faces.size()); // (other face, required relative sign)
		for(auto& e: ef) if(e.second.size()==2) {
			int f1=e.second[0].first, f2=e.second[1].first, r=-e.second[0].second*e.second[1].second;
			nbr[f1].push_back(make_pair(f2,r)); nbr[f2].push_back(make_pair(f1,r));
		}
		vector<int> sgn(faces.size(),0); long conflicts=0; double dmin=1e300;
		for(size_t f0=0;f0<faces.size();++f0) {
			if(sgn[f0]) continue;
			sgn[f0]=1; vector<int> q(1,(int)f0);
			while(!q.empty()) {
				int f=q.back(); q.pop_back();
				for(auto& nb: nbr[f]) {
					int want=sgn[f]*nb.second;
					if(!sgn[nb.first]) { sgn[nb.first]=want; q.push_back(nb.first); }
					else if(sgn[nb.first]!=want) {
						++conflicts;
						if(P) { pp cn=corner_point(); for(int x: faces_in[f]) dmin=min(dmin,pdist((*P)[x],cn)); }
					}
				}
			}
		}
		cout<<"orientation conflicts: "<<conflicts/2;
		if(conflicts && P) cout<<" (nearest to (inf,inf): "<<dmin<<")";
		cout<<(conflicts?"  -> NOT orientable":"  -> orientable")<<endl;
	}
	// Complex orientation: on a holomorphic curve the Kahler form omega is the
	// area form of the complex orientation (Wirtinger), so a polygon close to
	// the curve has omega/area close to +-1, and its sign orients it. Orient
	// every polygon by that sign; a shared edge traversed the same way by both
	// faces is then a place where the extraction disagrees with the complex
	// orientation. Binned by distance to (inf,inf).
	if(P) {
		const vector<pp>& Q=*P;
		auto tang=[](const c2& a, const c2& b){ c2 e=align2(a,b); e.v[0]-=a.v[0]; e.v[1]-=a.v[1];
			cx c=dot2(e,a); e.v[0]-=c*a.v[0]; e.v[1]-=c*a.v[1]; return e; };
		vector<int> osgn(faces_in.size()); vector<double> ratio(faces_in.size());
		for(size_t f=0;f<faces_in.size();++f) {
			const vector<int>& F=faces_in[f]; const pp& a=Q[F[0]];
			double om=0, ar=0;
			for(size_t i=1;i+1<F.size();++i) {
				c2 e1w=tang(a.W,Q[F[i]].W), e1z=tang(a.Z,Q[F[i]].Z), e2w=tang(a.W,Q[F[i+1]].W), e2z=tang(a.Z,Q[F[i+1]].Z);
				cx h=dot2(e2w,e1w)+dot2(e2z,e1z); // <e2,e1>; Im = omega(e1,e2) for the standard structure
				double n1=dot2(e1w,e1w).real()+dot2(e1z,e1z).real(), n2=dot2(e2w,e2w).real()+dot2(e2z,e2z).real();
				om+=0.5*h.imag(); ar+=0.5*std::sqrt(std::max(0.0,n1*n2-h.real()*h.real()));
			}
			osgn[f]=om>=0?1:-1; ratio[f]=ar>0? om/ar : 0;
		}
		map<pair<int,int>,vector<pair<int,int> > > ef;
		for(size_t f=0;f<faces_in.size();++f) { int n=(int)faces_in[f].size();
			for(int i=0;i<n;++i) { int a=faces_in[f][i], bb=faces_in[f][(i+1)%n];
				ef[make_pair(min(a,bb),max(a,bb))].push_back(make_pair((int)f, (a<bb?1:-1)*osgn[f])); } }
		long nbin[NCUSP_G]={0}, badbin[NCUSP_G]={0}, lowr[NCUSP_G]={0}, fbin[NCUSP_G]={0};
		pp cn=corner_point();
		for(auto& e: ef) if(e.second.size()==2) {
			double d=min(pdist(Q[e.first.first],cn),pdist(Q[e.first.second],cn));
			int b=cusp_bin_g(d); ++nbin[b];
			if(e.second[0].second==e.second[1].second) ++badbin[b];
		}
		for(size_t f=0;f<faces_in.size();++f) {
			double d=1e300; for(int x: faces_in[f]) d=min(d,pdist(Q[x],cn));
			int b=cusp_bin_g(d); ++fbin[b]; if(std::fabs(ratio[f])<0.9) ++lowr[b];
		}
		cout<<"complex orientation, by distance to (inf,inf): interior edges, edges inconsistent with it; faces, faces with |omega/area|<0.9"<<endl;
		double lo=0;
		for(int b=0;b<NCUSP_G;++b) { if(nbin[b]||fbin[b]) cout<<"  ["<<lo<<", "<<CUSP_EDGES_G[b]<<"): edges "<<nbin[b]<<" inconsistent "<<badbin[b]
			<<"; faces "<<fbin[b]<<" low ratio "<<lowr[b]<<endl; lo=CUSP_EDGES_G[b]; }
	}
}
//! Newton onto the curve from p: minimum-norm step in the affine chart of
//! each factor's largest coordinate.
pp newton_onto_curve(pp p, bool& ok) {
	int kw=std::abs(p.W.v[1])>=std::abs(p.W.v[0]) ? 1 : 0, kz=std::abs(p.Z.v[1])>=std::abs(p.Z.v[0]) ? 1 : 0;
	c2 W=p.W, Z=p.Z;
	cx sw=W.v[kw], sz=Z.v[kz];
	W.v[0]/=sw; W.v[1]/=sw; Z.v[0]/=sz; Z.v[1]/=sz;
	ok=false;
	for(int it=0;it<30;++it) {
		cx g[4]; cx f=peval(W,Z,g);
		cx ga=g[1-kw], gb=g[2+(1-kz)];
		double n2=std::norm(ga)+std::norm(gb);
		if(n2<1e-300) break;
		W.v[1-kw]-=f*std::conj(ga)/n2; Z.v[1-kz]-=f*std::conj(gb)/n2;
		pp q; q.W=W; q.Z=Z;
		if(presid(q)<1e-13) { ok=true; break; }
	}
	pp r; r.W=W; r.Z=Z; return pp_normalize(r);
}
void close_surface(vector<pp>& P, vector<vector<int> >& faces) {
	int nv=(int)P.size();
	vector<vector<pair<int,int> > > corners(nv);
	for(size_t f=0;f<faces.size();++f) for(size_t i=0;i<faces[f].size();++i) corners[faces[f][i]].push_back(make_pair((int)f,(int)i));
	vector<vector<int> > nf(faces); vector<pp> NP; int pinched=0;
	for(int v=0;v<nv;++v) {
		const vector<pair<int,int> >& cs=corners[v]; if(cs.empty()) continue;
		map<int,int> li; vector<int> par;
		auto idx=[&](int x){ auto it=li.find(x); if(it!=li.end()) return it->second; int k=(int)par.size(); li[x]=k; par.push_back(k); return k; };
		vector<int> ca(cs.size());
		for(size_t j=0;j<cs.size();++j) {
			const vector<int>& F=faces[cs[j].first]; int n=(int)F.size(), i=cs[j].second;
			int a=idx(F[(i+n-1)%n]), b=idx(F[(i+1)%n]); int ra=dsu_find(par,a), rb=dsu_find(par,b); if(ra!=rb) par[ra]=rb; ca[j]=a;
		}
		map<int,int> comp;
		for(size_t j=0;j<cs.size();++j) { int r=dsu_find(par,ca[j]); if(!comp.count(r)) { int k=(int)comp.size(); comp[r]=k; } }
		if(comp.size()>1) ++pinched;
		int base=(int)NP.size(); for(size_t k=0;k<comp.size();++k) NP.push_back(P[v]);
		for(size_t j=0;j<cs.size();++j) nf[cs[j].first][cs[j].second]=base+comp[dsu_find(par,ca[j])];
	}
	int n2=(int)NP.size(); vector<int> ps(n2); for(int i=0;i<n2;++i) ps[i]=i;
	for(auto& f: nf) for(size_t i=1;i<f.size();++i) { int a=dsu_find(ps,f[0]), b=dsu_find(ps,f[i]); if(a!=b) ps[a]=b; }
	map<int,int> csize; for(auto& f: nf) csize[dsu_find(ps,f[0])]+=1;
	int mainr=-1, best=-1; for(auto& c: csize) if(c.second>best) { best=c.second; mainr=c.first; }
	vector<vector<int> > kept; int dropped=0;
	for(auto& f: nf) { if(dsu_find(ps,f[0])==mainr) kept.push_back(f); else ++dropped; }
	map<pair<int,int>,int> ec;
	for(auto& f: kept) for(size_t i=0;i<f.size();++i) { int a=f[i], b=f[(i+1)%f.size()]; ec[make_pair(min(a,b),max(a,b))]++; }
	vector<pair<int,int> > bnd;
	for(auto& f: kept) for(size_t i=0;i<f.size();++i) { int a=f[i], b=f[(i+1)%f.size()]; if(ec[make_pair(min(a,b),max(a,b))]==1) bnd.push_back(make_pair(b,a)); }
	vector<int> pb(n2); for(int i=0;i<n2;++i) pb[i]=i;
	for(auto& e: bnd) { int a=dsu_find(pb,e.first), b=dsu_find(pb,e.second); if(a!=b) pb[a]=b; }
	map<int,vector<int> > loopv; for(auto& e: bnd) loopv[dsu_find(pb,e.first)].push_back(e.first);
	map<int,int> centre; int onc=0;
	for(auto& L: loopv) {
		const pp& ref=NP[L.second[0]]; c2 sw=mk2(0,0), sz=mk2(0,0);
		for(int x: L.second) { c2 a=align2(ref.W,NP[x].W), b=align2(ref.Z,NP[x].Z); for(int k=0;k<2;++k) { sw.v[k]+=a.v[k]; sz.v[k]+=b.v[k]; } }
		pp c; c.W=normalize2(sw); c.Z=normalize2(sz);
		bool ok; pp cc=newton_onto_curve(c,ok);
		double r=0; for(int x: L.second) r=max(r,pdist(c,NP[x]));
		if(ok && pdist(c,cc)<=2*r+1e-12) { c=cc; ++onc; }
		centre[L.first]=(int)NP.size(); NP.push_back(c);
	}
	for(auto& e: bnd) kept.push_back(vector<int>{e.first,e.second,centre[dsu_find(pb,e.first)]});
	vector<int> remap(NP.size(),-1); vector<pp> out;
	for(auto& f: kept) for(int& x: f) { if(remap[x]<0) { remap[x]=(int)out.size(); out.push_back(NP[x]); } x=remap[x]; }
	P.swap(out); faces.swap(kept);
	cout<<endl<<"--- close ---"<<endl<<"pinched vertices split: "<<pinched<<", island faces dropped: "<<dropped
		<<", holes capped: "<<loopv.size()<<" (cone apex moved onto the curve: "<<onc<<")"<<endl;
}

// --kahler-filter X: omega/area of an extracted polygon, omega the Kahler
// form (product of the factors' Fubini-Study forms), measured on the fan at
// its first vertex. By Wirtinger it is +-1 exactly on complex lines, so a
// polygon close to the curve has |ratio| near 1; a cell with a polygon
// below X is treated as bad.
double g_kahler_filter=0;
double polygon_kahler_ratio(const vector<crossing_node>& nodes, const vector<int>& cyc) {
	auto tang=[](const c2& a, const c2& b){ c2 e=align2(a,b); e.v[0]-=a.v[0]; e.v[1]-=a.v[1];
		cx c=dot2(e,a); e.v[0]-=c*a.v[0]; e.v[1]-=c*a.v[1]; return e; };
	const pp& a=nodes[cyc[0]].p; double om=0, ar=0;
	for(size_t i=1;i+1<cyc.size();++i) {
		const pp& b=nodes[cyc[i]].p; const pp& c=nodes[cyc[i+1]].p;
		c2 e1w=tang(a.W,b.W), e1z=tang(a.Z,b.Z), e2w=tang(a.W,c.W), e2z=tang(a.Z,c.Z);
		cx h=dot2(e2w,e1w)+dot2(e2z,e1z);
		double n1=dot2(e1w,e1w).real()+dot2(e1z,e1z).real(), n2=dot2(e2w,e2w).real()+dot2(e2z,e2z).real();
		om+=0.5*h.imag(); ar+=0.5*std::sqrt(std::max(0.0,n1*n2-h.real()*h.real()));
	}
	return ar>0 ? om/ar : 1.0;
}

// --- --cusp-diag: cells and area by distance to the corner (inf,inf) -------
static const double CUSP_EDGES[]={0.01,0.02,0.05,0.1,0.2,0.5,1e300};
static const int NCUSP=7;
pp corner_point() { pp c; c.W=mk2(1,0); c.Z=mk2(1,0); return c; }
int cusp_bin(double d) { int b=0; while(b<NCUSP-1 && d>=CUSP_EDGES[b]) ++b; return b; }

int main(int argc, char* argv[]) {
	int function_idx=0, max_depth=10, repair_rounds=0, repair_extra_depth=4, tangency_extra_depth=4;
	double tangency_threshold=0;
	bool generic=false, genus=false, close=false, cusp_diag=false, neighbor_test=false, timing=false;
	unsigned generic_seed=12345;
	string obj_path;
	for(int i=1;i<argc;++i) {
		string a=argv[i];
		if(a=="--function" && i+1<argc) function_idx=atoi(argv[++i]);
		else if(a=="--depth" && i+1<argc) max_depth=atoi(argv[++i]);
		else if(a=="--generic") generic=true;
		else if(a=="--generic-seed" && i+1<argc) { generic=true; generic_seed=(unsigned)atoi(argv[++i]); }
		else if(a=="--repair-rounds" && i+1<argc) repair_rounds=atoi(argv[++i]);
		else if(a=="--repair-extra-depth" && i+1<argc) repair_extra_depth=atoi(argv[++i]);
		else if(a=="--tangency-threshold" && i+1<argc) tangency_threshold=atof(argv[++i]);
		else if(a=="--tangency-extra-depth" && i+1<argc) tangency_extra_depth=atoi(argv[++i]);
		else if(a=="--genus-check") genus=true;
		else if(a=="--close") { close=true; genus=true; }
		else if(a=="--cusp-diag") cusp_diag=true;
		else if(a=="--kahler-filter" && i+1<argc) g_kahler_filter=atof(argv[++i]);
		else if(a=="--neighbor-test") neighbor_test=true;
		else if(a=="--timing") timing=true;
		else if(a=="--obj" && i+1<argc) obj_path=argv[++i];
		else if(a=="--list") { print_function_catalog(cout); return 0; }
		else {
			cerr<<"usage: "<<argv[0]<<" [--function N] [--depth N] [--generic] [--generic-seed N] [--repair-rounds N] [--repair-extra-depth N] "
				<<"[--tangency-threshold X] [--tangency-extra-depth N] [--genus-check] [--close] [--cusp-diag] [--kahler-filter X] [--neighbor-test] [--timing] [--obj PATH] [--list]"<<endl;
			return (a=="--help"||a=="-h") ? 0 : 1;
		}
	}
	vector<catalog_entry>& cat=function_catalog();
	if(function_idx<0 || function_idx>=(int)cat.size()) { cerr<<"bad --function"<<endl; return 1; }
	g_P=bihomogenize(cat[function_idx].F);
	cout<<"curve: "<<cat[function_idx].name<<" -- "<<cat[function_idx].description<<endl;
	cout<<"bidegree (n_w,n_z) = ("<<g_P.nw<<","<<g_P.nz<<")  generic="<<(generic?"on":"off");
	if(generic) cout<<" (seed="<<generic_seed<<")";
	cout<<endl;

	make_s2s2_seed();
	cx Uw[2][2]={{1,0},{0,1}}, Uz[2][2]={{1,0},{0,1}};
	if(generic) { std::mt19937 rng(generic_seed); random_unitary2(Uw,rng); random_unitary2(Uz,rng); }
	vector<pp> gp=s2s2_points(Uw,Uz);
	glpt_edge_cache id_cache;
	glpt_tree tree;
	tree.seed_all_roots();
	cout<<"seeded "<<tree.leaf_count()<<" root cells (octahedron x octahedron, 36 vertices)"<<endl;

	// --- Phase 1: continuation (riemann_cp2_glpt's default Phase 1)
	double t0=now_s();
	face_crossing_cache fcache;
	long nsub=0, ntested=0, ntang=0;
	int seed_root=-1;
	for(int s=0;s<glpt_seed_count() && seed_root<0;++s) {
		pp pts[5]; int ids[5]; cell_points_and_ids(glpt(s),gp,id_cache,pts,ids);
		cell_extraction_result res; extract_cell(pts,ids,fcache,res); ++ntested;
		if(res.any_edge) seed_root=s;
	}
	if(seed_root<0) { cout<<"curve doesn't touch any root cell"<<endl; return 1; }
	{
		set<glpt> visited; vector<glpt> frontier(1,glpt(seed_root));
		while(!frontier.empty()) {
			glpt c=frontier.back(); frontier.pop_back();
			if(!tree.exists(c) || visited.count(c)) continue;
			visited.insert(c);
			pp pts[5]; int ids[5]; cell_points_and_ids(c,gp,id_cache,pts,ids);
			cell_extraction_result res; extract_cell(pts,ids,fcache,res); ++ntested;
			if(!res.any_edge) continue;
			bool tang=tangency_threshold>0 && res.min_transversality<tangency_threshold && c.simplex_level()<max_depth+tangency_extra_depth;
			if(tang && c.simplex_level()>=max_depth) ++ntang;
			if(c.simplex_level()<max_depth || tang) {
				tree.clear_recent(); tree.compat_bisect(c); ++nsub;
				for(const glpt& r: tree.recent_leaves()) frontier.push_back(r);
			} else {
				for(int i=0;i<=4;++i) { glpt nb; if(tree.neighbor_leaf(c,i,nb)==glpt_tree::FOUND && !visited.count(nb)) frontier.push_back(nb); }
			}
		}
	}
	cout<<"seed root "<<seed_root<<"; subdivisions "<<nsub<<", cells tested "<<ntested;
	if(tangency_threshold>0) cout<<", forced deeper by --tangency-threshold "<<ntang;
	cout<<endl<<"leaf count after continuation: "<<tree.leaf_count()<<endl;

	if(neighbor_test) {
		struct T { glpt_tree* tr; const vector<pp>* gp; glpt_edge_cache* ic; long *nq,*nbad,*npos;
			void operator()(const glpt& c) const {
				pp P[5]; int I[5]; cell_points_and_ids(c,*gp,*ic,P,I);
				for(int i=0;i<5;++i) {
					glpt nb; if(tr->neighbor_leaf(c,i,nb)!=glpt_tree::FOUND) continue;
					if(nb.simplex_level()!=c.simplex_level()) continue;
					++*nq; pp Q[5]; int J[5]; cell_points_and_ids(nb,*gp,*ic,Q,J);
					for(int k=0;k<5;++k) { if(k==i) continue; int m=-1; for(int l=0;l<5;++l) if(J[l]==I[k]) m=l;
						if(m<0) { ++*nbad; break; } if(pdist(P[k],Q[m])>1e-9) { ++*npos; break; } }
				} } };
		long nq=0,nbad=0,npos=0; T t; t.tr=&tree; t.gp=&gp; t.ic=&id_cache; t.nq=&nq; t.nbad=&nbad; t.npos=&npos;
		tree.for_each_leaf(t);
		cout<<"neighbor test: "<<nq<<" same-level neighbour queries, "<<nbad<<" with a facet vertex missing, "<<npos<<" with a position mismatch"<<endl;
	}

	// --- Phase 2: extraction (+ repair rounds)
	double t1=now_s();
	struct Ext {
		const vector<pp>* gp; glpt_edge_cache* idc; face_crossing_cache* fc;
		long ok, bad, kahler, reason[4], okb[NCUSP], badb[NCUSP];
		vector<glpt> bad_list; vector<pp> vpts; map<node_key,int> vidx; vector<vector<int> > faces; vector<int> face_bin;
		void clear() { ok=bad=kahler=0; for(int i=0;i<4;++i) reason[i]=0; for(int b=0;b<NCUSP;++b) okb[b]=badb[b]=0;
			bad_list.clear(); vpts.clear(); vidx.clear(); faces.clear(); face_bin.clear(); }
		int emit(const crossing_node& nd) {
			node_key k=key_of(nd); auto it=vidx.find(k); if(it!=vidx.end()) return it->second;
			int i=(int)vpts.size(); vpts.push_back(nd.p); vidx[k]=i; return i;
		}
		void operator()(const glpt& c) {
			pp pts[5]; int ids[5]; cell_points_and_ids(c,*gp,*idc,pts,ids);
			cell_extraction_result res; extract_cell(pts,ids,*fc,res);
			if(!res.any_edge) return;
			double dc=1e300; pp cn=corner_point(); for(int k=0;k<5;++k) dc=min(dc,pdist(pts[k],cn));
			int b=cusp_bin(dc);
			// Orient every polygon by the complex orientation: reverse the cycle
			// when omega/area < 0 (see polygon_kahler_ratio), so faces sharing an
			// edge traverse it in opposite directions wherever the extraction is
			// consistent with the curve.
			if(res.decompose_ok)
				for(auto& cyc: res.cycles) {
					double r=polygon_kahler_ratio(res.nodes,cyc);
					if(g_kahler_filter>0 && std::fabs(r)<g_kahler_filter) { res.decompose_ok=false; res.fail_reason=0; ++kahler; break; }
					if(r<0) std::reverse(cyc.begin()+1,cyc.end()); // keeps cyc[0], the fan apex
				}
			if(!res.decompose_ok) { ++bad; ++reason[res.fail_reason]; ++badb[b]; bad_list.push_back(c); return; }
			++ok; ++okb[b];
			for(auto& cyc: res.cycles) { vector<int> f; for(int x: cyc) f.push_back(emit(res.nodes[x])); faces.push_back(f); face_bin.push_back(b); }
		}
	};
	Ext ext; ext.gp=&gp; ext.idc=&id_cache; ext.fc=&fcache; ext.clear();
	{ struct W { Ext* e; void operator()(const glpt& c) const { (*e)(c); } }; W w; w.e=&ext; tree.for_each_leaf(w);
	for(int round=0; round<repair_rounds && !ext.bad_list.empty(); ++round) {
		int nrep=0;
		for(const glpt& c: ext.bad_list) { if(!tree.exists(c) || c.simplex_level()>=max_depth+repair_extra_depth) continue; tree.compat_bisect(c); ++nrep; }
		cout<<"repair round "<<round+1<<": "<<ext.bad_list.size()<<" bad cells, "<<nrep<<" bisected further"<<endl;
		if(!nrep) break;
		ext.clear(); tree.for_each_leaf(w);
	} }
	cout<<"final leaf count: "<<tree.leaf_count()<<", distinct faces solved: "<<fcache.size()<<endl;
	long touched=ext.ok+ext.bad;
	cout<<"extracted cells: ok="<<ext.ok<<" bad="<<ext.bad<<"  (bad rate "<<(touched?100.0*ext.bad/touched:0.0)<<"%)"<<endl;
	cout<<"bad cells: irregular facet="<<ext.reason[1]<<" node degree="<<ext.reason[2]<<" bad cycle="<<ext.reason[3]
		<<" Kahler filter="<<ext.kahler<<endl;
	if(cusp_diag) {
		double ab[NCUSP]={0};
		for(size_t f=0;f<ext.faces.size();++f) { vector<vector<int> > one(1,ext.faces[f]); ab[ext.face_bin[f]]+=surface_area(ext.vpts,one); }
		cout<<"cells meeting C, by min distance of their vertices to (inf,inf) (product FS metric):"<<endl;
		double lo=0;
		for(int b=0;b<NCUSP;++b) {
			if(ext.okb[b]||ext.badb[b]) cout<<"  ["<<lo<<", "<<(CUSP_EDGES[b]>1e299?string("inf"):to_string(CUSP_EDGES[b]))<<"): ok="<<ext.okb[b]<<" bad="<<ext.badb[b]
				<<"  area="<<ab[b]<<endl;
			lo=CUSP_EDGES[b];
		}
	}
	report_area(ext.vpts,ext.faces);
	// one machine-readable line (figures/s2s2_sweep.sh)
	long near_ok=0, near_bad=0; for(int b=0;b<NCUSP;++b) if(CUSP_EDGES[b]<=0.1) { near_ok+=ext.okb[b]; near_bad+=ext.badb[b]; }
	double area_raw=surface_area(ext.vpts,ext.faces)/(M_PI*(g_P.nw+g_P.nz));
	if(genus) {
		genus_check(ext.faces,(int)ext.vpts.size(),&ext.vpts);
		if(close) { close_surface(ext.vpts,ext.faces); genus_check(ext.faces,(int)ext.vpts.size(),&ext.vpts); report_area(ext.vpts,ext.faces); }
	}
	cout<<"SUMMARY leaves="<<tree.leaf_count()<<" ok="<<ext.ok<<" bad="<<ext.bad<<" kahler="<<ext.kahler
		<<" near_ok="<<near_ok<<" near_bad="<<near_bad<<" area_raw="<<area_raw
		<<" area_out="<<surface_area(ext.vpts,ext.faces)/(M_PI*(g_P.nw+g_P.nz))<<endl;
	if(timing) cout<<"timing: refinement "<<t1-t0<<" s, extraction "<<now_s()-t1<<" s, peak RSS "<<peak_rss_kb()/1024.0<<" MiB"<<endl;
	if(!obj_path.empty()) { // (Re w, Im w, Re z), w=W0/W1, z=Z0/Z1
		ofstream out(obj_path.c_str());
		for(const pp& p: ext.vpts) { cx w=p.W.v[0]/p.W.v[1], z=p.Z.v[0]/p.Z.v[1]; out<<"v "<<w.real()<<" "<<w.imag()<<" "<<z.real()<<"\n"; }
		for(auto& f: ext.faces) { out<<"f"; for(int x: f) out<<" "<<x+1; out<<"\n"; }
		cout<<"wrote "<<obj_path<<endl;
	}
	return 0;
}
