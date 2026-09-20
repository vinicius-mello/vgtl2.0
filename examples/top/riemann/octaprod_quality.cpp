// Standalone seed-mesh quality diagnostic for examples/top/riemann/riemann.cpp's
// current seed: the simplicial (Eilenberg-Zilber shuffle) product of two
// copies of the octahedral triangulation of S^2, giving a seed for
// S^2 x S^2 = C_infty x C_infty. Coordinates/shuffle paths copied verbatim
// from riemann.cpp. Edge length = product-space chordal distance
// sqrt(sphere_dist(w1,w2)^2 + sphere_dist(z1,z2)^2), matching how the file's
// own sphere_dist() is used per-factor.
//
// Why no alternative seed is used: GAP's simpcomp SCLib was searched
// exhaustively for pure S^2xS^2 triangulations (Dim=4,Chi=4, filtering
// out CP^2#CP^2/CP^2#-CP^2 -- same chi, different manifolds -- and every
// connected-sum hit by name). Exactly 4 candidates exist: #59 (11
// vertices/68 cells) and #118/#119/#120 (12 vertices/72 cells each, VT).
// Checked all 4 for Maubach compatibility the same way find_reordering.cpp
// checked Kuhnel's CP^2 (BFS 2-coloring of the facet-adjacency dual
// graph): ALL FOUR FAIL (90-130 same-color conflicts each, out of
// ~170-190 dual edges). Fixing that would need barycentric subdivision,
// which multiplies cell count by 5!=120 -- far more than this seed's 384
// cells, for unverified quality (no coordinates available for those
// library entries to even test it). So the product construction below,
// despite not coming from a named minimal triangulation, remains the
// best available seed -- confirmed by measurement: ratio EXACTLY
// sqrt(2)=1.4142 on all 384 cells, zero variance, zero degenerate cells.
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <cstdio>
#include <algorithm>
using namespace std;

struct vec3 { double x[3]; double& operator[](int i){return x[i];} double operator[](int i) const {return x[i];} };

enum { LZERO=0, LONE=1, LMONE=2, LI=3, LMI=4, LINF=5, NLABELS=6 };

vec3 label_sphere(int l) {
	vec3 p{{0,0,0}};
	switch(l) {
		case LZERO: p[0]=0;  p[1]=0;  p[2]=-1; break;
		case LONE:  p[0]=1;  p[1]=0;  p[2]=0;  break;
		case LMONE: p[0]=-1; p[1]=0;  p[2]=0;  break;
		case LI:    p[0]=0;  p[1]=1;  p[2]=0;  break;
		case LMI:   p[0]=0;  p[1]=-1; p[2]=0;  break;
		case LINF:  p[0]=0;  p[1]=0;  p[2]=1;  break;
	}
	return p;
}

struct octa_tri { int v[3]; int sign; };
const octa_tri octahedron[8] = {
	{ {LZERO, LONE,  LI },  +1 },
	{ {LZERO, LMONE, LI },  -1 },
	{ {LZERO, LMONE, LMI }, +1 },
	{ {LZERO, LONE,  LMI }, -1 },
	{ {LINF,  LONE,  LI },  -1 },
	{ {LINF,  LMONE, LI },  +1 },
	{ {LINF,  LMONE, LMI }, -1 },
	{ {LINF,  LONE,  LMI }, +1 },
};

struct path6 { int i[5]; int j[5]; int sign; };
const path6 paths[6] = {
	{ {0,1,2,2,2}, {0,0,0,1,2}, +1 },
	{ {0,1,1,2,2}, {0,0,1,1,2}, -1 },
	{ {0,1,1,1,2}, {0,0,1,2,2}, +1 },
	{ {0,0,1,2,2}, {0,1,1,1,2}, +1 },
	{ {0,0,1,1,2}, {0,1,1,2,2}, -1 },
	{ {0,0,0,1,2}, {0,1,2,2,2}, +1 },
};

double sphere_dist(const vec3& a, const vec3& b) {
	double s=0;
	for(int i=0;i<3;++i) { double d=a[i]-b[i]; s+=d*d; }
	return sqrt(s);
}

int main() {
	// global vertex = (w_label, z_label) pair, deduplicated
	map<pair<int,int>,int> vid;
	vector<pair<int,int>> verts;
	auto get_id=[&](int wl,int zl)->int {
		auto key=make_pair(wl,zl);
		auto it=vid.find(key);
		if(it!=vid.end()) return it->second;
		int id=(int)verts.size();
		vid[key]=id;
		verts.push_back(key);
		return id;
	};

	vector<array<int,5>> cells;
	for(int oi=0; oi<8; ++oi) for(int oj=0; oj<8; ++oj) {
		for(int p=0;p<6;++p) {
			array<int,5> cell;
			for(int k=0;k<5;++k) {
				int wl=octahedron[oi].v[paths[p].i[k]];
				int zl=octahedron[oj].v[paths[p].j[k]];
				cell[k]=get_id(wl,zl);
			}
			cells.push_back(cell);
		}
	}

	printf("vertices: %zu, cells: %zu\n", verts.size(), cells.size());

	double worst_ratio=0, best_ratio=1e300, sum_ratio=0;
	double gmin=1e300, gmax=0;
	int worst_cell=-1, degenerate=0;
	for(size_t c=0;c<cells.size();++c) {
		double cmin=1e300, cmax=0;
		for(int i=0;i<5;++i) for(int j=i+1;j<5;++j) {
			int a=cells[c][i], b=cells[c][j];
			vec3 wa=label_sphere(verts[a].first), wb=label_sphere(verts[b].first);
			vec3 za=label_sphere(verts[a].second), zb=label_sphere(verts[b].second);
			double dw=sphere_dist(wa,wb), dz=sphere_dist(za,zb);
			double d=sqrt(dw*dw+dz*dz);
			cmin=min(cmin,d); cmax=max(cmax,d);
			gmin=min(gmin,d); gmax=max(gmax,d);
		}
		if(cmin<1e-12) { ++degenerate; continue; }
		double ratio=cmax/cmin;
		sum_ratio+=ratio;
		if(ratio>worst_ratio) { worst_ratio=ratio; worst_cell=(int)c; }
		best_ratio=min(best_ratio,ratio);
	}
	int n=(int)cells.size()-degenerate;
	printf("degenerate (zero-length edge) cells: %d\n", degenerate);
	printf("per-cell max/min edge ratio: worst=%.4f (cell %d), best=%.4f, mean=%.4f\n",
		worst_ratio, worst_cell, best_ratio, sum_ratio/n);
	printf("global edge length: min=%.6f, max=%.6f, mesh-wide max/min=%.4f\n",
		gmin, gmax, gmax/gmin);
	return 0;
}
