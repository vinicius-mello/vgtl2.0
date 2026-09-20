// Standalone seed-mesh quality diagnostic: compares Kuhnel's 9-vertex
// CP^2 (fails Maubach compatibility, kept only as a numeric reference)
// against Gaifullin's 15-vertex CP^2 (the seed actually used by
// riemann_pc2.cpp), using max/min Fubini-Study edge-length ratio per
// cell as the quality metric.
//
// Coordinates and cell lists copied verbatim from riemann_pc2.cpp
// (Gaifullin, current HEAD) and from commit 4b0e61b (Kuhnel, since the
// project since removed that seed from the live file).
//
// Why only these two: GAP's simpcomp SCLib was searched exhaustively,
// not just by name -- SCLibSearchByAttribute(SCLib,"Dim=4 and Chi=3")
// over the full 648-entry library (chi=3 at dim 4 is necessary for any
// manifold with CP^2's Betti numbers) returns EXACTLY these two entries
// (#16 "CP^2 (VT)", #397 "Gaifullin CP^2") -- there is no third CP^2
// triangulation anywhere in the library to compare against. Measured
// result: Kuhnel is exactly edge-regular (ratio 1.0 on all 36 cells) but
// unusable directly (dual graph not bipartite, see find_reordering.cpp);
// Gaifullin measures exactly 1.4636 on all 108 cells, zero variance, no
// degenerate cells.
#include <complex>
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <algorithm>
using namespace std;

typedef complex<double> cx;
typedef array<cx,3> pt3;

pt3 mkpt(cx a, cx b, cx c) { return pt3{a,b,c}; }

cx hdot(const pt3& a, const pt3& b) {
	cx r(0,0);
	for(int i=0;i<3;++i) r += a[i]*conj(b[i]);
	return r;
}
double hnorm(const pt3& a) { return sqrt(max(0.0, hdot(a,a).real())); }
pt3 normalize3(pt3 a) {
	double n=hnorm(a);
	for(int i=0;i<3;++i) a[i]/=n;
	return a;
}
double fs_dist(const pt3& a, const pt3& b) {
	double ip=abs(hdot(a,b));
	if(ip>1.0) ip=1.0;
	return acos(ip);
}

typedef array<int,5> cell5;

vector<pt3> hesse_points() {
	double s3=sqrt(3.0);
	cx one(1,0), zero(0,0);
	cx omega(-0.5, s3/2.0), omega2(-0.5,-s3/2.0);
	vector<pt3> P(9);
	P[0]=mkpt(zero, one, -one);
	P[1]=mkpt(zero, one, -omega);
	P[2]=mkpt(zero, one, -omega2);
	P[3]=mkpt(one, zero, -one);
	P[4]=mkpt(one, zero, -omega);
	P[5]=mkpt(one, zero, -omega2);
	P[6]=mkpt(one, -one, zero);
	P[7]=mkpt(one, -omega, zero);
	P[8]=mkpt(one, -omega2, zero);
	for(int i=0;i<9;++i) P[i]=normalize3(P[i]);
	return P;
}

int apply_S(int v) { static const int m[9]={3,4,5,6,7,8,0,1,2}; return m[v]; }
cell5 apply_S(const cell5& c) {
	cell5 r; for(int i=0;i<5;++i) r[i]=apply_S(c[i]);
	sort(r.begin(),r.end());
	return r;
}
const int base_raw[12][5] = {
	{0,4,1,7,8}, {0,1,2,7,8}, {0,2,5,7,8}, {3,4,1,7,8}, {3,1,2,7,8}, {3,2,5,7,8},
	{0,3,1,4,5}, {0,3,2,4,5}, {0,3,1,4,8}, {0,3,2,5,7},
	{0,3,6,1,5}, {0,3,6,5,7},
};
vector<cell5> kuhnel_cells() {
	vector<cell5> cells;
	for(int b=0;b<12;++b) {
		cell5 c; for(int i=0;i<5;++i) c[i]=base_raw[b][i];
		sort(c.begin(),c.end());
		cell5 cur=c;
		for(int step=0; step<3; ++step) {
			bool dup=false;
			for(size_t k=0;k<cells.size();++k) if(cells[k]==cur) { dup=true; break; }
			if(!dup) cells.push_back(cur);
			cur=apply_S(cur);
		}
	}
	return cells;
}

vector<pt3> gaifullin_points() {
	double s3=sqrt(3.0);
	cx one(1,0), zero(0,0);
	cx omega(-0.5, s3/2.0), omega2(-0.5,-s3/2.0);
	vector<pt3> P(15);
	P[0]=mkpt(zero,zero,one);
	P[1]=mkpt(zero,one,zero);
	P[2]=mkpt(one,zero,zero);
	P[3]=mkpt(-one,omega,omega2);
	P[4]=mkpt(-one,omega2,omega);
	P[5]=mkpt(-one,one,one);
	P[6]=mkpt(one,-omega,omega2);
	P[7]=mkpt(one,-omega2,omega);
	P[8]=mkpt(one,-one,one);
	P[9]=mkpt(one,omega,-omega2);
	P[10]=mkpt(one,omega2,-omega);
	P[11]=mkpt(one,one,-one);
	P[12]=mkpt(one,omega,omega2);
	P[13]=mkpt(one,omega2,omega);
	P[14]=mkpt(one,one,one);
	for(int i=0;i<15;++i) P[i]=normalize3(P[i]);
	return P;
}
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

struct Stats {
	double worst_ratio=0, best_ratio=1e300;
	double global_min=1e300, global_max=0;
	double sum_ratio=0;
	int n=0;
	int worst_cell=-1;
};

void analyze(const char* label, const vector<pt3>& pts, const vector<cell5>& cells) {
	Stats s;
	for(size_t c=0;c<cells.size();++c) {
		double cmin=1e300, cmax=0;
		for(int i=0;i<5;++i) for(int j=i+1;j<5;++j) {
			double d=fs_dist(pts[cells[c][i]], pts[cells[c][j]]);
			cmin=min(cmin,d); cmax=max(cmax,d);
			s.global_min=min(s.global_min,d);
			s.global_max=max(s.global_max,d);
		}
		double ratio = cmin>1e-12 ? cmax/cmin : 1e300;
		s.sum_ratio += ratio;
		++s.n;
		if(ratio>s.worst_ratio) { s.worst_ratio=ratio; s.worst_cell=(int)c; }
		s.best_ratio=min(s.best_ratio,ratio);
	}
	printf("=== %s ===\n", label);
	printf("  cells: %d\n", s.n);
	printf("  per-cell max/min edge ratio: worst=%.4f (cell %d), best=%.4f, mean=%.4f\n",
		s.worst_ratio, s.worst_cell, s.best_ratio, s.sum_ratio/s.n);
	printf("  global edge length (FS dist, radians): min=%.6f, max=%.6f, mesh-wide max/min=%.4f\n",
		s.global_min, s.global_max, s.global_max/s.global_min);
	printf("\n");
}

int main() {
	analyze("Kuhnel 9-vertex / 36 cells (fails Maubach compatibility -- reference only)",
		hesse_points(), kuhnel_cells());
	analyze("Gaifullin 15-vertex / 108 cells (current riemann_pc2.cpp seed)",
		gaifullin_points(), gaifullin_cells());
	return 0;
}
