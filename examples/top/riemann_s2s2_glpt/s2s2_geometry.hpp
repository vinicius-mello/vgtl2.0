#ifndef S2S2_GEOMETRY_HPP
#define S2S2_GEOMETRY_HPP

/*! \file
 * \brief CP^1 x CP^1 geometry for the glpt pipeline: the product
 * counterpart of riemann_cp2_glpt's glpt_points/glpt_crossing/
 * glpt_extraction headers.
 *
 * A point is a pair (W,Z) of unit vectors of C^2, each defined up to its
 * own phase. Each factor carries the Fubini-Study metric
 * d(a,b)=arccos|<a,b>| (a round sphere of radius 1/2), and the product
 * carries the product metric.
 *
 *  - New vertices: the Fubini-Study midpoint in each factor.
 *  - Faces: the geodesic cone from the lowest-id vertex A over the edge
 *    BC, built factor by factor exactly as in CP^2 (glpt_crossing.hpp,
 *    "Geodesic faces"): q(t)=(1-t)B+t align(B,C), v(s,t)=(1-s)A+s
 *    align(A,q(t)), in each factor separately. An edge is then the pair
 *    of factor geodesics, each traversed by its linear parameter from the
 *    lower-id endpoint, so it is the same curve with the same parameter in
 *    every face that contains it, and the midpoint lies on it. No affine
 *    chart is involved anywhere.
 *  - The curve: F(w,z) = w^n + sum_i f_i(z) w^i from examples/top/riemann's
 *    catalog, bihomogenized with w=W0/W1, z=Z0/Z1 into P(W,Z) of bidegree
 *    (n_w,n_z) = (n, max deg f_i). The points at infinity are ordinary
 *    zeros of P. By Wirtinger's theorem the area of the curve is
 *    pi*(n_w+n_z), also for the singular closures.
 *
 * The extraction (crossing nodes, facet arcs, cycles) is
 * glpt_extraction.hpp's logic verbatim, with this point type.
 */

#include <cmath>
#include <complex>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include "glpt.hpp"
#include "glpt_vertex_ids.hpp"
#include "glpt_hash_util.hpp"
#include "../riemann/functions.hpp"

// --- points ---------------------------------------------------------------
struct c2 { cx v[2]; };
struct pp { c2 W, Z; };

inline cx dot2(const c2& a, const c2& b) { return a.v[0]*std::conj(b.v[0])+a.v[1]*std::conj(b.v[1]); }
inline double norm2(const c2& a) { return std::sqrt(std::max(0.0, dot2(a,a).real())); }
inline c2 normalize2(c2 a) { double n=norm2(a); a.v[0]/=n; a.v[1]/=n; return a; }
inline c2 mk2(cx a, cx b) { c2 r; r.v[0]=a; r.v[1]=b; return r; }
//! Rephase b so that <a,b> is real and non-negative.
inline c2 align2(const c2& a, c2 b) {
	cx ip=dot2(a,b); double m=std::abs(ip);
	if(m<1e-14) return b;
	cx ph=ip/m; b.v[0]*=ph; b.v[1]*=ph;
	return b;
}
inline double fs1(const c2& a, const c2& b) {
	double ip=std::abs(dot2(a,b))/(norm2(a)*norm2(b));
	return std::acos(std::min(1.0,ip));
}
//! Product distance sqrt(d_W^2 + d_Z^2).
inline double pdist(const pp& a, const pp& b) {
	double dw=fs1(a.W,b.W), dz=fs1(a.Z,b.Z);
	return std::sqrt(dw*dw+dz*dz);
}
inline c2 fs_mid2(const c2& a, const c2& b) {
	c2 bb=align2(a,b), s; s.v[0]=a.v[0]+bb.v[0]; s.v[1]=a.v[1]+bb.v[1];
	return normalize2(s);
}
inline pp pp_mid(const pp& a, const pp& b) { pp m; m.W=fs_mid2(a.W,b.W); m.Z=fs_mid2(a.Z,b.Z); return m; }
inline pp pp_normalize(pp a) { a.W=normalize2(a.W); a.Z=normalize2(a.Z); return a; }

// --- the seed: product of two octahedra ------------------------------------
// Octahedron on CP^1: w = 0, inf, 1, -1, i, -i (indices 0..5), balanced by
// the three antipodal pairs, colours 0,0,1,1,2,2. Its 8 triangles take one
// vertex of each pair, ordered by colour. A product cell is a staircase
// (a[i_k], b[j_k]), k=0..4, i_k+j_k=k, through two triangles a, b: vertex
// (a,b) gets colour col(a)+col(b) = k, so every cell is colour-ordered and
// the seed is balanced (glpt_seed_init asserts it). 36 vertices (id
// 6*a+b), 8*8*6 = 384 cells.
static int s2s2_cells[384][5];
inline std::vector<c2> octahedron_points() {
	double r=1.0/std::sqrt(2.0);
	std::vector<c2> P(6);
	P[0]=mk2(0,1); P[1]=mk2(1,0);
	P[2]=mk2(r,r); P[3]=mk2(-r,r);
	P[4]=mk2(cx(0,r),r); P[5]=mk2(cx(0,-r),r);
	return P;
}
inline void make_s2s2_seed() {
	int tri[8][3]; int nt=0;
	for(int a=0;a<2;++a) for(int b=2;b<4;++b) for(int c=4;c<6;++c) { tri[nt][0]=a; tri[nt][1]=b; tri[nt][2]=c; ++nt; }
	int nc=0;
	for(int s=0;s<8;++s) for(int t=0;t<8;++t)
		for(int mask=0;mask<16;++mask) {
			if(__builtin_popcount(mask)!=2) continue; // bit k set: step k advances the w-triangle
			int i=0,j=0;
			s2s2_cells[nc][0]=tri[s][0]*6+tri[t][0];
			for(int k=0;k<4;++k) {
				if(mask>>k&1) ++i; else ++j;
				s2s2_cells[nc][k+1]=tri[s][i]*6+tri[t][j];
			}
			++nc;
		}
	assert(nc==384);
	glpt_set_seed(s2s2_cells, 384, 36, false);
}
//! Seed vertex positions; Uw, Uz: unitary 2x2 moves of each factor
//! (identity = aligned seed, with (inf,inf) a seed vertex).
inline std::vector<pp> s2s2_points(const cx Uw[2][2], const cx Uz[2][2]) {
	std::vector<c2> O=octahedron_points();
	auto ap=[](const cx U[2][2], const c2& a){ return mk2(U[0][0]*a.v[0]+U[0][1]*a.v[1], U[1][0]*a.v[0]+U[1][1]*a.v[1]); };
	std::vector<pp> P(36);
	for(int a=0;a<6;++a) for(int b=0;b<6;++b) { P[6*a+b].W=normalize2(ap(Uw,O[a])); P[6*a+b].Z=normalize2(ap(Uz,O[b])); }
	return P;
}

inline void s2s2_root_points(int seed, const std::vector<pp>& gp, pp pts[glpt::DIM+1]) {
	for(int k=0;k<=glpt::DIM;++k) pts[k]=gp[glpt_seed_cells[seed][k]];
}
//! Child vertex points (glpt_points.hpp's glpt_child_vertex_points, with
//! the product midpoint).
inline void s2s2_child_points(const glpt& parent, const pp parent_points[glpt::DIM+1], int zo, pp child_points[glpt::DIM+1]) {
	const int DIM=glpt::DIM;
	glpt child=parent.child(zo);
	int pw[DIM+1][DIM+1], cw[DIM+1][DIM+1], p_shift, c_shift;
	parent.vertex_weights_exact(pw,p_shift);
	child.vertex_weights_exact(cw,c_shift);
	int scale=1<<(c_shift-p_shift), nnew=0;
	for(int k=0;k<=DIM;++k) {
		int m=-1;
		for(int j=0;j<=DIM && m<0;++j) {
			bool eq=true;
			for(int c=0;c<=DIM && eq;++c) if(cw[k][c]!=pw[j][c]*scale) eq=false;
			if(eq) m=j;
		}
		if(m>=0) child_points[k]=parent_points[m];
		else { ++nnew; child_points[k]=pp_mid(parent_points[parent.level()], parent_points[DIM]); }
	}
	assert(nnew==1); (void)nnew;
}

// --- the curve: bihomogeneous P(W,Z) ---------------------------------------
struct bterm { int i, j; cx c; }; // c * W0^i W1^(nw-i) Z0^j Z1^(nz-j)
struct bihom { int nw, nz; std::vector<bterm> t; };
bihom g_P;
inline bihom bihomogenize(const poly_F& F) {
	bihom P; P.nw=F.n; P.nz=0;
	for(int i=0;i<F.n;++i) P.nz=std::max(P.nz,(int)F.c[i].size()-1);
	bterm lead; lead.i=F.n; lead.j=0; lead.c=1; P.t.push_back(lead);
	for(int i=0;i<F.n;++i) for(size_t j=0;j<F.c[i].size();++j)
		if(std::abs(F.c[i][j])>0) { bterm b; b.i=i; b.j=(int)j; b.c=F.c[i][j]; P.t.push_back(b); }
	return P;
}
inline cx ipow(cx b, int e) { cx r(1,0); for(int k=0;k<e;++k) r*=b; return r; }
//! P and its four partial derivatives (dW0, dW1, dZ0, dZ1) at (W,Z).
inline cx peval(const c2& W, const c2& Z, cx g[4]=0) {
	cx val(0,0);
	if(g) for(int k=0;k<4;++k) g[k]=0;
	for(const bterm& b: g_P.t) {
		int a0=b.i, a1=g_P.nw-b.i, b0=b.j, b1=g_P.nz-b.j;
		cx pw0=ipow(W.v[0],a0), pw1=ipow(W.v[1],a1), pz0=ipow(Z.v[0],b0), pz1=ipow(Z.v[1],b1);
		val+=b.c*pw0*pw1*pz0*pz1;
		if(g) {
			if(a0) g[0]+=b.c*double(a0)*ipow(W.v[0],a0-1)*pw1*pz0*pz1;
			if(a1) g[1]+=b.c*double(a1)*pw0*ipow(W.v[1],a1-1)*pz0*pz1;
			if(b0) g[2]+=b.c*double(b0)*pw0*pw1*ipow(Z.v[0],b0-1)*pz1;
			if(b1) g[3]+=b.c*double(b1)*pw0*pw1*pz0*ipow(Z.v[1],b1-1);
		}
	}
	return val;
}
//! |P| at the point with unit representatives.
inline double presid(const pp& p) {
	return std::abs(peval(p.W,p.Z))/(std::pow(norm2(p.W),(double)g_P.nw)*std::pow(norm2(p.Z),(double)g_P.nz));
}
inline cx pdir(const cx g[4], const pp& e) { return g[0]*e.W.v[0]+g[1]*e.W.v[1]+g[2]*e.Z.v[0]+g[3]*e.Z.v[1]; }

// --- geodesic-cone faces ----------------------------------------------------
struct geo_face { pp A,B,C; };
//! One factor of the cone and its derivatives (glpt_crossing.hpp's
//! geo_point_derivs, in C^2).
inline void cone2(const c2& A, const c2& B, const c2& C, double s, double t, c2& v, c2& vs, c2& vt) {
	c2 Ct=align2(B,C), q, dq;
	for(int k=0;k<2;++k) { q.v[k]=(1.0-t)*B.v[k]+t*Ct.v[k]; dq.v[k]=Ct.v[k]-B.v[k]; }
	cx z=dot2(A,q), dz=dot2(A,dq);
	double az=std::abs(z);
	cx u=(az>1e-14)? z/az : cx(1,0);
	cx du=(az>1e-14)? (dz-u*std::real(std::conj(u)*dz))/az : cx(0,0);
	for(int k=0;k<2;++k) {
		cx qa=q.v[k]*u;
		v.v[k]=(1.0-s)*A.v[k]+s*qa;
		vs.v[k]=qa-A.v[k];
		vt.v[k]=s*(dq.v[k]*u+q.v[k]*du);
	}
}
inline void geo_point_derivs(const geo_face& g, double s, double t, pp& v, pp& vs, pp& vt) {
	cone2(g.A.W,g.B.W,g.C.W,s,t,v.W,vs.W,vt.W);
	cone2(g.A.Z,g.B.Z,g.C.Z,s,t,v.Z,vs.Z,vt.Z);
}
inline pp geo_point(const geo_face& g, double s, double t) { pp v,a,b; geo_point_derivs(g,s,t,v,a,b); return v; }
inline bool geo_newton(const geo_face& g, double& s, double& t) {
	pp v,vs,vt; cx gr[4];
	for(int it=0; it<40; ++it) {
		geo_point_derivs(g,s,t,v,vs,vt);
		cx G=peval(v.W,v.Z,gr);
		if(presid(v)<1e-15) break;
		cx Gs=pdir(gr,vs), Gt=pdir(gr,vt);
		double j00=Gs.real(), j01=Gt.real(), j10=Gs.imag(), j11=Gt.imag();
		double det=j00*j11-j01*j10;
		if(std::fabs(det)<1e-300) return false;
		s-=( G.real()*j11-j01*G.imag())/det;
		t-=(-G.real()*j10+j00*G.imag())/det;
		if(std::fabs(s)>10||std::fabs(t)>10) return false;
	}
	if(presid(geo_point(g,s,t))>1e-11) return false;
	const double tol=1e-6;
	return s>=-tol && s<=1+tol && t>=-tol && t<=1+tol;
}
inline void geo_sort(const int id[3], int ord[3]) {
	ord[0]=0; ord[1]=1; ord[2]=2;
	for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) if(id[ord[j]]<id[ord[i]]) std::swap(ord[i],ord[j]);
}
inline geo_face make_geo_face(const pp& A, const pp& B, const pp& C) {
	geo_face g; g.A=pp_normalize(A); g.B=pp_normalize(B); g.C=pp_normalize(C); return g;
}
inline pp align_pp(const pp& a, const pp& b) { pp r; r.W=align2(a.W,b.W); r.Z=align2(a.Z,b.Z); return r; }
//! Up to two crossings of the curve with a face; bary relative to p[] in
//! the given order (glpt_crossing.hpp's compute_face_crossing_geo).
inline int compute_face_crossing(const pp p[3], const int id[3], pp out_pts[2], double out_bary[2][3]) {
	int ord[3]; geo_sort(id,ord);
	geo_face g=make_geo_face(p[ord[0]],p[ord[1]],p[ord[2]]);
	double seeds[5][2]; int nseeds=0;
	{
		pp Bt=align_pp(g.A,g.B), Ct=align_pp(g.A,g.C);
		cx F0=peval(g.A.W,g.A.Z), F1=peval(Bt.W,Bt.Z), F2=peval(Ct.W,Ct.Z);
		cx a1=F1-F0, a2=F2-F0;
		double m00=a1.real(), m01=a2.real(), m10=a1.imag(), m11=a2.imag();
		double det=m00*m11-m01*m10;
		if(std::fabs(det)>1e-14) {
			seeds[nseeds][0]=(-F0.real()*m11+m01*F0.imag())/det;
			seeds[nseeds][1]=(-m00*F0.imag()+F0.real()*m10)/det;
			++nseeds;
		}
	}
	static const double extra[4][2]={ {1/3.,1/3.},{0.1,0.1},{0.8,0.1},{0.1,0.8} };
	for(int k=0;k<4;++k) { seeds[nseeds][0]=extra[k][0]; seeds[nseeds][1]=extra[k][1]; ++nseeds; }
	int nfound=0; double found_l[2][3];
	for(int k=0; k<nseeds && nfound<2; ++k) {
		double lB=seeds[k][0], lC=seeds[k][1];
		double s=lB+lC, t=(s>1e-12)? lC/s : 0.5;
		if(!geo_newton(g,s,t)) continue;
		s=std::min(1.0,std::max(0.0,s)); t=std::min(1.0,std::max(0.0,t));
		double l[3]={1.0-s, s*(1.0-t), s*t};
		bool dup=false;
		for(int r=0;r<nfound;++r) if(std::fabs(found_l[r][1]-l[1])<1e-7 && std::fabs(found_l[r][2]-l[2])<1e-7) dup=true;
		if(dup) continue;
		for(int q=0;q<3;++q) found_l[nfound][q]=l[q];
		out_pts[nfound]=pp_normalize(geo_point(g,s,t));
		for(int q=0;q<3;++q) out_bary[nfound][ord[q]]=l[q];
		++nfound;
	}
	return nfound;
}
//! Transversality score on the cone: |Im(conj(w1) w2)| / (|e1||e2||grad P|^2)
//! in [0,1], with e_k the face's tangent vectors at the root, made
//! horizontal in each factor at unit representatives; 0 exactly when the
//! face's tangent plane is not transverse to the curve.
inline double crossing_transversality(const pp face_pts[3], const int face_ids[3], const double bary[3]) {
	int ord[3]; geo_sort(face_ids,ord);
	geo_face g=make_geo_face(face_pts[ord[0]],face_pts[ord[1]],face_pts[ord[2]]);
	double lA=bary[ord[0]], lC=bary[ord[2]];
	double s=1.0-lA; if(s<1e-9) return 1.0;
	double t=std::min(1.0,std::max(0.0,lC/s));
	pp v,Ts,Tt; geo_point_derivs(g,s,t,v,Ts,Tt);
	auto horiz=[](c2& e, const c2& x){ double xx=dot2(x,x).real(); cx c=dot2(e,x)/xx; e.v[0]-=c*x.v[0]; e.v[1]-=c*x.v[1];
		double n=std::sqrt(xx); e.v[0]/=n; e.v[1]/=n; };
	horiz(Ts.W,v.W); horiz(Tt.W,v.W); horiz(Ts.Z,v.Z); horiz(Tt.Z,v.Z);
	pp u=pp_normalize(v); cx gr[4]; peval(u.W,u.Z,gr);
	cx w1=pdir(gr,Ts), w2=pdir(gr,Tt);
	double g2=0; for(int k=0;k<4;++k) g2+=std::norm(gr[k]);
	double n1=std::sqrt(dot2(Ts.W,Ts.W).real()+dot2(Ts.Z,Ts.Z).real());
	double n2=std::sqrt(dot2(Tt.W,Tt.W).real()+dot2(Tt.Z,Tt.Z).real());
	if(n1<1e-14||n2<1e-14) return 1.0;
	if(g2<1e-300) return 0.0;
	return std::abs(std::imag(std::conj(w1)*w2))/(n1*n2*g2);
}

// --- face cache (glpt_crossing.hpp's face_crossing_cache, with pp) ----------
struct face_key { uint64_t hi; uint32_t lo; };
inline face_key make_face_key(int a, int b, int c) {
	int v[3]={a,b,c};
	for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) if(v[j]<v[i]) std::swap(v[i],v[j]);
	face_key k; k.hi=(uint64_t(uint32_t(v[0]))<<32)|uint32_t(v[1]); k.lo=uint32_t(v[2]); return k;
}
inline bool operator==(const face_key& x, const face_key& y) { return x.hi==y.hi && x.lo==y.lo; }
inline uint64_t hash_face_key(const face_key& k) { return glpt_hash64(glpt_hash64(k.hi)^(uint64_t(k.lo)*0x9E3779B97F4A7C15ULL)); }
struct face_result { unsigned char nroots; int root_idx; };
class face_crossing_cache {
	public:
		face_crossing_cache() : keys_(0), present_(0), results_(0), nb_(0), count_(0) { alloc_(glpt_next_prime(1031)); }
		~face_crossing_cache() { std::free(keys_); std::free(present_); std::free(results_); }
		const face_result& get(const pp p[3], const int id[3]) {
			face_key key=make_face_key(id[0],id[1],id[2]);
			size_t h=hash_face_key(key)%nb_;
			while(present_[h]) { if(keys_[h]==key) return results_[h]; h=(h+1)%nb_; }
			pp out[2]; double bary[2][3];
			int nr=compute_face_crossing(p,id,out,bary);
			face_result fr; fr.nroots=(unsigned char)nr; fr.root_idx=nr? (int)pts_.size() : -1;
			for(int r=0;r<nr;++r) { pts_.push_back(out[r]); b3 b; for(int q=0;q<3;++q) b.l[q]=bary[r][q]; bary_.push_back(b); }
			return insert_(key,fr);
		}
		pp root_point(const face_result& fr, int r) const { return pts_[fr.root_idx+r]; }
		const double* root_bary(const face_result& fr, int r) const { return bary_[fr.root_idx+r].l; }
		size_t size() const { return count_; }
	private:
		struct b3 { double l[3]; };
		face_key* keys_; bool* present_; face_result* results_; size_t nb_, count_;
		std::vector<pp> pts_; std::vector<b3> bary_;
		void alloc_(size_t n) {
			keys_=(face_key*)std::calloc(n,sizeof(face_key)); present_=(bool*)std::calloc(n,sizeof(bool));
			results_=(face_result*)std::calloc(n,sizeof(face_result)); nb_=n;
			assert(keys_ && present_ && results_);
		}
		const face_result& insert_(const face_key& key, const face_result& val) {
			if(double(count_+1)>0.7*double(nb_)) {
				face_key* ok=keys_; bool* op=present_; face_result* orr=results_; size_t on=nb_;
				alloc_(glpt_next_prime(2*on)); count_=0;
				for(size_t i=0;i<on;++i) if(op[i]) insert_(ok[i],orr[i]);
				std::free(ok); std::free(op); std::free(orr);
			}
			size_t h=hash_face_key(key)%nb_;
			while(present_[h]) h=(h+1)%nb_;
			keys_[h]=key; present_[h]=true; results_[h]=val; ++count_;
			return results_[h];
		}
		face_crossing_cache(const face_crossing_cache&);
		face_crossing_cache& operator=(const face_crossing_cache&);
};

// --- per-cell extraction (glpt_extraction.hpp, with pp) ---------------------
static const int CELL_FACES[10][3] = {
	{0,1,2},{0,1,3},{0,1,4},{0,2,3},{0,2,4},{0,3,4},{1,2,3},{1,2,4},{1,3,4},{2,3,4} };
static const int FACET_SUBFACES[5][4] = { {6,7,8,9},{3,4,5,9},{1,2,5,8},{0,2,4,7},{0,1,3,6} };
const double CROSSING_SNAP_TOL=1e-6, CROSSING_MERGE_TOL=1e-9;
struct crossing_node { int dim; int a,b; double param; face_key fkey; int sub; pp p; };
inline bool same_crossing_node(const crossing_node& x, const crossing_node& y) {
	if(x.dim!=y.dim) return false;
	if(x.dim==0) return x.a==y.a;
	if(x.dim==1) return x.a==y.a && x.b==y.b && std::fabs(x.param-y.param)<CROSSING_MERGE_TOL;
	return x.fkey==y.fkey && x.sub==y.sub;
}
inline void compute_crossing_node(const int fid[3], const double bary[3], const face_key& fk, int sub, const pp& p, crossing_node& nd) {
	int zeros=0, zi[3];
	for(int b=0;b<3;++b) if(std::fabs(bary[b])<CROSSING_SNAP_TOL) zi[zeros++]=b;
	nd.p=p;
	if(zeros>=2) { int keep=(zeros==3)?0:(3-zi[0]-zi[1]); nd.dim=0; nd.a=fid[keep]; }
	else if(zeros==1) {
		int lo=(zi[0]==0)?1:0, hi=(zi[0]==2)?1:2;
		int va=fid[lo], vb=fid[hi]; double t=bary[hi]/(bary[lo]+bary[hi]);
		if(va>vb) { std::swap(va,vb); t=1.0-t; }
		nd.dim=1; nd.a=va; nd.b=vb; nd.param=t;
	} else { nd.dim=2; nd.fkey=fk; nd.sub=sub; }
}
struct cell_extraction_result {
	std::vector<crossing_node> nodes;
	std::vector<std::vector<int> > cycles;
	int ntouching_tets, nbad_tets;
	bool any_edge, decompose_ok;
	double min_transversality;
	int fail_reason; // 0 ok/no arc, 1 irregular facet, 2 node degree, 3 bad cycle
};
inline void extract_cell(const pp pts[5], const int ids[5], face_crossing_cache& fc, cell_extraction_result& out) {
	const int DIM=4;
	out.nodes.clear(); out.cycles.clear();
	out.ntouching_tets=0; out.nbad_tets=0; out.min_transversality=1.0; out.fail_reason=0;
	std::vector<std::vector<int> > adj;
	bool tet_unhandled=false;
	for(int j=0;j<=DIM;++j) {
		int tn[4], ntn=0;
		for(int s=0;s<4;++s) {
			int f=FACET_SUBFACES[j][s];
			pp fp[3]; int fid[3];
			for(int k=0;k<3;++k) { fp[k]=pts[CELL_FACES[f][k]]; fid[k]=ids[CELL_FACES[f][k]]; }
			const face_result& fr=fc.get(fp,fid);
			for(int r=0;r<fr.nroots;++r) {
				double tv=crossing_transversality(fp,fid,fc.root_bary(fr,r));
				if(tv<out.min_transversality) out.min_transversality=tv;
				crossing_node nd;
				compute_crossing_node(fid,fc.root_bary(fr,r),make_face_key(fid[0],fid[1],fid[2]),r,fc.root_point(fr,r),nd);
				int idx=-1;
				for(size_t q=0;q<out.nodes.size();++q) if(same_crossing_node(out.nodes[q],nd)) { idx=(int)q; break; }
				if(idx<0) { idx=(int)out.nodes.size(); out.nodes.push_back(nd); adj.push_back(std::vector<int>()); }
				bool dup=false; for(int q=0;q<ntn;++q) if(tn[q]==idx) dup=true;
				if(!dup && ntn<4) tn[ntn++]=idx;
			}
		}
		if(ntn==2) { adj[tn[0]].push_back(tn[1]); adj[tn[1]].push_back(tn[0]); }
		else if(ntn==4) {
			static const int pr[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
			int best=0; double bc=1e300;
			for(int c=0;c<3;++c) {
				double cost=pdist(out.nodes[tn[pr[c][0]]].p,out.nodes[tn[pr[c][1]]].p)+pdist(out.nodes[tn[pr[c][2]]].p,out.nodes[tn[pr[c][3]]].p);
				if(cost<bc) { bc=cost; best=c; }
			}
			for(int e=0;e<2;++e) { int a=tn[pr[best][2*e]], b=tn[pr[best][2*e+1]]; adj[a].push_back(b); adj[b].push_back(a); }
		}
		else if(ntn==1) ++out.ntouching_tets;
		else if(ntn!=0) { ++out.nbad_tets; tet_unhandled=true; }
	}
	bool deg_ok=true; out.any_edge=false;
	for(size_t q=0;q<out.nodes.size();++q) { size_t dg=adj[q].size(); if(!dg) continue; out.any_edge=true; if(dg!=2) deg_ok=false; }
	out.decompose_ok=out.any_edge && deg_ok && !tet_unhandled;
	if(!out.decompose_ok) { if(out.any_edge) out.fail_reason=tet_unhandled?1:2; return; }
	std::vector<bool> seen(out.nodes.size(),false);
	for(size_t q=0;q<out.nodes.size();++q) {
		if(seen[q]||adj[q].empty()) continue;
		std::vector<int> cyc; int start=(int)q, prev=start, cur=adj[start][0];
		seen[q]=true; cyc.push_back(start);
		while(cur!=start) {
			if(seen[cur]) { out.decompose_ok=false; out.fail_reason=3; return; }
			seen[cur]=true; cyc.push_back(cur);
			int nxt=(adj[cur][0]==prev)?adj[cur][1]:adj[cur][0];
			prev=cur; cur=nxt;
			if(cyc.size()>(size_t)(DIM+1)) { out.decompose_ok=false; out.fail_reason=3; return; }
		}
		if(cyc.size()<3) { out.decompose_ok=false; out.fail_reason=3; return; }
		out.cycles.push_back(cyc);
	}
}

#endif // S2S2_GEOMETRY_HPP
