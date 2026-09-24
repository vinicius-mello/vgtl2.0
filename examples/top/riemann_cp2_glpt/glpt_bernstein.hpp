#ifndef GLPT_BERNSTEIN_HPP
#define GLPT_BERNSTEIN_HPP

/*! \file
 * \brief Certified Bernstein-Bezier bounds/pruning (--bernstein), the
 * one piece of riemann_cp2.cpp's pipeline glpt_crossing.hpp's own
 * header comment explicitly deferred ("(2) a certified Bernstein-Bezier
 * fallback ... DEFERS (2): if every Newton seed misses, this reports 0
 * roots for now ... Can be added later the same way it was there").
 *
 * Ports riemann_cp2.cpp's own bary_poly/compose_to_chart/
 * bernstein_leaf_bounds/bernstein_bounds_recursive (certified [min,max]
 * enclosure of Re(F)/Im(F) over a 2-face's own best-conditioned chart,
 * via Bernstein-Bezier coefficients -- see that file's comment on
 * compute_bernstein_bounds for the full derivation) and
 * bernstein_locate_root (a certified fallback seed for Newton, used
 * only when every fixed/dynamic seed misses) essentially verbatim.
 *
 * ============================================================
 * WHY THIS DEPENDS ON glpt_crossing.hpp, NOT THE OTHER WAY AROUND
 * ============================================================
 * bernstein_leaf_bounds/bernstein_locate_root need pick_chart/
 * dehomogenize/build_tri_frame/bary_to_xy/newton_on_triangle/g_F --
 * all already in glpt_crossing.hpp. Rather than have glpt_crossing.hpp
 * depend back on this file (a circular header dependency), it exposes
 * two small function-pointer HOOKS (g_bernstein_maybe_zero,
 * g_bernstein_fallback, declared there, both null -- i.e. zero
 * overhead -- by default) that this file's own free functions
 * (bernstein_maybe_zero_hook/bernstein_fallback_hook, below) implement
 * and that main() wires up only when --bernstein is actually on. This
 * mirrors riemann_cp2.cpp's own forward-declaration trick (dehomogenize/
 * rehomogenize/pick_chart forward-declared right where Bernstein needs
 * them, defined for real further down) adapted to two separate files
 * instead of one .cpp's forward-declare-then-define ordering.
 *
 * ============================================================
 * bernstein_bounds_cache: THE PER-FACE ANALOGUE OF face_crossing_cache
 * ============================================================
 * riemann_cp2.cpp caches Bernstein bounds in the SAME extra_data<2>
 * struct as crossing roots, but with an independent idempotency flag
 * (bernstein_computed vs computed) -- both attached for free to nmt<4>'s
 * own per-simplex storage. glpt has no such attached storage (see
 * glpt_tree.hpp's own header comment on why), so -- exactly like
 * face_crossing_cache already does for crossing roots -- this needs its
 * own explicit cache, keyed the SAME way (make_face_key's sorted
 * global-vertex-id triple), so a 2-face shared by many 4-cells only
 * gets its bounds computed once. Kept SEPARATE from face_crossing_cache
 * itself (not folded into face_result) because bounds are needed for
 * EVERY candidate face during refinement (cell_priority_bernstein), the
 * overwhelming majority of which never go on to have their actual
 * crossing roots computed at all.
 */

#include "glpt_crossing.hpp"
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <cassert>

// --bernstein/--bernstein-level: see bernstein_bounds_recursive below.
bool g_bernstein=false;
int g_bernstein_level=0;

// A degree-d homogeneous polynomial in (l0,l1,l2), stored as its
// monomial coefficients c[a][b] for l0^a l1^b l2^(d-a-b) (a+b<=d<=4 for
// every catalog curve, so a plain 5x5 grid suffices) -- copied verbatim
// from riemann_cp2.cpp's own bary_poly.
struct bary_poly {
	int d;
	cx c[5][5];
	bary_poly() : d(0) { for(int a=0;a<5;++a) for(int b=0;b<5;++b) c[a][b]=cx(0,0); }
};

inline double fact3(int n) { double r=1; for(int i=2;i<=n;++i) r*=i; return r; }
inline double multinom3(int d,int a,int b,int cc) { return fact3(d)/(fact3(a)*fact3(b)*fact3(cc)); }

// (l0*X0+l1*X1+l2*X2)^n via the trinomial theorem.
inline bary_poly pow_linear3(cx X0,cx X1,cx X2,int n) {
	bary_poly r; r.d=n;
	for(int p=0;p<=n;++p)
		for(int q=0;q<=n-p;++q) {
			int rr=n-p-q;
			cx coeff=multinom3(n,p,q,rr)*ipow(X0,p)*ipow(X1,q)*ipow(X2,rr);
			r.c[p][q]+=coeff;
		}
	return r;
}
inline bary_poly mul_bary(const bary_poly& A, const bary_poly& B) {
	bary_poly r; r.d=A.d+B.d;
	for(int a1=0;a1<=A.d;++a1)
		for(int b1=0;a1+b1<=A.d;++b1) {
			if(A.c[a1][b1]==cx(0,0)) continue;
			for(int a2=0;a2<=B.d;++a2)
				for(int b2=0;a2+b2<=B.d;++b2)
					r.c[a1+a2][b1+b2]+=A.c[a1][b1]*B.c[a2][b2];
		}
	return r;
}

// F restricted to the face's own best-conditioned CHART (the exact
// same one compute_face_crossing picks via pick_chart), as a uniform
// degree-F.d homogeneous polynomial in barycentric (l0,l1,l2) of the
// chart's 2 free coordinates (a,b) -- NOT yet Bernstein coefficients
// (see bernstein_leaf_bounds). Copied verbatim from riemann_cp2.cpp's
// own compose_to_chart (see that file's own comment for why this must
// dehomogenize FIRST and build the barycentric patch in the CHART, not
// in the raw ambient (X,Y,Z) affine combination -- a real domain-
// mismatch bug found and fixed there).
inline bary_poly compose_to_chart(const poly_F3& F, int chart, cx a0,cx a1,cx a2, cx b0,cx b1,cx b2) {
	bary_poly acc; acc.d=F.d;
	for(size_t k=0;k<F.t.size();++k) {
		const term3& tm=F.t[k];
		int ea,eb;
		switch(chart) {
			case 0: ea=tm.e[1]; eb=tm.e[2]; break; // X=1: a=Y,b=Z
			case 1: ea=tm.e[0]; eb=tm.e[2]; break; // Y=1: a=X,b=Z
			default: ea=tm.e[0]; eb=tm.e[1]; break; // Z=1: a=X,b=Y
		}
		bary_poly pa=pow_linear3(a0,a1,a2,ea);
		bary_poly pb=pow_linear3(b0,b1,b2,eb);
		bary_poly term=mul_bary(pa,pb);
		int termdeg=ea+eb;
		if(termdeg<acc.d) {
			bary_poly elev=pow_linear3(cx(1,0),cx(1,0),cx(1,0),acc.d-termdeg); // (l0+l1+l2)^(F.d-termdeg)
			term=mul_bary(term,elev);
		}
		for(int a=0;a<=acc.d;++a)
			for(int b=0;a+b<=acc.d;++b)
				acc.c[a][b]+=tm.c*term.c[a][b];
	}
	return acc;
}

// A point in the face's chart (the 2 free affine coordinates, e.g.
// a=Y/X,b=Z/X for chart 0) -- the recursion below works entirely in
// this 2D space, matching compute_face_crossing's own domain.
struct cpt { cx a,b; };
inline cpt lerp_cpt(const cpt& p, const cpt& q) { cpt r; r.a=0.5*(p.a+q.a); r.b=0.5*(p.b+q.b); return r; }

inline void bernstein_leaf_bounds(const poly_F3& F, int chart, const cpt& P0, const cpt& P1, const cpt& P2,
		double& reLo, double& reHi, double& imLo, double& imHi) {
	bary_poly m=compose_to_chart(F,chart,P0.a,P1.a,P2.a,P0.b,P1.b,P2.b);
	reLo=1e300; reHi=-1e300; imLo=1e300; imHi=-1e300;
	for(int a=0;a<=m.d;++a) {
		for(int b=0;a+b<=m.d;++b) {
			int cc=m.d-a-b;
			cx beta=m.c[a][b]/multinom3(m.d,a,b,cc); // Bernstein coefficient
			double re=beta.real(), im=beta.imag();
			if(re<reLo) reLo=re; if(re>reHi) reHi=re;
			if(im<imLo) imLo=im; if(im>imHi) imHi=im;
		}
	}
}

inline void get_subtriangle(int s, const cpt& P0,const cpt& P1,const cpt& P2,
		const cpt& M01,const cpt& M12,const cpt& M20, cpt& A, cpt& B, cpt& C) {
	switch(s) {
		case 0: A=P0;  B=M01; C=M20; break;
		case 1: A=M01; B=P1;  C=M12; break;
		case 2: A=M20; B=M12; C=P2;  break;
		default: A=M01; B=M12; C=M20; break; // the middle, "upside-down" piece
	}
}

// Bernstein bounds for F restricted to the chart-triangle (P0,P1,P2),
// optionally tightened by `level` rounds of standard 1-to-4 triangular
// subdivision before taking the enclosure -- see riemann_cp2.cpp's own
// comment on this same function for why this certifiably only ever
// tightens the box, never loosens it.
inline void bernstein_bounds_recursive(const poly_F3& F, int chart, const cpt& P0, const cpt& P1, const cpt& P2,
		int level, double& reLo, double& reHi, double& imLo, double& imHi) {
	if(level<=0) { bernstein_leaf_bounds(F,chart,P0,P1,P2,reLo,reHi,imLo,imHi); return; }
	cpt M01=lerp_cpt(P0,P1), M12=lerp_cpt(P1,P2), M20=lerp_cpt(P2,P0);
	reLo=1e300; reHi=-1e300; imLo=1e300; imHi=-1e300;
	for(int s=0;s<4;++s) {
		cpt A,B,C; get_subtriangle(s,P0,P1,P2,M01,M12,M20,A,B,C);
		double srl,srh,sil,sih;
		bernstein_bounds_recursive(F,chart,A,B,C,level-1,srl,srh,sil,sih);
		if(srl<reLo) reLo=srl; if(srh>reHi) reHi=srh;
		if(sil<imLo) imLo=sil; if(sih>imHi) imHi=sih;
	}
}

// A certified fallback seed for Newton -- ported from riemann_cp2.cpp's
// own bernstein_locate_root (see that file's comment for the full
// rationale): only tried when every fixed/dynamic seed misses entirely.
// Descends the same 1-to-4 subdivision bernstein_bounds_recursive uses,
// following any sub-triangle whose own enclosure still straddles 0 in
// both Re and Im, carrying the corresponding barycentric corners
// (relative to the ORIGINAL, undivided face) in lockstep, and returns
// the leaf's own centroid as a certified-plausible seed.
struct bary3 { double l[3]; };
inline bary3 lerp_bary3(const bary3& a, const bary3& b) {
	bary3 r; for(int i=0;i<3;++i) r.l[i]=0.5*(a.l[i]+b.l[i]); return r;
}
inline void get_subtriangle_bary3(int s, const bary3& P0,const bary3& P1,const bary3& P2,
		const bary3& M01,const bary3& M12,const bary3& M20, bary3& A, bary3& B, bary3& C) {
	switch(s) {
		case 0: A=P0;  B=M01; C=M20; break;
		case 1: A=M01; B=P1;  C=M12; break;
		case 2: A=M20; B=M12; C=P2;  break;
		default: A=M01; B=M12; C=M20; break;
	}
}
inline bool bernstein_locate_root(const poly_F3& F, int chart,
		const cpt& P0, const cpt& P1, const cpt& P2,
		const bary3& B0, const bary3& B1, const bary3& B2,
		int level, double& out_l1, double& out_l2) {
	double reLo,reHi,imLo,imHi;
	bernstein_leaf_bounds(F,chart,P0,P1,P2,reLo,reHi,imLo,imHi);
	if(!(reLo<=0 && reHi>=0 && imLo<=0 && imHi>=0)) return false; // certified empty here
	if(level<=0) {
		out_l1=(B0.l[1]+B1.l[1]+B2.l[1])/3.0;
		out_l2=(B0.l[2]+B1.l[2]+B2.l[2])/3.0;
		return true;
	}
	cpt M01=lerp_cpt(P0,P1), M12=lerp_cpt(P1,P2), M20=lerp_cpt(P2,P0);
	bary3 N01=lerp_bary3(B0,B1), N12=lerp_bary3(B1,B2), N20=lerp_bary3(B2,B0);
	for(int s=0;s<4;++s) {
		cpt A,Bp,Cp; get_subtriangle(s,P0,P1,P2,M01,M12,M20,A,Bp,Cp);
		bary3 Ab,Bb,Cb; get_subtriangle_bary3(s,B0,B1,B2,N01,N12,N20,Ab,Bb,Cb);
		if(bernstein_locate_root(F,chart,A,Bp,Cp,Ab,Bb,Cb,level-1,out_l1,out_l2)) return true;
	}
	return false; // extremely rare: level ran out before any leaf confirmed
}
const int BFALLBACK_LEVEL=6; // narrows the seed to within 2^-6 of the face's own size

// --- Per-face Bernstein bounds cache (see this file's own header
// comment): open addressing (glpt_tree.hpp's own design, via
// glpt_hash_util.hpp), keyed by the SAME exact sorted-vertex-id triple
// face_crossing_cache uses (glpt_crossing.hpp's face_key/make_face_key
// -- a genuine two-field key, not a single lossy hash, so ids can range
// over the full 32-bit space with no collision risk; see face_key's own
// comment for why this file was reverted off a 20-bit-per-id packing).
// Values narrowed to float, rounded OUTWARD (nextafterf away from
// 0-straddling) so the stored box stays a certified superset of the
// double-precision enclosure -- same technique and reasoning as
// riemann_cp2.cpp's own compute_bernstein_bounds.
class bernstein_bounds_cache {
	public:
		struct bbounds { float reLo,reHi,imLo,imHi; };

		explicit bernstein_bounds_cache(size_t initial_buckets = 1031)
			: keys_(0), present_(0), vals_(0), nbuckets_(0), count_(0)
		{
			alloc_(glpt_next_prime(initial_buckets));
		}
		~bernstein_bounds_cache() { std::free(keys_); std::free(present_); std::free(vals_); }

		const bbounds& get(const pt3 p[3], const int id[3]) {
			face_key key = make_face_key(id[0],id[1],id[2]);
			size_t slot = find_slot_(key);
			if(slot!=size_t(-1)) return vals_[slot];

			int chart=pick_chart(p);
			cpt Q0,Q1,Q2;
			dehomogenize(chart,p[0],Q0.a,Q0.b);
			dehomogenize(chart,p[1],Q1.a,Q1.b);
			dehomogenize(chart,p[2],Q2.a,Q2.b);
			double reLo,reHi,imLo,imHi;
			bernstein_bounds_recursive(g_F,chart,Q0,Q1,Q2,g_bernstein_level,reLo,reHi,imLo,imHi);
			bbounds b;
			b.reLo=nextafterf((float)reLo,-std::numeric_limits<float>::infinity());
			b.reHi=nextafterf((float)reHi, std::numeric_limits<float>::infinity());
			b.imLo=nextafterf((float)imLo,-std::numeric_limits<float>::infinity());
			b.imHi=nextafterf((float)imHi, std::numeric_limits<float>::infinity());
			return insert_(key,b);
		}
		static bool maybe_zero(const bbounds& b) { return b.reLo<=0 && b.reHi>=0 && b.imLo<=0 && b.imHi>=0; }
		size_t size() const { return count_; }

	private:
		face_key* keys_;
		bool* present_;
		bbounds* vals_;
		size_t nbuckets_, count_;

		void alloc_(size_t n) {
			keys_ = (face_key*)std::calloc(n,sizeof(face_key));
			present_ = (bool*)std::calloc(n,sizeof(bool));
			vals_ = (bbounds*)std::calloc(n,sizeof(bbounds));
			assert(keys_!=0 && present_!=0 && vals_!=0 && "bernstein_bounds_cache: out of memory");
			nbuckets_ = n;
		}
		size_t find_slot_(const face_key& key) const {
			size_t h = glpt_hash_face_key(key) % nbuckets_;
			while(present_[h]) {
				if(keys_[h]==key) return h;
				h=(h+1)%nbuckets_;
			}
			return size_t(-1);
		}
		void grow_if_needed_() {
			if(double(count_+1) <= 0.7*double(nbuckets_)) return;
			face_key* old_k=keys_; bool* old_p=present_; bbounds* old_v=vals_;
			size_t old_n=nbuckets_;
			alloc_(glpt_next_prime(2*old_n));
			count_=0;
			for(size_t i=0;i<old_n;++i) if(old_p[i]) insert_(old_k[i], old_v[i]);
			std::free(old_k); std::free(old_p); std::free(old_v);
		}
		const bbounds& insert_(const face_key& key, const bbounds& val) {
			grow_if_needed_();
			size_t h = glpt_hash_face_key(key) % nbuckets_;
			while(present_[h]) h=(h+1)%nbuckets_;
			keys_[h]=key; present_[h]=true; vals_[h]=val;
			++count_;
			return vals_[h];
		}

		bernstein_bounds_cache(const bernstein_bounds_cache&);
		bernstein_bounds_cache& operator=(const bernstein_bounds_cache&);
};

// --- Hooks wired into glpt_crossing.hpp's g_bernstein_maybe_zero/
// g_bernstein_fallback (see that file's own comment on those two
// pointers, and this file's header comment for why they're pointers
// rather than a direct #include). g_bernstein_cache is set once, by
// main(), to the single bernstein_bounds_cache instance shared by
// cell_priority_bernstein and these hooks -- the same "curve lives in a
// global" precedent set_curve()/g_F already established.
bernstein_bounds_cache* g_bernstein_cache = 0;

inline bool bernstein_maybe_zero_hook(const pt3 p[3], const int id[3]) {
	const bernstein_bounds_cache::bbounds& b = g_bernstein_cache->get(p,id);
	return bernstein_bounds_cache::maybe_zero(b);
}
inline bool bernstein_fallback_hook(const pt3 p[3], int chart, double& out_l1, double& out_l2) {
	cpt Q0,Q1,Q2;
	dehomogenize(chart,p[0],Q0.a,Q0.b);
	dehomogenize(chart,p[1],Q1.a,Q1.b);
	dehomogenize(chart,p[2],Q2.a,Q2.b);
	bary3 B0,B1,B2;
	B0.l[0]=1; B0.l[1]=0; B0.l[2]=0;
	B1.l[0]=0; B1.l[1]=1; B1.l[2]=0;
	B2.l[0]=0; B2.l[1]=0; B2.l[2]=1;
	return bernstein_locate_root(g_F,chart,Q0,Q1,Q2,B0,B1,B2,BFALLBACK_LEVEL,out_l1,out_l2);
}

//! Call once, after --bernstein/--bernstein-level parsing and set_curve(),
//! before building any mesh: wires the hooks up and hands the cache to
//! the caller (main() keeps it alive and also uses it directly from
//! cell_priority_bernstein for box-width ranking).
inline void enable_bernstein(bernstein_bounds_cache& cache) {
	g_bernstein_cache = &cache;
	g_bernstein_maybe_zero = &bernstein_maybe_zero_hook;
	g_bernstein_fallback = &bernstein_fallback_hook;
}

#endif // GLPT_BERNSTEIN_HPP
