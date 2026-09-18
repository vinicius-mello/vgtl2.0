// Tests the thesis claim: one barycentric subdivision of a simplicial
// pseudomanifold is enough to get a global vertex ordering per cell
// satisfying Maubach's face-compatibility condition
//   d_i(sigma1) = d_j(sigma2)  =>  i = j
// WITHOUT needing an edge-coloring search (unlike find_reordering.cpp,
// which needed one because it worked on the ORIGINAL 36-cell complex).
//
// --- Why this should work, in outline (verified computationally below,
// not just argued) ---------------------------------------------------
// Every top-dimensional simplex of Sd(K), the barycentric subdivision of
// a pure d-dimensional simplicial complex K, corresponds to a maximal
// flag of faces of some single top-simplex sigma of K:
//   F_0 (subset) F_1 (subset) F_2 (subset) F_3 (subset) F_4 = sigma
// i.e. to a permutation of sigma's own 5 vertices (F_j = the first j+1
// vertices in that order). Its 5 vertices in Sd(K) are the barycenters
// b(F_0),...,b(F_4) -- and they have a CANONICAL global order: by the
// dimension (rank) of F_j, 0..4. This rank is intrinsic to each
// barycenter (dimension of the face it's the center of), not a per-cell
// choice -- so there's no coloring/search needed, unlike the original
// 36-cell complex where "position" had to be solved for.
//
// The real question -- what this program actually checks -- is whether
// two flags sharing a facet (dropping one rank-level from the chain)
// always drop the SAME rank on both sides. This is true for the two
// easy cases (dropping the top or bottom of the chain), but for the
// "middle" cases needs the actual combinatorics of which two original
// top-simplices share which tetrahedron in K -- so it's worth an actual
// computational check rather than taking it on faith, which is what
// this program does, on the real (verified) Kuhnel CP^2_9 facet list.
//
// As a bonus, also checks whether this reordering, together with an
// orientation sign derived from a FIXED GLOBAL vertex order (0..8) --
// there's no more "ascending order of the cell's own vertices" once
// facets straddle different original cells -- gives a globally coherent
// orientation. Unlike the compatibility question, there's no clean
// abstract argument either way for this part (worked through by hand
// first: the "interior" flag-graph is bipartite by permutation parity,
// but the edges crossing between original top-simplices don't have an
// a-priori-guaranteed parity behavior), so this part is genuinely only
// answered by running it.

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <cassert>
#include <cstdint>

using namespace std;

typedef array<int,5> Cell;

int apply_S(int v) { static const int m[9]={3,4,5,6,7,8,0,1,2}; return m[v]; }
Cell apply_S(const Cell& c) { Cell r; for(int i=0;i<5;++i) r[i]=apply_S(c[i]); sort(r.begin(),r.end()); return r; }

const int base_raw[12][5] = {
	{0,4,1,7,8}, {0,1,2,7,8}, {0,2,5,7,8}, {3,4,1,7,8}, {3,1,2,7,8}, {3,2,5,7,8},
	{0,3,1,4,5}, {0,3,2,4,5}, {0,3,1,4,8}, {0,3,2,5,7},
	{0,3,6,1,5}, {0,3,6,5,7},
};

int main() {
	// --- Regenerate and verify the 36 Kuhnel cells (same as find_reordering.cpp) ---
	vector<Cell> cells;
	for(int b=0;b<12;++b) {
		Cell c; for(int i=0;i<5;++i) c[i]=base_raw[b][i]; sort(c.begin(),c.end());
		Cell cur=c;
		for(int step=0; step<3; ++step) {
			bool dup=false;
			for(size_t k=0;k<cells.size();++k) if(cells[k]==cur) { dup=true; break; }
			if(!dup) cells.push_back(cur);
			cur=apply_S(cur);
		}
		assert(cur==c);
	}
	assert(cells.size()==36);
	int N=(int)cells.size();
	cout<<"base 36-cell complex regenerated and verified (see find_reordering.cpp "
			<<"for the full checks -- facet-sharing, vertex-transitivity, f-vector, "
			<<"Euler characteristic -- skipped here for brevity, unchanged)"<<endl;

	// --- Generate all 36*120 = 4320 maximal flags -------------------------
	// flag[f] = the permutation tuple (5 vertex labels, in rank order 0..4)
	vector<Cell> flag;
	flag.reserve(N*120);
	for(int c=0;c<N;++c) {
		array<int,5> perm=cells[c]; // start sorted ascending, then permute
		sort(perm.begin(),perm.end());
		do {
			flag.push_back(perm);
		} while(next_permutation(perm.begin(),perm.end()));
	}
	cout<<"generated flags: "<<flag.size()<<" (expected 36*120=4320)"<<endl;
	assert(flag.size()==(size_t)N*120);

	// For each flag, F[j] = bitmask of the first j+1 elements (j=0..4).
	auto prefix_masks=[&](const Cell& t)->array<uint32_t,5> {
		array<uint32_t,5> F; uint32_t m=0;
		for(int j=0;j<5;++j) { m|=(1u<<t[j]); F[j]=m; }
		return F;
	};

	int M=(int)flag.size();
	vector<array<uint32_t,5> > Fm(M);
	for(int i=0;i<M;++i) Fm[i]=prefix_masks(flag[i]);

	// --- Build facet ownership: facet key = the 4 surviving F_j's (as a
	// sorted-by-size sequence of bitmasks -- NOT "delete this raw tuple
	// entry", which is wrong for interior ranks: two different fillings
	// of a dropped middle rank reorder the remaining raw entries but keep
	// the same SET at every surviving rank -- see the file comment above). --
	typedef array<uint32_t,4> FacetKey;
	map<FacetKey, vector<pair<int,int> > > owner; // facet -> [(flag idx, dropped rank), ...]
	for(int i=0;i<M;++i) {
		for(int k=0;k<5;++k) {
			FacetKey key; int p=0;
			for(int j=0;j<5;++j) if(j!=k) key[p++]=Fm[i][j];
			owner[key].push_back(make_pair(i,k));
		}
	}
	cout<<"distinct facets: "<<owner.size()<<" (expected 4320*5/2=10800 if closed)"<<endl;
	{
		map<int,int> hist;
		for(map<FacetKey,vector<pair<int,int> > >::iterator it=owner.begin(); it!=owner.end(); ++it)
			hist[(int)it->second.size()]++;
		cout<<"facet-owner-count histogram (expect {2: 10800}): ";
		for(map<int,int>::iterator it=hist.begin();it!=hist.end();++it) cout<<"{"<<it->first<<": "<<it->second<<"} ";
		cout<<endl;
		assert(hist.size()==1 && hist.count(2) && hist[2]==10800);
	}

	// --- THE ACTUAL TEST: does every shared facet drop the SAME rank on
	// both sides, with the ranks used AS-IS (no search, no coloring)? ---
	int nmismatch=0, nchecked=0;
	for(map<FacetKey,vector<pair<int,int> > >::iterator it=owner.begin(); it!=owner.end(); ++it) {
		const vector<pair<int,int> >& own=it->second;
		++nchecked;
		if(own[0].second != own[1].second) ++nmismatch;
	}
	cout<<endl<<"=== THESIS CLAIM TEST ==="<<endl;
	cout<<"shared facets checked: "<<nchecked<<endl;
	cout<<"rank mismatches (should be 0 if the claim holds): "<<nmismatch<<endl;
	if(nmismatch==0)
		cout<<"CONFIRMED: using each flag's own rank-order (0..4) directly as its "
				<<"Maubach vertex order satisfies d_i(sigma1)=d_j(sigma2) => i=j "
				<<"for ALL 10800 shared facets of the barycentric subdivision, with "
				<<"NO coloring/search needed -- the claim holds for CP^2_9."<<endl;
	else
		cout<<"REFUTED for this complex: rank alone is not enough, "<<nmismatch
				<<" facets need a real reordering choice even after subdivision."<<endl;

	// --- Bonus: does this also give a coherent global orientation, using
	// a sign derived from a single FIXED global vertex order (0..8)? ---
	vector<int> sign(M);
	for(int i=0;i<M;++i) {
		int inv=0;
		for(int a=0;a<5;++a) for(int b=a+1;b<5;++b) if(flag[i][a]>flag[i][b]) ++inv;
		sign[i]=(inv%2==0)?+1:-1;
	}
	int nincoherent=0;
	for(map<FacetKey,vector<pair<int,int> > >::iterator it=owner.begin(); it!=owner.end(); ++it) {
		const vector<pair<int,int> >& own=it->second;
		int i=own[0].first, ki=own[0].second, j=own[1].first, kj=own[1].second;
		int s0=(ki%2==0)?1:-1, s1=(kj%2==0)?1:-1;
		if(sign[i]*s0 + sign[j]*s1 != 0) ++nincoherent;
	}
	cout<<endl<<"=== BONUS: global-order permutation sign as orientation ==="<<endl;
	cout<<"incoherently oriented facets: "<<nincoherent<<"/"<<nchecked<<endl;
	if(nincoherent==0)
		cout<<"Barycentric subdivision ALSO gives a coherent orientation for free "
				<<"(via global permutation parity) -- resolves the bipartiteness "
				<<"obstruction found on the un-subdivided complex."<<endl;
	else
		cout<<"Global permutation parity is NOT a coherent orientation here -- "
				<<"the compatibility claim (ranks match) and orientation coherence "
				<<"are separate questions; subdivision fixed the first but not "
				<<"(with this particular sign choice) the second."<<endl;

	// Independent BFS derivation (same algorithm as find_reordering.cpp /
	// riemann.cpp), in case a *different* coherent orientation exists even
	// though the naive global-parity one isn't it.
	{
		vector<int> owner_i, owner_j, owner_k;
		for(map<FacetKey,vector<pair<int,int> > >::iterator it=owner.begin(); it!=owner.end(); ++it) {
			const vector<pair<int,int> >& own=it->second;
			owner_i.push_back(own[0].first);
			owner_j.push_back(own[1].first);
			owner_k.push_back(own[0].second); // ranks match, verified above (if nmismatch==0)
		}
		int E=(int)owner_i.size();
		vector<vector<int> > adj(M);
		for(int e=0;e<E;++e) { adj[owner_i[e]].push_back(e); adj[owner_j[e]].push_back(e); }

		vector<int> bfsSign(M,0);
		vector<char> seen(M,false);
		vector<int> q;
		bfsSign[0]=+1; seen[0]=true; q.push_back(0);
		int nconflict=0;
		for(size_t qi=0; qi<q.size(); ++qi) {
			int u=q[qi];
			for(size_t t=0;t<adj[u].size();++t) {
				int e=adj[u][t];
				int v=(owner_i[e]==u)?owner_j[e]:owner_i[e];
				int k=owner_k[e];
				int s=(k%2==0)?1:-1;
				int required=-bfsSign[u]*s*s; // s*s==1 always (matching ranks)
				if(seen[v]) { if(bfsSign[v]!=required) ++nconflict; continue; }
				bfsSign[v]=required; seen[v]=true; q.push_back(v);
			}
		}
		cout<<endl<<"=== BFS coherent-orientation search (independent of the sign choice above) ==="<<endl;
		cout<<"BFS reached "<<q.size()<<"/"<<M<<" flags, "<<nconflict
				<<" conflicts (0 = a coherent orientation EXISTS for the barycentric "
				<<"subdivision, regardless of which sign convention finds it)"<<endl;
	}

	return 0;
}
