import pickle, itertools, numpy as np
from fractions import Fraction
found = pickle.load(open("hybrid.pkl","rb"))
d = 4
def sig(P):
    E = [P[a]-P[b] for a,b in itertools.combinations(range(d+1),2)]
    L = [int(e@e) for e in E]  # integer coords scaled
    m = max(L)
    return tuple(Fraction(x, m) for x in L)
def children(P, F, l):
    F = sorted(F); v = sum(P[p] for p in F)
    k = len(F)
    Q = [p*k for p in P]  # scale so centroid integer
    out = []
    for p in F:
        rest = [Q[r] for r in range(d+1) if r != p]; rest.insert(l, v); out.append(rest)
    return out
def qual(P):
    E = [np.linalg.norm((P[a]-P[b]).astype(float)) for a,b in itertools.combinations(range(d+1),2)]
    M = np.array([(P[i]-P[0]).astype(float) for i in range(1,d+1)])
    return abs(np.linalg.det(M))/max(E)**d
def run(T, levels):
    P0 = [np.array([1 if i < m else 0 for i in range(d)], dtype=object) for m in range(d+1)]
    cur = {sig(P0): P0}; q0 = qual([p for p in P0]); hist = []
    for lev in range(levels):
        F, l = T[lev % len(T)]
        nxt = {}
        for P in cur.values():
            for C in children(P, F, l):
                g = C[0]*0
                # normalize to keep integers small: divide by gcd
                import math
                gg = 0
                for c in C:
                    for x in c: gg = math.gcd(gg, int(x))
                C = [c // gg for c in C] if gg > 1 else C
                s = sig(C)
                if s not in nxt: nxt[s] = C
        cur = nxt
        qs = min(qual(C) for C in cur.values())/q0
        hist.append((len(cur), qs))
        if len(cur) > 400: break
    return hist
res = []
for T in found:
    h = run(T, 4*len(T) if len(T)==3 else 12)
    res.append((h[-1][0], -h[-1][1], T, h))
res.sort(key=lambda r: (r[0], r[1]))
good = [r for r in res if r[0] <= 50]
print("candidates:", len(res), " with <=50 similarity classes at last level:", len(good))
for n, mq, T, h in res[:25]:
    print(n, f"minq={-mq:.3g}", [(tuple(sorted(F)), l) for F, l in T], [x[0] for x in h])
