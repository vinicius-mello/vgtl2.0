"""Search cyclic schemes in d=4 whose types mix edges and 2-faces:
type = (F, l). Invariants: (A) eq level: i=j or {i,j} subset F; (B) fine i = l(coarse), coarse j in F(coarse)."""
import itertools
d = 4
def children(cell, F, l, v):
    out = []
    for p in sorted(F):
        rest = [cell[r] for r in range(d+1) if r != p]; rest.insert(l, v); out.append(tuple(rest))
    return out
def pair(A, B, T):
    (va, la, ta), (vb, lb, tb) = A, B
    common = set(va) & set(vb)
    if len(common) < d: return True
    i = next(p for p in range(d+1) if va[p] not in common); j = next(p for p in range(d+1) if vb[p] not in common)
    if va[:i]+va[i+1:] != vb[:j]+vb[j+1:]: return False
    if la == lb: return i == j or {i, j} <= T[ta][0]
    if la < lb: i, j, ta, tb = j, i, tb, ta
    return i == T[tb][1] and j in T[tb][0]
def check(T):
    n = len(T); W = tuple(f"w{p}" for p in range(d+1))
    for t in range(n):
        F, l = T[t]; c = (t+1) % n
        if l in T[c][0]: return False           # alternating / direction
        ch = [(x, 1, c) for x in children(W, F, l, "v")]
        if not all(pair(a, b, T) for a, b in itertools.combinations(ch, 2)): return False
        for i in range(d+1):
            for j in range(d+1):
                cfg = []
                if i == j or {i, j} <= F: cfg.append(0)
                if i == l and j in F: cfg.append(1)
                for lev in cfg:
                    fac = W[:j]+W[j+1:]; S = fac[:i]+("x",)+fac[i:]
                    Sc = (S, lev, c if lev else t); other = [Sc]
                    face = {W[p] for p in F}
                    if face <= set(S):
                        if lev or {S[p] for p in F} != face: return False
                        other = [(x, 1, c) for x in children(S, F, l, "v")]
                    for a in ch:
                        for b in other:
                            if not pair(a, b, T): return False
    return True
types = [(frozenset(F), l) for s in (2, 3) for F in itertools.combinations(range(d+1), s) for l in range(d+1)]
found = []
for n in (2, 3):
    for T in itertools.product(types, repeat=n):
        if T[0] != min(T[r:] + T[:r] for r in range(n))[0] and False: pass
        sizes = {len(F) for F, _ in T}
        if sizes != {2, 3}: continue
        if check(T): found.append(T)
    print("cycle length", n, "hybrid schemes passing:", len(found))
import pickle; pickle.dump(found, open("hybrid.pkl", "wb"))
for T in found[:40]: print([ (tuple(sorted(F)), l) for F, l in T])
