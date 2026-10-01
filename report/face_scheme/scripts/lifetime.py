"""Combinatorial edge lifetime: max number of levels an edge of a cell survives
(stays an edge of some descendant) under uniform application of the type cycle.
Descendants of a cell = recursive children; an edge {x,y} of the root 'dies' in a
descendant branch when a split of an edge-type cell has split face {x,y}."""
import pickle, itertools
d = 4
found = pickle.load(open("hybrid.pkl","rb"))
def children(cell, F, l, v):
    out = []
    for p in sorted(F):
        rest = [cell[r] for r in range(d+1) if r != p]; rest.insert(l, v); out.append(tuple(rest))
    return out
def lifetime(T, start, maxlev=40):
    """max level L such that some level-L descendant of a root of type `start`
    still contains an edge between two root vertices."""
    root = tuple(range(d+1))
    cur = {root}; cnt = 0; last = 0
    for lev in range(maxlev):
        F, l = T[(start+lev) % len(T)]
        nxt = set()
        for c in cur:
            for ch in children(c, F, l, ('v', lev, c)):
                # keep only descendants that still contain >= 2 root vertices
                if sum(1 for x in ch if isinstance(x, int)) >= 2:
                    # canonicalize: only root-vertex positions matter + types
                    nxt.add(tuple(x if isinstance(x, int) else None for x in ch))
        # canonical form collapses new-vertex identities (fine: only positions of root verts matter)
        cur = {tuple(x if x is not None else ('n',) for x in c) for c in nxt}
        if not cur: return lev
    return None  # survives maxlev levels
if __name__ == "__main__":
    ok = []
    for T in found:
        L = [lifetime(T, s) for s in range(len(T))]
        if all(x is not None for x in L): ok.append((max(L), T))
    ok.sort(key=lambda r: r[0])
    print(len(found), "valid hybrids;", len(ok), "with finite edge lifetime")
    from collections import Counter
    print("lifetime distribution:", sorted(Counter(x for x,_ in ok).items()))
    for L, T in ok[:30]: print(L, [(tuple(sorted(F)), l) for F, l in T])
    pickle.dump(ok, open("finite.pkl","wb"))
    # sanity: Maubach and user's scheme
    maub = [(frozenset({0,t}), t) for t in (4,3,2,1)]
    user = [(frozenset({0,1,t}), t) for t in (4,3,2)]
    print("Maubach lifetime:", [lifetime(maub, s) for s in range(4)])
    print("user scheme lifetime:", [lifetime(user, s) for s in range(3)])
