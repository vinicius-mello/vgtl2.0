"""Exhaustive symbolic check of the 'two cells' + 'children' lemmas for the
prefix-face scheme: F(t) = {0..k-1} u {t}, t in {k..d}, new vertex at position t,
child type = t-1 (wrapping k-1 -> d)."""
import itertools, sys

def F(k, t): return set(range(k)) | {t}
def nxt(k, d, t): return t - 1 if t > k else d

def children(k, d, cell, t, v):
    P = sorted(F(k, t))
    out = []
    for p in P:
        rest = [cell[r] for r in range(d + 1) if r != p]
        rest.insert(t, v)
        out.append(tuple(rest))
    return out

def pair_ok(k, d, A, B, weak):
    (va, la, ta), (vb, lb, tb) = A, B
    common = set(va) & set(vb)
    if len(common) < d: return None           # not adjacent
    i = next(p for p in range(d+1) if va[p] not in common)
    j = next(p for p in range(d+1) if vb[p] not in common)
    if va[:i]+va[i+1:] != vb[:j]+vb[j+1:]: return "unordered facet"
    if abs(la-lb) > 1: return "OneLevel"
    if la == lb:
        assert ta == tb
        if i == j: return True
        ok = {i, j} <= F(k, ta) if weak else None
        return True if ok else f"A: t={ta} i={i} j={j}"
    if la < lb: (i, j, ta, tb) = (j, i, tb, ta)   # now a is finer
    # a finer (type ta), b coarser (type tb)
    if i == tb and j in F(k, tb): return True
    return f"B: coarse t={tb} fine i={i} coarse j={j}"

def check(d, k, weak=True):
    bad = []; ncases = 0
    types = list(range(k, d + 1))
    W = tuple(f"w{p}" for p in range(d + 1))
    for t in types:
        # siblings
        ch = children(k, d, W, t, "v")
        cs = [(c, 1, nxt(k, d, t)) for c in ch]
        for a, b in itertools.combinations(cs, 2):
            r = pair_ok(k, d, a, b, weak)
            if r is not True: bad.append(("siblings", t, r))
        configs = []
        for i in range(d + 1):
            for j in range(d + 1):
                if i == j or {i, j} <= F(k, t):
                    configs.append(("eq", i, j))
                if i == t and j in F(k, t):
                    configs.append(("fine", i, j))
        for kind, i, j in configs:
            ncases += 1
            facet = W[:j] + W[j+1:]
            S = facet[:i] + ("x",) + facet[i:]
            if kind == "eq":
                Scell = (S, 0, t)
            else:
                Scell = (S, 1, nxt(k, d, t))
            # sanity: original pair satisfies invariants
            assert pair_ok(k, d, (W, 0, t), Scell, weak) is True, (kind, i, j)
            face = {W[p] for p in F(k, t)}
            new = [(c, 1, nxt(k, d, t)) for c in children(k, d, W, t, "v")]
            other = [Scell]
            if face <= set(S):
                if kind != "eq" or {S[p] for p in F(k, t)} != face:
                    bad.append((kind, t, i, j, "S contains face of W but split face differs")); continue
                other = [(c, 1, nxt(k, d, t)) for c in children(k, d, S, t, "v")]
            for a in new:
                for b in other:
                    r = pair_ok(k, d, a, b, weak)
                    if r is None or r is True: continue
                    bad.append((kind, t, i, j, r))
    return ncases, bad

for d in range(2, 8):
    for k in range(1, d):
        n, bad = check(d, k)
        print(f"d={d} k={k}: {n} configurations, {'OK' if not bad else 'FAIL '+str(bad[:3])}")
