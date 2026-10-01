"""Geometry on Kuhn's simplex, midpoint/centroid placement, uniform refinement.
Dedup by similarity (ordered edge-length signature), keeping the largest scale:
exact for h_max and q_min as long as the class count stays under CAP."""
import pickle, itertools, numpy as np, random
d = 4; CAP = 6000
PAIRS = list(itertools.combinations(range(d+1), 2))
def stats(P):
    E = np.array([np.linalg.norm(P[a]-P[b]) for a, b in PAIRS]); h = E.max()
    vol = abs(np.linalg.det(P[1:]-P[0]))
    return h, vol/h**d, tuple(np.round(E/h, 9))
def run(T, levels, start=0):
    P0 = np.array([[1.0 if i < m else 0.0 for i in range(d)] for m in range(d+1)])
    h0, q0, s0 = stats(P0)
    cur = {s0: P0}; out = []; exact = True
    for lev in range(levels):
        F, l = T[(start+lev) % len(T)]; F = sorted(F)
        nxt = {}
        for P in cur.values():
            v = P[F].mean(axis=0)
            for p in F:
                C = np.insert(np.delete(P, p, axis=0), l, v, axis=0)
                h, q, s = stats(C)
                if s not in nxt or stats(nxt[s])[0] < h: nxt[s] = C
        if len(nxt) > CAP:
            exact = False; keys = random.Random(0).sample(list(nxt), CAP); nxt = {k: nxt[k] for k in keys}
        cur = nxt
        H = max(stats(P)[0] for P in cur.values()); Q = min(stats(P)[1] for P in cur.values())
        out.append((H/h0, Q/q0))
    return out, exact
if __name__ == "__main__":
    import sys
    ok = pickle.load(open("finite.pkl", "rb"))
    maub = [(frozenset({0,t}), t) for t in (4,3,2,1)]
    LEV = 12
    m, _ = run(maub, LEV)
    print("Maubach  h:", " ".join(f"{h:.3f}" for h, q in m)); print("         q:", " ".join(f"{q:.3f}" for h, q in m))
    res = []
    for L, T in ok:
        if L > 5: continue
        o, ex = run(T, LEV)
        res.append((o[-1][0], -o[-1][1], L, T, o, ex))
    # rank by cells-per-halving would be fairer; report h and q at level 12 plus growth factor
    res.sort(key=lambda r: (round(r[0], 3), r[1]))
    for hL, mq, L, T, o, ex in res[:12]:
        g = np.prod([len(F) for F, _ in T])**(1/len(T))
        print(f"life={L} growth/level={g:.2f} h12={hL:.3f} q12={-mq:.3g} {'exact' if ex else 'sampled'}",
              [(tuple(sorted(F)), l) for F, l in T])
        print("    q:", " ".join(f"{q:.3f}" for h, q in o))
    pickle.dump(res, open("geo2.pkl", "wb"))
