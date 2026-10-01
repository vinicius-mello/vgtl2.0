import itertools, math, numpy as np, pickle
from multiprocessing import Pool
import charac_rule as cr, geo2, lifetime as lt
geo2.CAP = 3000
maub = [(frozenset({0,t}), t) for t in (4,3,2,1)]
faces = [(frozenset(F), l) for F in itertools.combinations(range(5), 3) for l in range(5)]
cands = set()
for pos in range(5):
    for f in faces:
        T = tuple(maub[:pos] + [f] + maub[pos:]); cands.add(T)
    for rep in range(4):
        for f in faces:
            T = list(maub); T[rep] = f; cands.add(tuple(T))
cands = [T for T in cands if cr.rule(T)]
print("valid Maubach+face cycles:", len(cands))
def job(T):
    L = [lt.lifetime(list(T), s) for s in range(len(T))]
    if any(x is None for x in L): return None
    LEV = 3 * len(T)
    o, ex = geo2.run(T, LEV)
    N = np.prod([len(T[i % len(T)][0]) for i in range(LEV)])
    h, q = o[-1]
    return (-math.log(h)/math.log(N), q, max(L), T, ex)
if __name__ == "__main__":
    with Pool() as p: res = [r for r in p.map(job, cands) if r]
    print("with finite edge lifetime:", len(res))
    res.sort(key=lambda r: -r[0])
    for a, q, L, T, ex in res[:12]:
        print(f"alpha={a:.3f} qmin={q:.4f} life={L} {'exact' if ex else 'sampled'}", [(tuple(sorted(F)), l) for F, l in T])
    print("--- best qmin:")
    for a, q, L, T, ex in sorted(res, key=lambda r: -r[1])[:8]:
        print(f"alpha={a:.3f} qmin={q:.4f} life={L} {'exact' if ex else 'sampled'}", [(tuple(sorted(F)), l) for F, l in T])
    pickle.dump(res, open("insert.pkl", "wb"))
