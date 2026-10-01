import pickle, numpy as np, math, sys
from multiprocessing import Pool
import geo2
LEV = int(sys.argv[1]); geo2.CAP = int(sys.argv[2])
def job(item):
    L, T = item
    o, ex = geo2.run(T, LEV)
    N = np.prod([len(T[i % len(T)][0]) for i in range(LEV)])   # cells per root after LEV levels
    h, q = o[-1]
    alpha = -math.log(h) / math.log(N) if h < 1 else 0.0        # h ~ N^-alpha ; ideal 1/d = 0.25
    return (alpha, q, L, T, o, ex)
if __name__ == "__main__":
    ok = [x for x in pickle.load(open("finite.pkl", "rb"))]
    maub = [(frozenset({0,t}), t) for t in (4,3,2,1)]
    ok.append((4, maub))
    with Pool() as p: res = p.map(job, ok, chunksize=4)
    pickle.dump(res, open(f"screen_{LEV}.pkl", "wb"))
    m = [r for r in res if r[3] == maub][0]
    print(f"Maubach: alpha={m[0]:.3f} qmin={m[1]:.3f}")
    res.sort(key=lambda r: -(r[0] + 0.0*r[1]))
    for a, q, L, T, o, ex in res[:15]:
        print(f"alpha={a:.3f} qmin={q:.4f} life={L} {'exact' if ex else 'sampled'}", [(tuple(sorted(F)), l) for F, l in T])
    print("--- best quality among alpha >= 0.2:")
    for a, q, L, T, o, ex in sorted([r for r in res if r[0] >= 0.2], key=lambda r: -r[1])[:15]:
        print(f"alpha={a:.3f} qmin={q:.4f} life={L}", [(tuple(sorted(F)), l) for F, l in T])
