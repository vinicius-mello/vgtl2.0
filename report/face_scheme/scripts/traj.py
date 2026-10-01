import geo2, math, numpy as np
geo2.CAP = 20000
S = {
 "Maubach": [(frozenset({0,4}),4),(frozenset({0,3}),3),(frozenset({0,2}),2),(frozenset({0,1}),1)],
 "M+face(l=2)": [(frozenset({0,4}),4),(frozenset({0,3}),3),(frozenset({0,1,2}),2),(frozenset({0,1}),1)],
 "M+face(l=3)": [(frozenset({0,4}),4),(frozenset({0,3}),3),(frozenset({0,1,2}),3),(frozenset({0,1}),1)],
}
for name, T in S.items():
    o, ex = geo2.run(T, 24)
    N = 1; row = []
    for i, (h, q) in enumerate(o):
        N *= len(T[i % 4][0])
        if (i+1) % 4 == 0: row.append(f"c{(i+1)//4}: h={h:.3f} q={q:.4f} N={N}")
    print(name, "exact" if ex else "sampled"); print("   " + "\n   ".join(row))
