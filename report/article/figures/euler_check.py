"""Genus check of an extracted Riemann-surface mesh (alpha projection,
vertices shared by node identity).

Bad cells leave holes. Before counting, every non-manifold ("pinched")
vertex, i.e. a vertex whose link consists of k>1 components, is split into
k copies, one per link component. On the resulting surface with boundary,
if every hole is a disc,
    chi = sum over components of (2 - 2 g_i) - b,
so for a connected surface g = (2 - b - chi) / 2.

Usage: python3 euler_check.py MESH.obj [expected_genus]
"""
import sys
from collections import defaultdict


def components(nodes, adj):
    seen, comps = set(), []
    for x in nodes:
        if x in seen:
            continue
        comp, st = [], [x]
        while st:
            y = st.pop()
            if y in seen:
                continue
            seen.add(y); comp.append(y); st.extend(adj[y])
        comps.append(comp)
    return comps


faces = []
for line in open(sys.argv[1]):
    if line.startswith("f "):
        faces.append([int(t.split("/")[0]) - 1 for t in line.split()[1:]])

# ---- split pinched vertices: corners (face, position) grouped by link component
link = defaultdict(list)            # v -> list of (a, b, face_idx, pos)
for fi, f in enumerate(faces):
    n = len(f)
    for i, v in enumerate(f):
        link[v].append((f[i - 1], f[(i + 1) % n], fi, i))
newid, pinched, nid = {}, 0, 0
for v, es in link.items():
    ladj = defaultdict(set)
    for a, b, _, _ in es:
        ladj[a].add(b); ladj[b].add(a)
    comps = components(list(ladj), ladj)
    which = {x: k for k, comp in enumerate(comps) for x in comp}
    if len(comps) > 1:
        pinched += 1
    base = nid; nid += len(comps)
    for a, b, fi, i in es:
        newid[(fi, i)] = base + which[a]
faces = [[newid[(fi, i)] for i in range(len(f))] for fi, f in enumerate(faces)]

# ---- counts on the split surface
edges = defaultdict(int)
for f in faces:
    for a, b in zip(f, f[1:] + f[:1]):
        edges[(min(a, b), max(a, b))] += 1
V, E, F = nid, len(edges), len(faces)
chi = V - E + F
nonmanifold = sum(1 for c in edges.values() if c > 2)
badj = defaultdict(list)
for (a, b), c in edges.items():
    if c == 1:
        badj[a].append(b); badj[b].append(a)
loops = components(list(badj), badj)
b = len(loops)
sadj = defaultdict(set)
for f in faces:
    for a in f[1:]:
        sadj[f[0]].add(a); sadj[a].add(f[0])
ncomp = len(components(list(range(V)), sadj))

g = (2 * ncomp - b - chi) / 2
print(f"pinched vertices split: {pinched}; non-manifold edges: {nonmanifold}; "
      f"components: {ncomp}")
print(f"V={V} E={E} F={F} chi={chi} boundary loops b={b} "
      f"(largest {sorted(len(l) for l in loops)[-3:]})")
print(f"genus estimate g = (2c - b - chi)/2 = {g}" +
      (f"   (expected {sys.argv[2]})" if len(sys.argv) > 2 else ""))
