"""Topological post-processing of an extracted Riemann-surface mesh.

1. split pinched vertices (one copy per link component);
2. keep the main component (islands cut off by pinches are dropped);
3. cap every boundary loop (hole left by bad cells) with a cone from a new
   centre vertex (a fan from an existing vertex could duplicate an edge);
4. verify: every edge in exactly 2 faces, every vertex link one cycle,
   connected, and chi = 2 - 2g with g = (n-1)(n-2)/2.

Usage: python3 close_surface.py IN.obj DEGREE [OUT.obj]
(The centre vertex is placed at the loop's centroid in the OBJ coordinates;
placing it on the curve by Newton's method in CP^2 is left to the program.)
"""
import sys
from collections import defaultdict


def comps(nodes, adj):
    seen, out = set(), []
    for x in nodes:
        if x in seen:
            continue
        st, c = [x], []
        while st:
            y = st.pop()
            if y in seen:
                continue
            seen.add(y); c.append(y); st.extend(adj[y])
        out.append(c)
    return out


pts, faces = [], []
for line in open(sys.argv[1]):
    if line.startswith("v "):
        pts.append([float(t) for t in line.split()[1:4]])
    elif line.startswith("f "):
        faces.append([int(t.split("/")[0]) - 1 for t in line.split()[1:]])
n = int(sys.argv[2]); genus = (n - 1) * (n - 2) // 2

# 1. split pinched vertices
corners = defaultdict(list)
for fi, f in enumerate(faces):
    for i, v in enumerate(f):
        corners[v].append((fi, i))
newid, P, npinch = {}, [], 0
for v, cs in corners.items():
    ladj = defaultdict(list)
    for fi, i in cs:
        f = faces[fi]; a, b = f[i - 1], f[(i + 1) % len(f)]
        ladj[a].append(b); ladj[b].append(a)
    cc = comps(list(ladj), ladj)
    which = {x: k for k, c in enumerate(cc) for x in c}
    npinch += len(cc) > 1
    base = len(P)
    P.extend([pts[v]] * len(cc))
    for fi, i in cs:
        newid[(fi, i)] = base + which[faces[fi][i - 1]]
F = [[newid[(fi, i)] for i in range(len(f))] for fi, f in enumerate(faces)]

# 2. keep the main component
sadj = defaultdict(list)
for f in F:
    for a in f[1:]:
        sadj[f[0]].append(a); sadj[a].append(f[0])
cc = comps(list(range(len(P))), sadj)
main = set(max(cc, key=len))
dropped_faces = sum(1 for f in F if f[0] not in main)
F = [f for f in F if f[0] in main]

# 3. cap boundary loops with cones.  Boundary edges are oriented as in
# their face, reversed, so the cone faces match the surface orientation
# when the input polygons are consistently oriented (otherwise only the
# unoriented topology is guaranteed).
ecount = defaultdict(int)
for f in F:
    for a, b in zip(f, f[1:] + f[:1]):
        ecount[(min(a, b), max(a, b))] += 1
bnd = []
for f in F:
    for a, b in zip(f, f[1:] + f[:1]):
        if ecount[(min(a, b), max(a, b))] == 1:
            bnd.append((b, a))
badj = defaultdict(list)
for a, b in bnd:
    badj[a].append(b); badj[b].append(a)
loops = comps(list(badj), badj)
loop_of = {x: k for k, c in enumerate(loops) for x in c}
centre = {}
for k, c in enumerate(loops):
    centre[k] = len(P)
    P.append([sum(P[x][j] for x in c) / len(c) for j in range(3)])
for a, b in bnd:
    F.append([a, b, centre[loop_of[a]]])

# 4. verify
ecount = defaultdict(int)
for f in F:
    for a, b in zip(f, f[1:] + f[:1]):
        ecount[(min(a, b), max(a, b))] += 1
used = sorted({x for f in F for x in f})
V, E, NF = len(used), len(ecount), len(F)
chi = V - E + NF
bad_edges = sum(1 for c in ecount.values() if c != 2)
link = defaultdict(list)
for f in F:
    for i, v in enumerate(f):
        link[v].append((f[i - 1], f[(i + 1) % len(f)]))
bad_links = 0
for v, es in link.items():
    ladj = defaultdict(list)
    for a, b in es:
        ladj[a].append(b); ladj[b].append(a)
    if len(comps(list(ladj), ladj)) != 1 or any(len(ladj[x]) != 2 for x in ladj):
        bad_links += 1
sadj = defaultdict(list)
for f in F:
    for a in f[1:]:
        sadj[f[0]].append(a); sadj[a].append(f[0])
ncomp = len(comps(used, sadj))
print(f"pinched split: {npinch}; islands dropped: {len(cc) - 1} "
      f"({dropped_faces} faces); holes capped: {len(loops)}")
print(f"closed surface check: edges not in exactly 2 faces = {bad_edges}, "
      f"vertices whose link is not one cycle = {bad_links}, components = {ncomp}")
print(f"V={V} E={E} F={NF} chi={chi}; expected chi = 2-2g = {2 - 2 * genus} (g={genus})  "
      + ("OK" if (chi == 2 - 2 * genus and bad_edges == 0 and bad_links == 0 and ncomp == 1) else "MISMATCH"))

if len(sys.argv) > 3:
    remap = {x: k for k, x in enumerate(used)}
    with open(sys.argv[3], "w") as out:
        for x in used:
            out.write("v %.9g %.9g %.9g\n" % tuple(P[x]))
        for f in F:
            out.write("f " + " ".join(str(remap[x] + 1) for x in f) + "\n")
