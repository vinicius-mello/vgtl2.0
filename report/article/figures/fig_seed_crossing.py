"""Figure for Section 6.4: two colour-ordered tetrahedra (balanced seed) glued
along facet 3, adaptively Maubach-refined with conforming closure near a
point of the shared facet.  The leaf c (in tau) with a facet on the shared
facet and its neighbour c' across it (in tau') are highlighted; the script
asserts that c' has the same bisection code (path) as c, i.e. Prop. 6.2.

Style follows the line drawings of Medeiros e Sa et al. (Adaptive Voids,
2015): thin black lines on white, one grey highlight, dashed auxiliary lines.
Output: seed_crossing.pdf
"""
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","mathtext.fontset":"cm","font.serif":["DejaVu Serif"]})
from matplotlib.patches import Polygon
from itertools import combinations

# ---------------------------------------------------------------- geometry
V = []  # vertex coordinates
def addv(p):
    V.append(np.asarray(p, float)); return len(V) - 1

s3 = np.sqrt(3.0)
v0 = addv([-1.0, -s3/3, 0.0])
v1 = addv([ 1.0, -s3/3, 0.0])
v2 = addv([ 0.0, 2*s3/3, 0.0])
h = 1.55
v3 = addv([0.15, 0.05,  h])     # tau  (colour 3)
v3p = addv([-0.1, 0.1, -h])     # tau' (colour 3)

M = 3  # dimension
# cell -> (level, seed, path)
cells = {}
cells[(v0, v1, v2, v3)] = (0, "tau", "")
cells[(v0, v1, v2, v3p)] = (0, "tau'", "")
midpoint = {}

def ctype(c):
    return M - (cells[c][0] % M)

def facets_of(c):
    return [(i, frozenset(c[:i] + c[i+1:])) for i in range(M + 1)]

def adjacent(c, i):
    f = frozenset(c[:i] + c[i+1:])
    for d in cells:
        if d != c and f <= set(d):
            return d
    return None

def split_cell(c, k):
    lvl, seed, path = cells.pop(c)
    a, b = c[0], c[k]
    key = frozenset((a, b))
    if key not in midpoint:
        midpoint[key] = addv((V[a] + V[b]) / 2)
    v = midpoint[key]
    c0 = c[:k] + (v,) + c[k+1:]
    c1 = c[1:k+1] + (v,) + c[k+1:]
    cells[c0] = (lvl + 1, seed, path + "0")
    cells[c1] = (lvl + 1, seed, path + "1")

def maubach_subdivide(sigma):
    """Algorithm 1 of the article (thesis Listing 4.2)."""
    k = ctype(sigma)
    st, visited, queue = [], set(), [sigma]
    while queue:
        c = queue.pop(0)
        if c in visited:
            continue
        visited.add(c); st.append(c)
        for i in range(M + 1):
            if i == 0 or i == k:
                continue
            ca = adjacent(c, i)
            if ca is None or ca in visited:
                continue
            if cells[ca][0] < cells[c][0]:
                maubach_subdivide(ca)
                ca = adjacent(c, i)
            queue.append(ca)
    edge = frozenset((sigma[0], sigma[k]))
    for c in st:
        kc = ctype(c)
        assert frozenset((c[0], c[kc])) == edge, "invariant violated"
        split_cell(c, kc)

def bary(c, p):
    A = np.array([V[c[i]] - V[c[3]] for i in range(3)]).T
    w = np.linalg.solve(A, p - V[c[3]])
    return np.append(w, 1 - w.sum())

def leaf_containing(p):
    for c in cells:
        if (bary(c, p) >= -1e-12).all():
            return c

# target: a point just above the shared facet, off-centre
p = np.array([0.22, 0.12, 0.02])
L = 6
while True:
    c = leaf_containing(p)
    if cells[c][0] >= L:
        break
    maubach_subdivide(c)

c = leaf_containing(p)
shared = {v0, v1, v2}
on_plane = [i for i in range(4) if all(abs(V[x][2]) < 1e-12 for j, x in enumerate(c) if j != i)]
assert on_plane, "chosen leaf has no facet on the shared facet; move p"
i_b = on_plane[0]
cp = adjacent(c, i_b)
print("c :", cells[c], " facet index", i_b)
print("c':", cells[cp])
assert cells[c][1] == "tau" and cells[cp][1] == "tau'"
assert cells[c][2] == cells[cp][2], "Prop. 6.2 violated: codes differ"
# same local index on both sides (ordered facets equal)
assert c[:i_b] + c[i_b+1:] == cp[:i_b] + cp[i_b+1:]
print("leaves:", len(cells), "  path of c = path of c' =", cells[c][2])

# ---------------------------------------------------------------- drawing
def proj(P, elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    R = np.array([[np.cos(a), np.sin(a), 0],
                  [-np.sin(a)*np.sin(e), np.cos(a)*np.sin(e), np.cos(e)]])
    return (R @ np.asarray(P).T).T

def view_dir(elev, azim):
    # unit vector pointing from the scene towards the viewer
    e, a = np.radians(elev), np.radians(azim)
    r0 = np.array([np.cos(a), np.sin(a), 0])
    r1 = np.array([-np.sin(a)*np.sin(e), np.cos(a)*np.sin(e), np.cos(e)])
    n = np.cross(r0, r1)
    return n if n[2] * np.sin(e) >= 0 else -n

def painter(ax, faces, elev, azim, lw, zbase=3):
    # faces: list of (3 points, facecolor); drawn far-to-near
    n = view_dir(elev, azim)
    faces = sorted(faces, key=lambda t: np.mean([p @ n for p in t[0]]))
    for k, (pts, fc) in enumerate(faces):
        ax.add_patch(Polygon(proj(pts, elev, azim), closed=True, fc=fc, ec="k",
                             lw=lw, alpha=0.95, zorder=zbase + 1e-3 * k))

def edges_of(cs):
    E = set()
    for c in cs:
        for a, b in combinations(c, 2):
            E.add((min(a, b), max(a, b)))
    return E

def draw(ax, offset_tau, offset_taup, elev, azim, labels=True):
    def pos(x, seed):
        o = offset_tau if seed == "tau" else offset_taup
        return V[x] + o
    # shared facet: light grey fill (drawn for each copy)
    for seed in (["tau", "tau'"] if np.any(offset_tau != offset_taup) else ["tau"]):
        tri = proj([pos(x, seed) for x in (v0, v1, v2)], elev, azim)
        ax.add_patch(Polygon(tri, closed=True, fc="0.93", ec="none", zorder=0))
    # highlighted cells: faces filled
    faces = []
    for cc, fc in ((c, "0.55"), (cp, "0.78")):
        seed = cells[cc][1]
        for f in combinations(cc, 3):
            faces.append(([pos(x, seed) for x in f], fc))
    painter(ax, faces, elev, azim, 0.5)
    # all leaf edges
    for seed, col in (("tau", "k"), ("tau'", "0.35")):
        cs = [x for x in cells if cells[x][1] == seed]
        for a, b in edges_of(cs):
            P = proj([pos(a, seed), pos(b, seed)], elev, azim)
            ax.plot(P[:, 0], P[:, 1], color=col, lw=0.25, zorder=1)
    # seed tetrahedra outlines, dashed
    for seed, top in (("tau", v3), ("tau'", v3p)):
        for a, b in combinations((v0, v1, v2, top), 2):
            P = proj([pos(a, seed), pos(b, seed)], elev, azim)
            ax.plot(P[:, 0], P[:, 1], "k--", lw=0.6, dashes=(4, 2), zorder=2)
    if labels:
        for x, name in ((v0, "$v_0$"), (v1, "$v_1$"), (v2, "$v_2$")):
            q = proj([pos(x, "tau")], elev, azim)[0]
            ax.text(q[0] + 0.05, q[1] + 0.03, name, fontsize=8)
        q = proj([pos(v3, "tau")], elev, azim)[0]; ax.text(q[0] + 0.05, q[1], "$v_3$", fontsize=8)
        q = proj([pos(v3p, "tau'")], elev, azim)[0]; ax.text(q[0] + 0.05, q[1], "$v_3'$", fontsize=8)
        q = proj([pos(v3, "tau") * 0.5 + pos(v0, "tau") * 0.5], elev, azim)[0]
        ax.text(q[0] - 0.25, q[1], r"$\tau$", fontsize=10)
        q = proj([pos(v3p, "tau'") * 0.5 + pos(v0, "tau'") * 0.5], elev, azim)[0]
        ax.text(q[0] - 0.25, q[1], r"$\tau'$", fontsize=10)
    ax.set_aspect("equal"); ax.axis("off"); ax.margins(0.06)

elev, azim = 26, -58
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.6), gridspec_kw=dict(width_ratios=[1, 1, 0.85]))
z = np.zeros(3)
draw(axs[0], z, z, elev, azim)
axs[0].set_title("(a)", fontsize=9, y=-0.08)
n = np.array([0, 0, 1.0])
draw(axs[1], 0.6 * n, -0.6 * n, elev, azim)
axs[1].set_title("(b)", fontsize=9, y=-0.08)

# (c): close-up of the two highlighted cells, exploded, with local indices
ax = axs[2]
def show_cell(cc, off, fc):
    painter(ax, [([V[x] + off for x in f], fc) for f in combinations(cc, 3)], elev, azim, 0.6, zbase=1)
    for j, x in enumerate(cc):
        q = proj([V[x] + off], elev, azim)[0]
        ax.plot(q[0], q[1], "ko", ms=2)
        ax.text(q[0] + 0.004, q[1] + 0.004, str(j), fontsize=7)
ctr = np.mean([V[x] for x in c], axis=0)
sc = max(np.ptp([V[x] for x in c], axis=0))
show_cell(c, 0.6 * sc * n, "0.55")
show_cell(cp, -0.6 * sc * n, "0.78")
f = [x for j, x in enumerate(c) if j != i_b]
for x in f:  # dashed correspondences between the shared facet copies
    P = proj([V[x] + 0.6 * sc * n, V[x] - 0.6 * sc * n], elev, azim)
    ax.plot(P[:, 0], P[:, 1], "k--", lw=0.4, dashes=(2, 2))
qa=proj([ctr + 0.6 * sc * n], elev, azim)[0]; qb=proj([ctr - 0.6 * sc * n], elev, azim)[0]
ax.text(qa[0]-0.45*sc, qa[1], "$c$", fontsize=11)
ax.text(qb[0]-0.45*sc, qb[1], "$c'$", fontsize=11)
ax.set_aspect("equal"); ax.axis("off"); ax.margins(0.06)
ax.set_title("(c)", fontsize=9, y=-0.08)

plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.06, wspace=0.02)
plt.savefig("seed_crossing.pdf")
print("wrote seed_crossing.pdf")
