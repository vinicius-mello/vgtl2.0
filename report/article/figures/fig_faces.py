"""Chart-flat vs. geodesic realization of faces (figure for Sec. "Realizing
the faces").  Drawn on a projective line CP^1 in CP^2, which with the
Fubini-Study metric is a round sphere of radius 1/2 (the Bloch sphere);
a triangle whose vertices lie on that line lies on it entirely, so it can
be drawn without loss.  Two triangles ABC and ABD share the edge AB.

(a) Each triangle is flat in its own best-conditioned chart, chosen by the
    same rule as the code (max over k of min |p_k| over the vertices):
    ABC in chart 0 (p_0 = 1), ABD in chart 1 (p_1 = 1).  An edge flat in
    chart k is an arc of the circle through that chart's point at
    infinity, so the two versions of AB differ and enclose a lens (the
    "seam").  (b) Geodesic cones: both triangles share the great-circle
    arc AB.  Asserts the chart choices and that the geodesic edge does not depend on
the endpoint it is traversed from.
Output: faces.pdf"""
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6})
C0, C1, CG = "#2a78d6", "#eb6834", "#1b8a5a"


def nrm(p):
    return p / np.linalg.norm(p)


def pt(t, phi):  # [cos t : e^{i phi} sin t]
    return nrm(np.array([np.cos(t), np.exp(1j * phi) * np.sin(t)]))


def hdot(a, b):  # <a,b> = sum a_k conj(b_k)
    return np.sum(a * np.conj(b))


def align(a, b):  # b rephased so that <a,b> >= 0
    ip = hdot(a, b)
    return b * ip / abs(ip)


def best_chart(ps):
    return int(np.argmax([min(abs(p[k]) for p in ps) for k in range(2)]))


def flat(p, q, k, lam):  # affine segment in chart k
    return nrm((1 - lam) * p / p[k] + lam * q / q[k])


def geo(p, q, lam):  # geodesic, linear parameter
    return nrm((1 - lam) * p + lam * align(p, q))


def bloch(p):  # CP^1 -> sphere of radius 1/2
    u, v = p
    w = 2 * u * np.conj(v)
    return 0.5 * np.array([w.real, w.imag, abs(u) ** 2 - abs(v) ** 2])


def fs(p, q):
    return np.arccos(min(1.0, abs(hdot(p, q))))


A = pt(np.pi / 4 + 0.02, -0.55)
B = pt(np.pi / 4 - 0.02, 0.55)
Cv = pt(np.pi / 4 - 0.45, 0.0)
D = pt(np.pi / 4 + 0.45, 0.05)
assert best_chart([A, B, Cv]) == 0 and best_chart([A, B, D]) == 1
lam = np.linspace(0, 1, 120)

# orthographic view centred on the two triangles
w = nrm(sum(bloch(p) for p in (A, B, Cv, D)))
e1 = nrm(np.cross([0, 0, 1], w))
e2 = np.cross(w, e1)


def proj(X):
    X = np.atleast_2d(X)
    return np.c_[X @ e1, X @ e2], X @ w


def curve(fun, p, q):
    return np.array([bloch(fun(p, q, l)) for l in lam])


def edges(fun, verts):
    out = []
    for i in range(3):
        out.append(curve(fun, verts[i], verts[(i + 1) % 3]))
    return out


def draw_sphere(ax):
    ax.add_patch(Circle((0, 0), 0.5, fill=False, lw=0.6, color="0.55"))
    for lat in np.linspace(-60, 60, 5):  # faint graticule, front half only
        t = np.radians(lat)
        ring = np.array([[0.5 * np.cos(t) * np.cos(s), 0.5 * np.cos(t) * np.sin(s), 0.5 * np.sin(t)]
                         for s in np.linspace(0, 2 * np.pi, 240)])
        xy, d = proj(ring)
        xy[d < 0] = np.nan
        ax.plot(xy[:, 0], xy[:, 1], lw=0.3, color="0.85")
    for lon in np.linspace(0, np.pi, 7)[:-1]:
        ring = np.array([[0.5 * np.cos(s) * np.cos(lon), 0.5 * np.cos(s) * np.sin(lon), 0.5 * np.sin(s)]
                         for s in np.linspace(0, 2 * np.pi, 240)])
        xy, d = proj(ring)
        xy[d < 0] = np.nan
        ax.plot(xy[:, 0], xy[:, 1], lw=0.3, color="0.85")


def draw_triangle(ax, E, color, fill=True, lw=1.2, ls="-"):
    poly = np.vstack(E)
    xy, _ = proj(poly)
    if fill:
        ax.add_patch(Polygon(xy, closed=True, fc=color, ec="none", alpha=0.18))
    for e in E:
        exy, _ = proj(e)
        ax.plot(exy[:, 0], exy[:, 1], color=color, lw=lw, ls=ls, solid_capstyle="round")


def label_vertices(ax):
    for name, p, off in (("A", A, (0.015, -0.005)), ("B", B, (-0.05, -0.005)),
                         ("C", Cv, (0.0, 0.02)), ("D", D, (0.015, -0.045))):
        xy, _ = proj(bloch(p))
        ax.plot(*xy[0], "o", ms=3.2, color="0.15", zorder=5)
        ax.annotate(f"${name}$", xy[0], xytext=(xy[0][0] + off[0], xy[0][1] + off[1]), fontsize=10)


def full_chart_line(p, q, k):  # the whole line through p,q in chart k (a circle on the sphere)
    ts = np.tan(np.linspace(-np.pi / 2 + 1e-3, np.pi / 2 - 1e-3, 600))
    return np.array([bloch(flat(p, q, k, t)) for t in ts])


fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.3))
for ax in axes:
    ax.set_aspect("equal")
    ax.axis("off")
    draw_sphere(ax)

# (a) chart-flat faces
ax = axes[0]
Eabc = edges(lambda p, q, l: flat(p, q, 0, l), [A, B, Cv])
Eabd = edges(lambda p, q, l: flat(p, q, 1, l), [A, B, D])
for k, col in ((0, C0), (1, C1)):  # the full chart lines carrying the two versions of AB
    xy, d = proj(full_chart_line(A, B, k))
    xy[d < 0] = np.nan
    ax.plot(xy[:, 0], xy[:, 1], color=col, lw=0.6, ls=(0, (1, 2)))
lens = np.vstack([Eabc[0], Eabd[0][::-1]])
ax.add_patch(Polygon(proj(lens)[0], closed=True, fc="0.35", ec="none", alpha=0.35, zorder=4))
draw_triangle(ax, Eabc, C0)
draw_triangle(ax, Eabd, C1)
label_vertices(ax)
for p, name, col in ((np.array([0, 1]), r"$\infty_0$", C0), (np.array([1, 0]), r"$\infty_1$", C1)):
    xy, d = proj(bloch(p.astype(complex)))
    if d[0] > -0.05:
        ax.plot(*xy[0], "x", ms=4, color=col)
        dy = 0.035 if xy[0][1] > 0 else -0.06
        ax.annotate(name, xy[0], xytext=(xy[0][0] + 0.025, xy[0][1] + dy), color=col, fontsize=9)
ax.set_title("(a) flat in each face's best chart", fontsize=9)

# (b) geodesic faces
ax = axes[1]
Gabc = edges(geo, [A, B, Cv])
Gabd = edges(geo, [A, B, D])
# the geodesic edge is the same curve, with the same parameter, from either end
assert max(fs(geo(A, B, l), geo(B, A, 1 - l)) for l in lam) < 1e-7
draw_triangle(ax, Gabc, CG)
draw_triangle(ax, Gabd, CG)
label_vertices(ax)
ax.set_title("(b) geodesic cones", fontsize=9)

for ax in axes:
    ax.set_xlim(-0.56, 0.56)
    ax.set_ylim(-0.6, 0.6)
fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01, wspace=0.05)
fig.savefig("faces.pdf")
dev = max(fs(flat(A, B, 0, l), flat(A, B, 1, l)) for l in lam)
print(f"wrote faces.pdf; |AB| = {fs(A, B):.3f} rad, max distance between the two flat versions of AB "
      f"at equal parameter = {dev:.3f} rad")
