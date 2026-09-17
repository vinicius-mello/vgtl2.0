# Riemann surface triangulation

An adaptive simplicial algorithm that triangulates the Riemann surface
of an algebraic function `w = w(z)`, given as the zero locus of a
polynomial

```
F(w,z) = w^n + f_{n-1}(z) w^{n-1} + ... + f_1(z) w + f_0(z)
```

monic in `w`, with `z, w ∈ C_∞ = C ∪ {∞}` and each `f_i` a polynomial
in `z`. The output is a piecewise-linear approximation (a triangle
mesh) of the surface `F(w,z) = 0`.

This is a VGTL example built on the library's `nmt<4>` pseudo-manifold
container and its Maubach (bisection) subdivision scheme. Everything
specific to this problem lives in two files:

- `functions.hpp` -- the catalog of test curves (`poly_F`, `cx`, and
  `function_catalog()`). Add a new curve here without touching the
  algorithm.
- `riemann.cpp` -- the algorithm itself: seed mesh, adaptive
  refinement, surface extraction, and OBJ export.

## Why this is hard: the dimension count

`C_∞^2` has real dimension 4. For a generic algebraic curve, `F = 0`
has real codimension 2 in it (`Re F` and `Im F` are two independent
real conditions), so the zero locus is a real surface (dimension 2)
sitting inside a 4-dimensional simplicial mesh. This is a genuine
codimension-2 isosurface extraction, not the classical codimension-1
case that marching-cubes-style algorithms handle.

The general rule for expected intersection dimension is

```
dim(probe) + dim(zero locus) - dim(ambient) = dim(intersection)
```

- An **edge** (dim 1) against `F=0` (dim 2) in the dim-4 ambient space
  gives `1+2-4 = -1`: generically, edges do **not** cross the surface,
  even though `Re F` alone (a weaker, codimension-1 condition) may
  change sign along them.
- A **triangle** (dim 2) gives `2+2-4 = 0`: isolated points. This is
  the right place to look for crossings -- solving `Re F = Im F = 0`
  for the triangle's two barycentric coordinates -- exactly analogous
  to the role an *edge* plays in classical marching tetrahedra (where
  the field is real-valued and codimension 1: `1+2-3=0`).

Every design choice below follows from this: bisection happens on
edges (to refine the mesh), but crossing detection happens on
triangles (to extract the surface).

## Method

### Phase 1 -- seed mesh

`C_∞` is triangulated as the standard octahedron via stereographic
projection: vertices `{0, ∞, 1, -1, i, -i}` (`0` = south pole,
`∞` = north pole), 8 triangles, with orientation signs

```
+<0,1,i>  -<0,-1,i>  +<0,-1,-i>  -<0,1,-i>      (south)
-<∞,1,i>  +<∞,-1,i>  -<∞,-1,-i>  +<∞,1,-i>      (north)
```

(the two hemispheres must carry opposite sign patterns around a
shared equatorial edge for the octahedron itself to be orientation
coherent -- verified by `count_incoherent()`).

`C_∞^2` is seeded as the simplicial product of two copies of this
octahedron: for each pair of triangles `(σ_w, σ_z)` and each of the 6
monotone paths across the 3x3 grid of their vertices, one 4-simplex is
produced (the classical prism-to-simplices decomposition). The
orientation of each resulting 4-simplex is the Eilenberg-Zilber
shuffle sign `ε_w · ε_z · sign(path)`.

This gives 36 vertices and 384 4-simplices. The code cross-checks the
hand-derived shuffle formula against an independent BFS re-derivation
of a coherent orientation over the cell-adjacency graph (both must
agree, and both are checked against `euler_characteristic(t) == 4`,
matching `S^2 x S^2`, and zero boundary facets).

### Phase 2 -- adaptive refinement

Refinement is Maubach bisection (`maubach_subdivide`) driven by a
priority queue, biased toward the branch locus of the curve.

- **Edge bisection rule**: the new vertex's `(w,z)` is the
  *spherical* midpoint of the edge's endpoints -- each of `w` and `z`
  is projected to its point on the unit sphere (via stereographic
  projection), the two unit vectors are summed and renormalized (this
  is exact SLERP at `t=1/2`), then mapped back to `C_∞`. Degenerates
  only for antipodal endpoints (`{0,∞}`), which never occurs in the
  seed mesh and is guarded with a fallback if it ever arose from
  refinement.
- **Refinement criterion**: each vertex caches `branch_gap(z)`, an
  estimate of the smallest pairwise gap between the `n` roots of
  `F(., z)` in `w`, computed from the resultant
  `Res(F(.,z), F_w(.,z))` via a Sylvester-matrix determinant (this
  *is* (up to sign and a normalizing root) the discriminant of `F` in
  `w`; it vanishes exactly at branch points and shrinks near them).
  A cell's priority is `wdiam / min(branch_gap over its vertices)`:
  large when the cell is wide in `w` relative to how close the sheets
  get there -- the actual risk factor for one of its triangles
  containing more than one root, not just "some vertex is near a
  branch point" (which alone says nothing about the cell's size).
  Cells pop off a max-priority queue and get bisected until the queue
  drains below `--threshold` or hits `--depth`.
- Vertices whose `w` (or `z`) sits near the north-pole proxy (see
  *Known limitations*) are excluded from refinement priority.

### Phase 3 -- surface extraction

Per 2-simplex (`triangle_intersection`): find where `F=0` crosses the
triangle, as up to **two** points (`F` has degree `n` in `w`, so a
triangle's `w`-span can still straddle more than one root even after
refinement -- `branch_gap` bounds how close two roots get, not
whether a given triangle spans two of them). Each root is sought by
2D Newton on the barycentric coordinates, seeded first from the
closed-form solution of the *linear* interpolation of `F` at the 3
vertices (exact for a linear field -- the affine analogue of the
edge-crossing linear interpolation used in marching tetrahedra), then
from 7 fixed seeds spread over the triangle, both to catch a second
sheet the linear seed misses and as a fallback when the linear seed
is degenerate. Distinct converged roots are cached per triangle.

Each crossing point is given a **canonical identity**
(`compute_crossing_node`) based on *where* it sits on the ambient
4-simplex, not on which triangle happened to find it: a mesh vertex,
a point on a mesh edge (with a scale-invariant parameter along it), or
a triangle interior point (with a root index, so two sheets through
the same triangle stay distinct). This is what makes the algorithm
robust to a crossing landing exactly on a lower-dimensional face
shared by several triangles -- a case that is not rare (see *Known
limitations*) and, if identified only by raw `(w,z)` distance, breaks
under refinement because a genuinely short segment and a true
duplicate shrink at the same rate.

Per 4-simplex `σ`: each of its 5 tetrahedral facets contributes an
edge between the (at most 2, rarely 4) crossing nodes found on its 4
triangular faces. Since each triangle of `σ` belongs to exactly 2 of
its 5 tetrahedra, these edges close into one or more cycles (checked
by walking the node-adjacency graph, not assumed) of length 3 to 5.
A `σ` with more than one cycle means more than one sheet threads it;
each cycle becomes one output polygon.

### Phase 4 -- projection and export

The extracted polygons are projected from `(w,z) ∈ C_∞^2` (real
dimension 4) down to 3D via `project_for_viz`: currently
`(Re w, Im w, Re z)`. This drops `Im z` entirely, so it is a
*projection*, not an embedding -- it will show self-intersections
that are not really there. The function is deliberately isolated so
other projections (stereographic re-embedding, a local PCA-based
projection, coloring by the dropped coordinate) are easy to try.

`obj_writer` writes one OBJ face per extracted polygon, with its own
unwelded vertices (no attempt to merge vertices shared with polygons
from adjacent 4-simplices -- see *Known limitations*).

## Known limitations / open threads

- **Output mesh is not welded.** Each polygon owns its vertices; nothing
  merges vertices shared between neighboring 4-simplices. Fine for a
  first visual check, not a topologically clean mesh.
- **No orientation propagation to the output.** Extracted polygons
  don't yet carry a consistent winding induced from `σ`'s orientation
  (or from the curve's intrinsic complex-tangent orientation
  `(v, iv)`, which is available since `F=0` is a genuine complex
  curve and this would double as a correctness check).
- **Coordinate-alignment sensitivity.** The dominant remaining failure
  mode (`bad_cells`, cells whose crossing graph doesn't close into
  clean cycles) is far more common when the curve's own symmetry lines
  up with the seed octahedron's (real coefficients, seed vertices on
  the real/imaginary axes) -- e.g. `w^2-z` has a branch point exactly
  at the seed vertex `z=0`. `--generic` breaks this alignment via a
  fixed complex-affine change of coordinates and consistently cuts the
  bad-cell rate by roughly half to an order of magnitude on the curves
  tested so far, without changing the curve. The remainder still
  shrinks with `--depth` (roughly halving every +2 levels in the runs
  so far) and concentrates one level below the cap, consistent with a
  sampling limit rather than a structural one -- not yet fully
  isolated.
- **Pole handling is a placeholder.** `w=∞` or `z=∞` are approximated
  by a large finite value (`1e8`) rather than treated projectively;
  cells/triangles that touch this proxy are excluded from refinement
  and extraction (`near_pole`). Fine as long as the curve of interest
  doesn't pass near a pole; not correct in general.
- **Memory grows monotonically with depth.** `nmt`'s cell containers
  are append-only -- Maubach subdivision marks parent cells
  not-current rather than freeing them -- so deep runs can exhaust
  memory (observed: `--depth 16` needs several GB; see *Usage* for
  numbers). No compaction/GC exists yet.
- **`branch_gap` is an order-of-magnitude proxy for `n>2`.** It's the
  `n(n-1)`-th root of the resultant, which is exact for a single pair
  of roots (`n=2`) but is a product over *all* pairs for `n>2`, not
  just the closest pair. Untested on cubic-or-higher `F`.

## Usage

Build (no project file yet; plain g++ against the VGTL headers):

```sh
g++ -I ../../../include -std=c++03 -O2 riemann.cpp -o riemann
```

(from `examples/top/riemann/`; adjust the include path if run from
elsewhere -- it must point at VGTL's top-level `include/` directory).

Run:

```sh
./riemann [--function N] [--generic] [--depth N] [--threshold X] [--list-functions]
```

| flag | default | meaning |
|---|---|---|
| `--function N` | `0` | index into the curve catalog (see below) |
| `--generic` | off | apply a fixed complex-affine change of coordinates to `F` (rotation + translation, `\|a\|=\|c\|=1`) to break curve/mesh symmetry alignment -- see *Known limitations* |
| `--depth N` | `12` | maximum Maubach subdivision level |
| `--threshold X` | `0.1` | refinement stops popping cells once the priority queue's max drops below this |
| `--list-functions` | -- | print the curve catalog and exit |

Output is `riemann_surface_<name>[_generic].obj` in the working
directory (e.g. `riemann_surface_elliptic.obj`,
`riemann_surface_parabola_generic.obj`), viewable in any standard OBJ
viewer (Blender, MeshLab, a three.js-based page, ...). Progress and
diagnostics (mesh invariants, refinement stats, extraction quality,
residuals) print to stdout -- worth reading; in particular
`extracted cells: ok=.. bad=..` and the `ok/bad cells by cell level`
breakdown are the headline numbers for judging a run.

### The curve catalog

```sh
./riemann --list-functions
```

```
  0: elliptic -- w^2 - z^3 + z : smooth elliptic curve (genus 1); branch points at z=0,+1,-1.
  1: parabola -- w^2 - z : smooth parabola; one branch point at (0,0), which lands exactly on the octahedron seed vertex z=0.
```

To add a curve, append an entry to `function_catalog()` in
`functions.hpp` -- a name, a one-line description (branch points,
degree, anything relevant to interpreting a run), and the `poly_F`
coefficient table (`c[i]` = coefficients of `f_i(z)`, lowest degree
first, `i=0..n-1`). Nothing in `riemann.cpp` is specific to degree 2;
higher-degree curves are untested but should work mechanically (see
the `branch_gap` caveat above).

### Reference numbers (elliptic curve `w^2=z^3-z`, `--generic`)

Depth-vs-quality trend from recent runs (`--threshold 0.1`), useful
for calibrating how deep to go:

| `--depth` | 4-simplices | wall time | cells closed |
|---|---|---|---|
| 12 | 1.35M | ~26s | 99.0% |
| 14 | 4.6M | ~89s | 99.4% |
| 16 | 16.2M | ~6.3min | 99.7% |

Memory is the binding constraint at depth 16 (observed OOM-kill on a
7GB-RAM machine at the first attempt; succeeded once other processes
were closed). Wall time and cell count both grow roughly ×3.4 per
+2 levels in this range.
