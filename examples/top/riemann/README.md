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
2D Newton (`newton_on_triangle`), seeded first from the closed-form
solution of the *linear* interpolation of `F` at the 3 vertices (exact
for a linear field -- the affine analogue of the edge-crossing linear
interpolation used in marching tetrahedra), then from 7 fixed seeds
spread over the triangle, both to catch a second sheet the linear seed
misses and as a fallback when the linear seed is degenerate. Distinct
converged roots are cached per triangle.

The Newton iteration itself runs in a **local orthogonal frame**
(`tri_frame`/`build_tri_frame`), not raw barycentric coordinates: `x`
along the triangle's longest side, `y` along the perpendicular dropped
from the opposite vertex (its foot `Po` lands strictly inside the
longest side, since the angles at its two endpoints are the triangle's
two smallest). Raw barycentric edge vectors `(w1-w0,z1-z0)` and
`(w2-w0,z2-z0)` can be near-parallel for a needle-shaped triangle,
ill-conditioning the Jacobian through the parametrization alone,
independent of `F`; the longest-side frame's two basis directions are
genuinely orthogonal (in the real inner product on `C^2 = R^4`) and
scaled to the triangle's own extent, removing that source of
ill-conditioning. Seeds are converted barycentric-to-frame
(`bary_to_xy`) before iterating, and the result is converted back to
barycentric coordinates (w.r.t. the triangle's original vertex order)
immediately after convergence, so the domain/boundary check and
`compute_crossing_node`'s use of `face_op` indexing are unaffected --
the frame is purely an internal conditioning device. (This is a real,
not complex/holomorphic, change of coordinates: an earlier idea to
reparametrize by a single complex variable `x+iy` so each triangle's
crossings could be found by a closed-form/Durand-Kerner polynomial
solve turned out not to work in general -- restricting `F` to an
arbitrary real 2-plane in `C^2` is holomorphic in `x+iy` only when that
plane happens to be a complex line, which is not the generic case for
a mesh triangle.)

For a triangle with a vertex far enough in `w` or `z` to be near a pole,
`Ffun`/`Fwfun`/`Fzfun` (passed into `newton_on_triangle`) switch from
`F`/`Fw`/`Fz` to the corner-chart polynomial `F_corner`/`Fwr_corner`/
`Fzr_corner` and the triangle's vertices are expressed in `(wr,zr)`
instead of `(w,z)` -- Newton itself is unchanged, it just iterates in
whichever chart it's given. See *Known limitations* below for what this
chart is and its (real, curve-dependent) limits.

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
dimension 4) down to 3D via `project_for_viz`, in one of two modes.
Both are lossy projections, not embeddings -- each can and will show
self-intersections that are not really there.

- **flat** (default): `(Re w, Im w, Re z)`. Drops `Im z` entirely.
- **onion** (`--onion`): direction from `z`'s own position on the
  sphere (`to_sphere(z)`, the forward stereographic projection --
  inverse of `from_sphere`), displaced *perpendicular to that sphere*
  (i.e. radially) by `1 + scale*onion_radial(w)`. `onion_radial(w)` is
  a fixed-angle projection of `w`'s own stereographic image onto a
  generic direction, bounded to `[-1,1]` -- deliberately not the raw
  `Re(w)` axis, so it isn't blind to some future curve's symmetry other
  than `w -> -w` (every degree-2 catalog curve so far is `w^2=f_0(z)`,
  whose two roots are always a `+-w` pair; a fixed generic angle
  separates that case and generalizes better to others). For a fixed
  `z`, the (up to `n`) sheets of the surface spread into concentric
  shells instead of overlapping -- the branch structure becomes visible
  as nested spheres that pinch together where sheets meet. `--onion-scale
  X` (default `0.3`) sets the radial spread; output radius stays within
  `[1-X, 1+X]` by construction.

`project_for_viz` is deliberately isolated so other projections
(a local PCA-based projection, coloring by a dropped coordinate) are
easy to add as a third mode.

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
  at the seed vertex `z=0`. `--generic` breaks this alignment for the
  five *finite* landmark vertices via a fixed complex-affine change of
  coordinates applied to `F` (`GA w+GB`, `GC z+GD`), and independently
  breaks it for the *infinite* ones (`w=∞`, `z=∞`) via a fixed 3D
  rotation applied to the seed octahedron's own vertex placement
  instead (`g_rot_w`/`g_rot_z`, `sphere_rot::apply`) -- these are two
  different mechanisms because an affine map always fixes infinity, so
  it can never move the seed's own pole vertices off of it, while a
  sphere rotation (itself a unitary Möbius transform of `C_∞`, just
  applied to mesh construction rather than to `F`) can. This matters:
  by Riemann-Hurwitz a degree-`n` covering of `P^1` needs a fixed total
  ramification count, so a curve with an odd number of *finite* branch
  points (e.g. `parabola`, one at `z=0`) is forced to have one at
  `z=∞` too -- confirmed directly (`branch_gap_corner(zr)->0` as
  `zr->0`). Before the rotation fix, that seed vertex sat exactly on
  this real branch point *forever*, at every depth (Maubach bisection
  never touches an original seed vertex), giving a ~38-40% failure rate
  among corner-adjacent cells that stayed flat instead of shrinking
  with `--depth` -- unlike the ordinary alignment-sensitivity case
  above, which does shrink (roughly halving every +2 levels). With the
  rotation, no seed vertex sits at a literal pole any more, so this
  class of failure is gone (parabola, `--depth 10`: bad cells corner
  went from 174/190 to 0; overall bad_cells 177->16).
  `cell_priority` also now force-refines any cell whose 5 vertices
  straddle the `sphere_is_far` cutoff in `w` or `z` (`straddles_far`),
  on the theory that a cell mixing ordinary- and corner-chart triangles
  is inherently risky -- kept as a real improvement for configurations
  where it matters, though at the depths/thresholds tested here every
  such cell was already being pushed to `max_depth` by `branch_gap`
  anyway, so it measured as a no-op on top of the rotation fix.
- **Pole handling is projective at the joint corner.** `F` is monic in
  `w`, so `w=∞` is never a root for finite `z` (the `w`-homogenized
  polynomial is `1` at `w'=0`) -- the only place a second chart is
  needed is where `w` and `z` diverge *together*, handled via the
  doubly-homogenized polynomial `K(wr,zr) = wr^n zr^dz F_raw(1/wr,1/zr)`
  (`wr=1/(GA w+GB)`, `zr=1/(GC z+GD)`, finite at `(wr,zr)=(0,0)`),
  switched to per-triangle once *either* vertex crosses `sphere_is_far`
  (an ordinary-chart `|w|` or `|z|` past `1e6`). It used to require
  *both* `w` and `z` far before switching charts, on the assumption
  that "far in only one" never occurs for a monic-in-`w` curve -- true
  only in the literal `z=∞` limit: on `parabola`, `|w|~sqrt(|z|)`, so
  `z` crosses the cutoff long before `w` does, and that whole band was
  silently excluded (confirmed: the excluded-triangle count *grew*
  with `--depth` instead of shrinking). Fixed by switching to the
  corner chart whenever *either* is far; `from_wr`/`from_zr` also now
  clamp their recovered magnitude to the same `1e8` practical-infinity
  scale used elsewhere, since a root can otherwise land arbitrarily far
  out and break the plain-`F` residual diagnostic (not the solve
  itself, which stays well-conditioned in `(wr,zr)`).

  This closes the hole completely for a curve whose point at infinity is
  a *regular* point of `K=0` there (verified on `parabola`: `K(wr,zr) =
  zr-wr^2`, `∂K/∂zr=1≠0` at the origin -- the residual gap shrinks with
  depth exactly like anywhere else). It does *not* fully close the hole
  for `elliptic`: there, `K(wr,zr) = zr^3+(zr^2-1)wr^2`, and `K`, `∂K/∂wr`,
  `∂K/∂zr` all vanish at `(wr,zr)=(0,0)` -- a genuine cusp (`wr^2≈zr^3`
  near the origin), because `w` grows like `z^1.5` along this curve, a
  power incompatible with the `C_∞ x C_∞ = P^1 x P^1` grid. This is a
  structural fact about this curve in *this* ambient compactification
  (its standard smooth model lives in `P^2` instead, where infinity is a
  single regular point), not a bug -- more depth shrinks the residual gap
  somewhat but a cusp's degenerate tangent keeps it from vanishing the
  way it does at a regular point, and the `--generic` rotation fix above
  doesn't touch this either: a cusp is a coordinate-independent property
  of the curve, not an artifact of where the mesh's vertices sit.
  Resolving that fully would need a curve-specific local (Puiseux)
  reparametrization at the singular point, well beyond the general
  per-triangle treatment here.
- **Memory grows monotonically with depth.** `nmt`'s cell containers
  are append-only -- Maubach subdivision marks parent cells
  not-current rather than freeing them -- so deep runs can exhaust
  memory (observed: `--depth 16` needs several GB; see *Usage* for
  numbers). No compaction/GC of the mesh itself exists yet, but the
  per-vertex/per-triangle auxiliary caches were shrunk (dead `fw`/`fz`/
  `w_label`/`z_label` fields removed, the 2-root Newton cache moved out
  of `extra_data<2>` into a flat pool, matching the memory work done in
  `examples/top/riemann_cp2/riemann_cp2.cpp`) -- measured ~52% lower
  peak RSS at the same `--depth`/`--threshold`, same extraction result,
  with no other flags involved.
- **`branch_gap` is an order-of-magnitude proxy for `n>2`.** It's the
  `n(n-1)`-th root of the resultant, which is exact for a single pair
  of roots (`n=2`) but is a product over *all* pairs for `n>2`, not
  just the closest pair. Untested on cubic-or-higher `F`.
- **`--threshold` is relative, not absolute, so far-from-any-branch-point
  cells can stay large in raw `(w,z)` terms even after "converging".**
  `cell_priority`/`cell_priority_bernstein` both compare a cell's own
  size against something LOCAL (`branch_gap`, or the Bernstein box
  width) -- a region with nothing driving refinement locally stops
  early, at whatever absolute size it happens to be at. Symptom seen and
  fixed once: giant spikes in a `--onion` render, traced to two
  genuinely different, individually-correct roots of `F`, found on
  different faces of the same (large, but "converged" per the ratio)
  tetrahedron, connected as if adjacent. Fixed a real, unconditional bug
  in `cell_priority` (it only ever measured `w`'s own diameter, never
  `z`'s -- now uses the combined `cell_diam`, same as `--proximity`) and
  in the 4-node pairing heuristic (now measured on the `w`/`z` spheres,
  not raw chart distance, the same lesson
  `examples/top/riemann_cp2/riemann_cp2.cpp`'s analogous pairing already
  applied via `fs_dist`). This measurably helps under the *default*
  `cell_priority`, and a stricter `--threshold` (e.g. `0.01` vs the
  default `0.1`) helps further there too (confirmed: worst edge length
  dropped roughly 10x at the same `--depth`). **Under `--bernstein`,
  though, its own priority formula has the same "relative, not
  absolute" blind spot in a different shape, and does NOT respond the
  same way to a stricter `--threshold`** (tested: made a case measurably
  *worse*, not better) -- this remains unresolved; `--bernstein`'s
  own size-insensitivity would need its own dedicated fix, not a
  borrowed one.

## Usage

Build (no project file yet; plain g++ against the VGTL headers):

```sh
g++ -I ../../../include -std=c++17 -O2 riemann.cpp -o riemann
```

(the `--generic`/`--generic-seed` rotation already needed `<random>`'s
`std::mt19937`, i.e. at least C++11, before this; `-std=c++03` as
originally documented here no longer builds at all on a current
toolchain.)

(from `examples/top/riemann/`; adjust the include path if run from
elsewhere -- it must point at VGTL's top-level `include/` directory).

Run:

```sh
./riemann [--function N] [--generic] [--generic-seed N] [--depth N] [--threshold X] [--onion] [--onion-scale X] [--cutoff X] [--proximity] [--bernstein] [--bernstein-level N] [--bernstein-selftest] [--list-functions]
```

| flag | default | meaning |
|---|---|---|
| `--function N` | `0` | index into the curve catalog (see below) |
| `--generic` | off | apply a random rotation to the seed octahedron's own vertex placement (independently for `w` and `z`) to break curve/mesh symmetry alignment -- see *Known limitations* |
| `--generic-seed N` | `12345` | implies `--generic`; picks which random rotation |
| `--depth N` | `12` | maximum Maubach subdivision level |
| `--threshold X` | `0.1` | refinement stops popping cells once the priority queue's max drops below this |
| `--onion` | off | export via the onion projection (nested spheres) instead of the flat `(Re w, Im w, Re z)` one -- see Phase 4 |
| `--onion-scale X` | `0.3` | radial spread for `--onion`; output radius stays within `[1-X, 1+X]` |
| `--cutoff X` | off | drop any extracted polygon with a vertex farther than `X` from the origin in the projected 3D point -- keeps pole-proxy outliers from blowing out the flat projection's bounding box; ignored under `--onion`, where the radius is already `~1` by construction |
| `--proximity` | off | weight `cell_priority` by `diam/(mind+diam)` (`mind` = a first-order distance-to-curve estimate at the cell's own vertices) to bias refinement toward the curve itself -- ported from `riemann_cp2.cpp`; cuts subdivisions/final cell count without changing the bad-cell rate |
| `--bernstein` | off | use a certified Bernstein-Bezier enclosure of `F` over each 2-simplex, instead of `branch_gap`, both to PRUNE cells that provably can't contain a crossing (dropped forever, never refined again) and to rank the survivors -- ported from `riemann_cp2.cpp` |
| `--bernstein-level N` | `0` | implies `--bernstein`; rounds of 1-to-4 triangular subdivision used to tighten the certified enclosure before caching it |
| `--bernstein-selftest` | -- | brute-force-verify the Bernstein enclosure machinery against a sample triangle (independent of the mesh/curve catalog) and exit |
| `--list-functions` | -- | print the curve catalog and exit |

Output is `riemann_surface_<name>[_generic][_onion].obj` in the working
directory (e.g. `riemann_surface_elliptic.obj`,
`riemann_surface_parabola_generic_onion.obj`), viewable in any standard
OBJ viewer (Blender, MeshLab, a three.js-based page, ...). Progress and
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
