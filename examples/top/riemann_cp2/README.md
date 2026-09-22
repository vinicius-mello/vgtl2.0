# Riemann surface triangulation, CP^2 approach

A parallel implementation of `examples/top/riemann/riemann.cpp`'s
algorithm (adaptive simplicial extraction of the zero locus of an
algebraic curve), working directly in the complex projective plane
`CP^2` instead of the product of two Riemann spheres `C_infty^2`. A
smooth plane curve `F(X,Y,Z)=0` (homogeneous, degree `n`) sits in
`CP^2` with real codimension 2, the same "codimension-2 isosurface
extraction" problem `riemann.cpp`'s own README explains in detail --
read that first for the shared rationale (why bisection happens on
edges but crossing detection on triangles, the per-tetrahedron
crossing-graph closure, canonical crossing-node identity, etc.). This
file's own comments focus on what's different here.

- `functions_cp2.hpp` -- the catalog of test curves, as homogeneous
  polynomials in `(X,Y,Z)`.
- `riemann_cp2.cpp` -- the algorithm itself.

## What's different from `riemann.cpp`

- **Seed mesh**: Gaifullin's 15-vertex/108-cell triangulation of
  `CP^2` (`gaifullin_points()`/`gaifullin_cells()`), already
  Maubach-compatible by construction (a backtracking 5-edge-coloring
  search over the seed's dual graph, verified offline against GAP's
  simpcomp) -- no barycentric subdivision step needed, unlike the
  project's earlier Kuhnel-seed attempt (`Triangulação PC2.txt` in
  this directory is that earlier, superseded design note).
- **Charts**: 3 symmetric affine charts (`X!=0`, `Y!=0`, `Z!=0`)
  instead of an ordinary/corner pair -- every `CP^2` point has at
  least one homogeneous coordinate with `|coordinate| >= 1/sqrt(3)`
  (unit-norm representative), so a best-conditioned chart always
  exists per triangle (`pick_chart`), no special-casing needed.
- **`--generic`/`--generic-seed`** rotates the Gaifullin seed's own 15
  vertex positions by a random unitary matrix (`random_unitary` +
  `apply_unitary`), the direct `CP^2` analogue of `riemann.cpp`'s
  sphere rotation -- most catalog curves hit several seed vertices
  exactly and need this (see the catalog table below).
- **Newton seeds and the Bernstein fallback**: identical in spirit to
  `riemann.cpp` (see that file's own comment on
  `bernstein_locate_root`) -- 1 dynamic (linear-interpolation) seed +
  4 fixed seeds (centroid, 3 near-vertex), then a certified
  Bernstein-Bezier bisection fallback tried only when all 5 miss. The
  3 near-edge-midpoint seeds this used to also try were dropped after
  measurement showed they never found a genuinely new root on any
  curve/flag combination tested there; the fallback covers the rare
  cases they occasionally (and only by luck) used to catch.
- **Crossing-node identity**: same `CROSSING_SNAP_TOL`/
  `CROSSING_MERGE_TOL` split as `riemann.cpp` (see that file's own
  comment for the root-caused bug this fixes -- two different
  triangles' genuinely distinct near-edge crossings getting wrongly
  merged by a single shared tolerance).
- **Projections**: `flat` (per-**polygon**, not per-point, best-
  conditioned chart -- fixes a real stitching bug, see
  `pick_chart_polygon`'s own comment), `--onion` (a `CP^2` analogue of
  `riemann.cpp`'s onion mode, direction from a *fixed* `Z/X` ratio on
  `S^2`, went through 3 rounds of real bugs -- see `project_for_viz`'s
  comment for the full history), and `--alpha-projection` (see below).

## `--alpha-projection`

`alpha-tilde([X:Y:Z]) = (|X|^2/S, Re(conj(X)Y)/S, Im(conj(X)Y)/S)`,
`S=|X|^2+|Y|^2+|Z|^2` -- from Seth Dutter, "Visualization of Complex
Projective Curves" (arXiv:2608.04323), attributing the map itself to
S. Kranich (2015). Unlike `flat`/`--onion`, this is a function of the
homogeneous point *itself*, not an affine chart of it -- no chart
choice, no stitching-bug class to have (`pt3` is already exactly this
map's domain, `CP^2`). Bounded by construction: the image is the
closed ball of radius `1/2` centered at `(1/2,0,0)` (confirmed
empirically: the max distance from that center over a real extracted
mesh was exactly `0.5`) -- no `--cutoff` needed, same as `--onion`.
Not injective -- `alpha(u,v)=alpha(u',v')` exactly when
`(u',v')=lambda*(u,v)` for a unit-**modulus** `lambda` (an overall
phase only, not a general rescaling), strictly less information
thrown away than `flat`'s (which drops an entire real coordinate,
`Im(b)`, outright). `--onion` and `--alpha-projection` are mutually
exclusive (checked, rejected with an error if both are passed).

## Usage

Build:

```sh
g++ -I ../../../include -std=c++17 -O2 riemann_cp2.cpp -o riemann_cp2
```

Run:

```sh
./riemann_cp2 [--function N] [--depth N] [--threshold X] [--cutoff X] [--onion] [--onion-scale X] [--alpha-projection] [--proximity] [--generic] [--generic-seed N] [--bernstein] [--bernstein-level N] [--bernstein-selftest] [--list-functions]
```

| flag | default | meaning |
|---|---|---|
| `--function N` | `1` (fermat_cubic) | index into the curve catalog (see below) |
| `--depth N` | `10` | maximum Maubach subdivision level |
| `--threshold X` | `0.1` | refinement stops popping cells once the priority queue's max drops below this -- the real lever for mesh density, see *Memory* below |
| `--cutoff X` | off | drop any extracted polygon with a vertex farther than `X` from the origin (`flat` mode only -- `--onion`/`--alpha-projection` are already bounded by construction) |
| `--onion` | off | export via the onion projection (nested shells) instead of `flat` |
| `--onion-scale X` | `0.5` | radial spread for `--onion` |
| `--alpha-projection` | off | export via Dutter's `alpha-tilde` map (see above); mutually exclusive with `--onion` |
| `--proximity` | off | bias `cell_priority` toward the curve itself, ported from `riemann.cpp`'s own `--proximity` |
| `--generic` | off | rotate the Gaifullin seed's own 15 vertices by a random unitary matrix, breaking curve/seed alignment |
| `--generic-seed N` | `12345` | implies `--generic`; picks which random rotation |
| `--bernstein` | off | use the certified Bernstein-Bezier enclosure to prune whole cells during refinement, not just as the per-triangle root-finding fallback (which runs regardless of this flag) |
| `--bernstein-level N` | `0` | implies `--bernstein`; rounds of 1-to-4 subdivision tightening the certified enclosure |
| `--bernstein-selftest` | -- | brute-force-verify the Bernstein enclosure machinery and exit |
| `--list-functions` | -- | print the curve catalog and exit |

Output is `riemann_cp2_<name>[_onion][_alpha].obj`. Progress and
diagnostics print to stdout -- `extracted cells: ok=.. bad=..` and the
`Bernstein fallback seeds used: ..` line are the headline numbers for
judging a run.

### The curve catalog

```sh
./riemann_cp2 --list-functions
```

```
  0: conic -- XY-Z^2=0 : smooth conic (genus 0); hits 8/15 Gaifullin seed vertices exactly -- needs --generic-seed
  1: fermat_cubic -- X^3+Y^3+Z^3=0 : smooth elliptic curve (genus 1); hits 0/15 exactly, but still needs --generic-seed (some resonance beyond exact vertex incidence, not yet understood)
  2: cubic_generic -- X^3+Y^3+Z^3-3XYZ=0 : Hesse pencil member away from Fermat; hits 3/15 exactly -- needs --generic-seed
  3: quartic -- X^4+Y^4+Z^4=0 : smooth quartic (genus 3); hits 8/15 exactly -- needs --generic-seed
  4: parabola -- X^2-YZ=0 : homogenization of riemann.cpp's w^2-z; one branch point at (w,z)=(0,0); hits 8/15 exactly -- needs --generic-seed
  5: elliptic -- X^2*Z-Y^3+Y*Z^2=0 : homogenization of riemann.cpp's w^2-z^3+z; genus 1, branch points at z=0,+1,-1; hits 2/15 exactly -- needs --generic-seed
```

To add a curve, append an entry to `function_catalog_cp2()` in
`functions_cp2.hpp`.

### Memory

`nmt`'s cell containers are append-only (Maubach subdivision marks
parent cells not-current rather than freeing them), so deep/low-
threshold runs grow memory monotonically with no ceiling of their own
-- see `riemann.cpp`'s README for the same caveat. Measured directly
on a 7.2GB machine (`--generic-seed 1`, both `elliptic` and
`parabola`): `--threshold 0.02 --depth 18` completed (peak RSS
~3.9-4.3GB, 0.84%/0.95% bad-cell rate); one step further in either
direction (`--threshold 0.01 --depth 19`, `--threshold 0.005 --depth
20`) exhausted memory before refinement even finished, on both curves
tried, twice each. Bad-cell rate dropped consistently and
substantially as threshold decreased (2.6% at the `0.1` default down
to 0.84%) -- lowering `--threshold` is the effective lever for mesh
quality here, not raising `--depth` alone: refinement's own priority
queue typically drains well before `--depth`'s cap binds (confirmed:
depth 18 vs 20 gave byte-identical output at `--threshold 0.1`, since
the queue had already emptied by depth ~18-19).
