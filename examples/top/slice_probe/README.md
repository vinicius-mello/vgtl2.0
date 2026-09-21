# slice_probe

A controlled stress test for the codim-2 crossing-graph extraction
algorithm shared by `examples/top/riemann/riemann.cpp` and
`examples/top/riemann_cp2/riemann_cp2.cpp` (the "Phase 3" logic: find
where a codim-2 locus crosses each triangle of a 4-simplex mesh,
canonically identify each crossing point, and check that a cell's
crossings close into a polygon), isolated from everything algebraic
(Newton iteration, poles/charts, branch points).

## Setup

Take the standard 4-simplex in `R^4` (vertices `0, e1, e2, e3, e4`).
Draw a uniformly random interior point `P` and two random unit
directions `u, v`. The affine subspace

```
L = { x : dot(x-P,u)=0, dot(x-P,v)=0 }
```

is codimension 2 and passes through `P`, so it always meets the
simplex. Maubach-refine the seed 4-simplex (only cells that could
possibly contain a crossing -- an EXACT filter here, since the field
is literally affine and an affine map's extrema over a simplex are
always at its vertices), then run the same per-tetrahedron/per-cell
crossing-graph closure `riemann.cpp` uses, with the per-triangle
crossing found by an exact 2x2 linear solve instead of Newton (the
field really is affine on a triangle, so the linear solve isn't just
a seed -- it's the answer).

## Why this is a strong test, not just a simpler one

Every 4-simplex cell is convex, so its intersection with `L` is
provably empty, a point, or a *single* convex polygon -- never more
than one polygon, never an open (non-closing) fan. So here, unlike the
algebraic-curve case, a "bad cell" (crossing graph doesn't close) or a
cell producing more than one cycle (`multi_cycle_cells`) cannot be a
real geometric outcome: it can only be a bug or a tolerance/precision
issue in the extraction code itself.

## Usage

```sh
g++ -I ../../../include -std=c++17 -O2 slice_probe.cpp -o slice_probe
./slice_probe [--trials N] [--depth N] [--seed N] [--dump-dir DIR] [--dump-budget N]
              [--corner-bias X] [--parallel-eps RAD] [--axis-align-p X]
```

| flag | default | meaning |
|---|---|---|
| `--trials N` | 2000 | number of independent random `(P,u,v)` draws |
| `--depth N` | 8 | max Maubach refinement level (uniform within the relevant region -- no threshold needed, since relevance is exact here) |
| `--seed N` | 12345 | base RNG seed; trial `i` uses a distinct derived seed, so a run is fully reproducible |
| `--dump-dir DIR` | off | for each bad cell (up to `--dump-budget`), write a `.txt` (P, u, v, node/adjacency dump, all in full precision) and a `.obj` wireframe (projected into an orthonormal frame adapted to `L`, so a genuinely-on-`L` node has zero residual on 2 of its 4 raw coordinates -- a real gap between two edges meant to meet at "the same" node is visible directly in Blender) |
| `--corner-bias X` | 0 | pull `P` this fraction of the way toward a random *original* simplex vertex -- see *Finding* below |
| `--parallel-eps RAD` | 0 | force the angle between `u` and `v` to exactly this many radians instead of a fully random one |
| `--axis-align-p X` | 0 | probability that `u` (independently, `v`) is replaced by a coordinate axis or simplex edge direction instead of a generic random unit vector |

## Finding (2026-09-21 run)

Pure random `(P,u,v)`: **zero bad cells and zero multi-cycle cells**
over ~700k extracted cells, depths 8 through 20 (`--trials 3000
--depth 10`, `--trials 100 --depth 20`, etc.) -- confirms the
crossing-graph algorithm itself is sound in the generic case, at
double-precision depths well past what's practical for the real
curves.

`--corner-bias 0.999999 --depth 16`: **62 bad cells / 94270 extracted
(0.066%), in 44/500 trials (8.8%)** -- reproducible with
`--seed`. Every dumped failure has the exact same shape: exactly 2
crossing nodes, both within ~1e-6 of the *same* physical point (an
**original, never-refined seed vertex** that `P` happens to land
extremely close to -- the direct codim-2 analogue of `riemann.cpp`'s
documented "branch point sitting exactly on a seed vertex forever"
bug), but classified inconsistently: e.g. `fail_t48_c76.txt` finds the
crossing on the *same edge* (`desc=166`) from two different incident
triangles, at parameters `6.49e-6` and `7.85e-6` -- both genuinely
tiny, but differing by `1.36e-6`, just over `same_crossing_node`'s
fixed `1e-6` merge tolerance (shared verbatim by `riemann.cpp`), so
two spurious nodes are created for what is really one crossing on one
edge. Root cause: each triangle's crossing is solved independently, so
two triangles sharing an edge can (and near a preserved low-order mesh
feature, will) compute a tiny but nonzero *disagreement* about exactly
where on that edge the crossing sits -- and `same_crossing_node`'s
tolerance is a fixed absolute number, not scaled to the actual
floating-point precision available at that depth/position. This is the
same tolerance, doing the same job, in the *exact* code both
`riemann.cpp` and `riemann_cp2.cpp` use -- so the fix (a looser and/or
depth-adaptive merge tolerance, or computing each shared edge's
crossing once and reusing it from both incident triangles instead of
re-deriving it twice) should transfer directly.

`--axis-align-p 1.0 --depth 16`: zero bad cells, but 12.2M degenerate
(rank-deficient) triangles out of ~480k actual crossings -- the
algorithm correctly *skips* triangles it can't solve rather than
misbehaving, but this leaves open whether skipped-degenerate triangles
ever silently drop a real crossing without tripping `bad_cells` (not
yet instrumented here -- `g_tri_in_plane` covers only the
whole-triangle-in-`L` extreme, not a triangle that's rank-deficient
but still genuinely near a crossing).

`--parallel-eps 0.0001 --depth 16`: zero bad cells -- an
ill-conditioned *definition* of `L` itself (nearly parallel `u,v`)
does not, on its own, reproduce the failure; it's specifically about a
crossing landing near a preserved mesh feature.

## Root cause, and a fix that works (2026-09-21)

Added a per-tetrahedron trace to the failure dump (which raw
`(l0,l1,l2)` each of a cell's 10 triangle-faces produced, and which
node index it got matched to) to see *why* a tetrahedron ends up with
the wrong node count. For `fail_t22_c58`, it shows the real mechanism:
two genuinely DIFFERENT triangles (desc `347` and `377`), each with its
own distinct crossing point (differing by `~4e-7` -- far above
floating-point noise, and the same order as how far each point
actually sits off the shared edge, so two real, separate near-edge
points, not one point measured twice), both land within the *same*
`1e-6` window of shared edge `147` and get **merged** by
`same_crossing_node`. That collapses what should be a proper 3-node
triangle cycle (`node0 -> A -> B -> node0`) into a spurious `node0 <->
A` digon: the tetrahedron that should contribute the `(A,B)` edge
instead sees `A` and `B` as the same node and reports a single
"touching" point, no edge.

This is the opposite of what a *looser* tolerance would fix (confirmed
by the failed experiment above) -- the actual bug is that one constant,
`1e-6`, was doing two different jobs: (1) **classification** --
"is this triangle's own crossing close enough to a face to treat it as
on that face?", and (2) **identity** -- "are these two
already-classified points actually the same point?". A genuinely
identical point found via two different incident triangles agrees far
more tightly than `1e-6` (when a point is truly *on* a shared edge,
both triangles' formulas reduce to the same edge-endpoint-only
computation, so the only disagreement left is ordinary floating-point
roundoff, `~1e-10`-ish here) -- so the identity check needs to be much
stricter than the classification check, not looser.

Fix: split the one constant into `SNAP_TOL=1e-6` (classification,
unchanged) and `MERGE_TOL=1e-9` (identity, three orders of magnitude
tighter). Result on the same `--corner-bias 0.999999 --depth 16 --seed
12345` sweep:

| | bad cells | trials affected | touching_tets |
|---|---|---|---|
| baseline (`1e-6` for both) | 62 / 94270 (0.066%) | 44/500 (8.8%) | 295 |
| **`SNAP_TOL=1e-6`, `MERGE_TOL=1e-9`** | **4 / 94309 (0.004%)** | **2/500 (0.4%)** | 3 |

A 94% reduction, with **zero regressions** on the plain-random battery
(re-verified: 0 bad cells over ~600k cells, depths 16 and 20). The 4
remaining failures are structurally different and more involved (one
has 5 candidate nodes and a tetrahedron with 3 distinct crossings,
`unhandled_tets=1`) -- plausibly a genuine lower-dimensional
degeneracy (`L` tangent to a sub-face of that specific cell) rather
than the same tolerance-conflation bug; not yet root-caused.

**This fix (splitting `same_crossing_node`'s bare `1e-6` into a
classification tolerance and a much tighter identity tolerance) should
transfer directly to `riemann.cpp` and `riemann_cp2.cpp`, which share
the exact same one-constant-for-two-jobs code.** Not yet ported there.
