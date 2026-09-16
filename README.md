# VGTL

A C++ template library for simplicial topology and geometry — generic
simplicial complexes, incidence/face-operator queries, subdivision schemes
(including Maubach's bisection), isosurface extraction, and a small
linear-algebra layer on top of GSL — written between 2008 and 2010 as part
of a PhD thesis:

> **Novos Métodos Simpliciais em Computação Gráfica**
> (*New Simplicial Methods in Computer Graphics*)
> Vinícius Moreira Mello, April 2006

The thesis itself is included as [`thesis-vinicius.pdf`](thesis-vinicius.pdf)
(in Portuguese).

## Layout

- **`include/vgtl/`** — the library itself (header-only, except for the
  Fortran routine below): `alg/` (vectors, points, quaternions, barycentric
  algebra, GSL wrappers), `comb/` (combinatorics: permutations, combinations,
  binomials), `geo/` (bounding boxes, distances, mesh I/O), `top/` (the
  simplicial-complex core: face operators, incidence, Euler characteristic,
  subdivision schemes, and the `model/` subdirectory with several concrete
  complex representations — `vc`, `lc`, `nmt`, ...), `opt/` (an L-BFGS-B
  wrapper and simplex-domain optimization), `cg/` (arcball/trackball camera
  controls, thin OpenGL wrappers).
- **`src/`** — `routines01.f`, the Fortran L-BFGS-B implementation the `opt/`
  module wraps; builds into `lib/libvgtl.a`.
- **`examples/`** — small, focused programs exercising each part of the
  library (`alg/`, `opt/`, `top/`).
- **`apps/isoview/`** — a Qt5/OpenGL viewer for interactively subdividing a
  tetrahedral mesh and extracting/inspecting isosurfaces.
- **`docsrc/`** — Doxygen main-page source; `vgtl.cfg` is the Doxygen config.

## Building

Dependencies: Qt5 (`opengl`, `widgets`, `concurrent`), GSL, and a Fortran
compiler (for the L-BFGS-B routine).

```sh
qmake -o Makefile vgtl2.0.pro
make
```

This builds the static library, every example, and the `isoview` app.

### A note on the vintage

This code was written for compilers and a Qt version current around
2008–2010. It has since been brought back to compile and run cleanly on a
modern toolchain (GCC 16, Qt 5.15, GSL 2.8); the two systemic issues were
`vgtl::array<T,N>` having drifted to a different index type than the rest of
the library uses (breaking template deduction almost everywhere), and
`vgtl::array` colliding with `std::array`, which didn't exist yet when this
was first written. Nothing about the underlying algorithms changed.

## Related

A from-scratch C reimplementation of the topological core (face operators
over oriented pseudomanifolds, fixed at dimension 4, plus Maubach's
subdivision scheme) is being developed separately, for a new application
that only needs that slice of the library.
