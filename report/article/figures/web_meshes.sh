#!/bin/bash
# Meshes for the web viewer (vinicius-mello.github.io/riemann): runs the
# pipeline with the flags of the article's figures and converts each OBJ
# with obj2web.py. Usage: web_meshes.sh OUTDIR   (e.g. .../riemann/meshes)
# The result tables of the viewer are the CSVs of this directory, copied
# to OUTDIR/../data.
set -e
cd "$(dirname "$0")"
OUT=$(realpath -m "${1:?usage: web_meshes.sh OUTDIR}")
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
S=../../../examples/top/riemann_s2s2_glpt/riemann_s2s2_glpt
PY=${PYTHON:-python3}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FLAGS="--generic --tangency-threshold 0.05 --repair-rounds 2 --close"

cp2() { # id fn depth title equation degree [extra program flags...] -- [extra obj2web flags...]
  local id=$1 fn=$2 d=$3 title=$4 eq=$5 n=$6; shift 6
  local pf=() cf=()
  while [ $# -gt 0 ] && [ "$1" != -- ]; do pf+=("$1"); shift; done; [ "$1" = -- ] && shift; cf=("$@")
  local proj=alpha; for f in "${pf[@]}"; do [ "$f" = --chart ] && proj=chart; done
  [ $proj = alpha ] && pf+=(--alpha-projection)
  $P --function $fn --depth $d $FLAGS "${pf[@]}" --obj "$TMP/$id.obj" > "$TMP/$id.log" 2>&1
  $PY obj2web.py "$TMP/$id.obj" "$TMP/$id.log" "$OUT" $id "${cf[@]}" --meta \
    "space=CP2" "projection=$proj" "title=$title" "equation=$eq" "degree=$n" "depth=$d"
}
s2s2() { # id fn depth title equation bidegree genus (geometric: the elliptic curve has a cusp at (inf,inf))
  $S --function $2 --depth $3 --generic-seed 12345 --close --kahler-filter 0.9 --obj "$TMP/$1.obj" > "$TMP/$1.log" 2>&1
  $PY obj2web.py "$TMP/$1.obj" "$TMP/$1.log" "$OUT" $1 --cutoff 4 --meta \
    "space=CP1xCP1" "projection=chart" "title=$4" "equation=$5" "bidegree=$6" "genus=$7" "depth=$3"
}

cp2 conic          0 16 "Conic"              "XY - Z^2 = 0"              2
cp2 elliptic       5 16 "Elliptic curve"     "X^2 Z - Y^3 + Y Z^2 = 0"   3
cp2 fermat_quartic 3 16 "Fermat quartic"     "X^4 + Y^4 + Z^4 = 0"       4
cp2 fermat_5       6 14 "Fermat quintic"     "X^5 + Y^5 + Z^5 = 0"       5
cp2 fermat_6       7 14 "Fermat sextic"      "X^6 + Y^6 + Z^6 = 0"       6
cp2 kostlan_6     11 14 "Kostlan random sextic" "random, degree 6 (seed 2032)" 6
cp2 elliptic_chart 5 16 "Elliptic curve (chart Z=1)" "w^2 = z^3 - z"     3 \
    --chart 2 --flat-swap --cutoff 4
s2s2 s2s2_elliptic 0 14 "Elliptic curve" "w^2 = z^3 - z" "[2,3]" 1
s2s2 s2s2_parabola 1 14 "Parabola"       "w^2 = z"       "[2,1]" 0

mkdir -p "$OUT/../data"
cp genus.csv degree.csv s2s2.csv badrate.csv "$OUT/../data/"
