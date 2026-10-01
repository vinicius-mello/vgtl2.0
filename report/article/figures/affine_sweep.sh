#!/bin/bash
# E1: intrinsic CP^2 mesh vs. a Kuhn box [-R,R]^4 in an affine chart
# (same code: continuation, extraction, verification), at equal
# Fubini-Study resolution h (every cell meeting C refined to FS
# diameter < h). Generic seed 12345 in both: the Gaifullin seed is
# rotated, or the curve is rotated relative to the axis-aligned box.
# Output: affine.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=affine.csv
echo "curve,function,mesh,R,h,leaves,touched,bad,c,b,chi,area_ratio,refine_s,extract_s,peak_mb" > $OUT
run() { name=$1; fn=$2; mesh=$3; R=$4; h=$5; shift 5
  log=$($P --function $fn --generic --depth 4 --fs-diam $h --genus-check --timing "$@" --obj /tmp/aff_$$.obj 2>&1)
  leaves=$(echo "$log" | grep "final leaf count" | tail -1 | awk '{print $4}')
  ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
  bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
  c=$(echo "$log" | grep "components c=" | head -1 | sed 's/.*components c=\([0-9]*\).*/\1/')
  v=$(echo "$log" | grep "^V=" | head -1)
  b=$(echo "$v" | sed 's/.*loops b=\([0-9]*\).*/\1/'); chi=$(echo "$v" | sed 's/.* chi=\([-0-9]*\).*/\1/')
  ar=$(echo "$log" | grep "^FS area" | head -1 | sed 's/.*ratio \([0-9.e-]*\)).*/\1/')
  t=$(echo "$log" | grep "^timing" | sed 's/timing: refinement \([0-9.e-]*\) s, extraction (incl. repair) + post-processing \([0-9.e-]*\) s, peak RSS \([0-9.e-]*\) MB/\1,\2,\3/')
  echo "$name,$fn,$mesh,$R,$h,$leaves,$((ok+bad)),$bad,$c,$b,$chi,$ar,$t" | tee -a $OUT
}
for spec in "elliptic 5" "fermat_quartic 3"; do
  set -- $spec; name=$1; fn=$2
  for h in 0.3 0.2; do
    run $name $fn cp2 - $h
    for R in 2 4 8 16; do run $name $fn box $R $h --box $R; done
  done
done
rm -f /tmp/aff_$$.obj
