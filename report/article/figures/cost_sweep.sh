#!/bin/bash
# Running time per phase and peak memory (one core).
# Output: cost.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=cost.csv
echo "curve,function,depth,mode,leaves,touched,bad,refine_s,extract_s,peak_mb" > $OUT
run() { name=$1; fn=$2; d=$3; mode=$4; shift 4
  log=$($P --function $fn --generic --depth $d --timing "$@" --obj /tmp/cost_$$.obj 2>&1)
  leaves=$(echo "$log" | grep "final leaf count" | tail -1 | awk '{print $4}')
  ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
  bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
  t=$(echo "$log" | grep "^timing" | sed 's/timing: refinement \([0-9.e-]*\) s, extraction (incl. repair) + post-processing \([0-9.e-]*\) s, peak RSS \([0-9.e-]*\) MB/\1,\2,\3/')
  echo "$name,$fn,$d,$mode,$leaves,$((ok+bad)),$bad,$t" | tee -a $OUT
}
for spec in "conic 0" "elliptic 5" "fermat_cubic 1" "fermat_quartic 3"; do
  set -- $spec
  for d in 12 14 16 18; do
    run $1 $2 $d continuation
    run $1 $2 $d both --tangency-threshold 0.05 --repair-rounds 2
  done
done
run elliptic 5 20 continuation
rm -f /tmp/cost_$$.obj
