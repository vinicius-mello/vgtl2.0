#!/bin/bash
# Elliptic curve, depths 14 and 20: effect of repair rounds and of the
# tangency threshold, and the old chart-flat faces for reference.
# Output: modes.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=modes.csv
echo "depth,mode,leaves,ok,bad" > $OUT
run() { d=$1; mode=$2; shift 2
  log=$($P --function 5 --generic --depth $d "$@" --obj /tmp/modes_$$.obj 2>&1)
  leaves=$(echo "$log" | grep "final leaf count" | tail -1 | awk '{print $4}')
  ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
  bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
  echo "$d,$mode,$leaves,$ok,$bad" | tee -a $OUT
}
for d in 14 20; do
  run $d chart_flat --chart-flat-faces
  run $d continuation
  run $d repair --repair-rounds 2
  run $d tangency --tangency-threshold 0.05
  run $d both --tangency-threshold 0.05 --repair-rounds 2
done
rm -f /tmp/modes_$$.obj
