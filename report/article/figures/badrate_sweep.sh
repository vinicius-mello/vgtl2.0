#!/bin/bash
# Bad-cell rate versus depth (data for the article's bad-rate figure).
# Output: badrate.csv  (curve,function,depth,mode,leaves,ok,bad)
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=badrate.csv
echo "curve,function,depth,mode,leaves,ok,bad" > $OUT
run() { # name fn depth mode flags...
  name=$1; fn=$2; d=$3; mode=$4; shift 4
  log=$($P --function $fn --generic --depth $d "$@" --obj /tmp/badrate_$$.obj 2>&1)
  leaves=$(echo "$log" | grep "final leaf count" | tail -1 | awk '{print $4}')
  ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
  bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
  echo "$name,$fn,$d,$mode,$leaves,$ok,$bad" | tee -a $OUT
}
for spec in "conic 0" "elliptic 5" "fermat_cubic 1" "fermat_quartic 3"; do
  set -- $spec; name=$1; fn=$2
  for d in 8 10 12 14 16; do
    run $name $fn $d baseline
    run $name $fn $d combined --tangency-threshold 0.05 --repair-rounds 2
  done
done
run elliptic 5 18 baseline
run elliptic 5 18 combined --tangency-threshold 0.05 --repair-rounds 2
rm -f /tmp/badrate_$$.obj
