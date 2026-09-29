#!/bin/bash
# Robustness of the Riemann-surface example across random rotations of the seed.
# For each run: bad-cell count, genus estimate of the raw extraction
# (--genus-check), and the closed surface produced by --close
# (chi, boundary loops b, components c).
# Output: rotations.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=rotations.csv
echo "curve,function,depth,seed,ok,bad,genus_est,genus_expected,closed_chi,closed_b,closed_c" > $OUT
for spec in "elliptic 5" "fermat_quartic 3"; do
  set -- $spec; name=$1; fn=$2
  for seed in 1 2 3 4 5; do
    log=$($P --function $fn --generic-seed $seed --depth 14 --tangency-threshold 0.05 \
          --repair-rounds 2 --close --obj /tmp/rot_$$.obj 2>&1)
    ok=$(echo "$log"  | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
    bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
    ge_line=$(echo "$log" | grep "^genus estimate" | head -1)
    g=$(echo "$ge_line"  | awk '{print $5}')
    gx=$(echo "$ge_line" | sed 's/.*(expected \([0-9]*\).*/\1/')
    last=$(echo "$log" | grep "^V=" | tail -1)
    chi=$(echo "$last" | sed 's/.* chi=\([-0-9]*\).*/\1/')
    b=$(echo "$last"   | sed 's/.*loops b=\([0-9]*\).*/\1/')
    c=$(echo "$log" | grep "components c=" | tail -1 | sed 's/.*components c=\([0-9]*\).*/\1/')
    echo "$name,$fn,14,$seed,$ok,$bad,$g,$gx,$chi,$b,$c" | tee -a $OUT
  done
done
rm -f /tmp/rot_$$.obj
