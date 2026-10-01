#!/bin/bash
# Fubini-Study area of the raw and closed extraction against Wirtinger's
# n*pi, as a function of depth (continuation only). Output: area.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=area.csv
echo "curve,function,depth,bad,area_ratio_raw,area_ratio_closed" > $OUT
for spec in "conic 0" "elliptic 5" "fermat_cubic 1" "fermat_quartic 3"; do
  set -- $spec
  for d in 10 12 14 16 18; do
    log=$($P --function $2 --generic --depth $d --close --obj /tmp/area_$$.obj 2>&1)
    bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
    ar=$(echo "$log" | grep "^FS area" | head -1 | sed 's/.*ratio \([0-9.e-]*\)).*/\1/')
    ac=$(echo "$log" | grep "^FS area" | tail -1 | sed 's/.*ratio \([0-9.e-]*\)).*/\1/')
    echo "$1,$2,$d,$bad,$ar,$ac" | tee -a $OUT
  done
done
rm -f /tmp/area_$$.obj
