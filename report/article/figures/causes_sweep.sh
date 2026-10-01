#!/bin/bash
# Why the remaining bad cells fail (default geodesic faces), depth 14,
# smooth curves only: failure reason and transversality of good vs bad.
# Output: causes.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=causes.csv
echo "curve,function,depth,ok,bad,irregular_facet,node_degree,bad_cycle,tv_good_median,tv_bad_median,bad_frac_tv_below_0.05" > $OUT
for spec in "conic 0" "elliptic 5" "fermat_cubic 1" "fermat_quartic 3"; do
  set -- $spec; name=$1; fn=$2
  log=$($P --function $fn --generic --depth 14 --cause-stats --obj /tmp/causes_$$.obj 2>&1)
  ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
  bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
  r=$(echo "$log" | grep "^bad cells: irregular" | sed 's/.*facet=\([0-9]*\) node degree=\([0-9]*\) bad cycle=\([0-9]*\).*/\1,\2,\3/')
  tg=$(echo "$log" | grep "^good cells: trans" | sed 's/.*median=\([^ ]*\).*/\1/')
  tb=$(echo "$log" | grep "^bad  cells: trans" | sed 's/.*median=\([^ ]*\).*/\1/')
  fb=$(echo "$log" | grep "^bad  cells: trans" | sed 's/.*fraction<0.05=\([^ ]*\).*/\1/')
  echo "$name,$fn,14,$ok,$bad,$r,$tg,$tb,$fb" | tee -a $OUT
done
rm -f /tmp/causes_$$.obj
