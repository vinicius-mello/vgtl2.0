#!/bin/bash
# Genus of the raw extraction, closing, Fubini-Study area (Wirtinger) and
# vertex residuals, for the curves of the article (tangency threshold 0.05,
# two repair rounds). Output: genus.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=genus.csv
echo "curve,function,depth,ok,bad,pinched,c,b,chi,genus_est,genus,closed_chi,closed_b,closed_c,area_ratio_raw,area_ratio_closed,F_median,F_p99,F_max" > $OUT
for spec in "conic 0 14" "elliptic 5 14" "elliptic 5 18" "fermat_cubic 1 14" "fermat_quartic 3 17"; do
  set -- $spec; name=$1; fn=$2; d=$3
  log=$($P --function $fn --generic --depth $d --tangency-threshold 0.05 --repair-rounds 2 --close --obj /tmp/genus_$$.obj 2>&1)
  ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
  bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
  pin=$(echo "$log" | grep "components c=" | head -1 | sed 's/.*pinched vertices split: \([0-9]*\).*/\1/')
  c=$(echo "$log" | grep "components c=" | head -1 | sed 's/.*components c=\([0-9]*\).*/\1/')
  first=$(echo "$log" | grep "^V=" | head -1); last=$(echo "$log" | grep "^V=" | tail -1)
  chi=$(echo "$first" | sed 's/.* chi=\([-0-9]*\).*/\1/'); b=$(echo "$first" | sed 's/.*loops b=\([0-9]*\).*/\1/')
  ge=$(echo "$log" | grep "^genus estimate" | head -1 | awk '{print $5}')
  gx=$(echo "$log" | grep "^genus estimate" | head -1 | sed 's/.*(expected \([0-9]*\).*/\1/')
  cchi=$(echo "$last" | sed 's/.* chi=\([-0-9]*\).*/\1/'); cb=$(echo "$last" | sed 's/.*loops b=\([0-9]*\).*/\1/')
  cc=$(echo "$log" | grep "components c=" | tail -1 | sed 's/.*components c=\([0-9]*\).*/\1/')
  ar=$(echo "$log" | grep "^FS area" | head -1 | sed 's/.*ratio \([0-9.e-]*\)).*/\1/')
  ac=$(echo "$log" | grep "^FS area" | tail -1 | sed 's/.*ratio \([0-9.e-]*\)).*/\1/')
  fl=$(echo "$log" | grep "^|F| at" | tail -1 | sed 's/.*median=\([^ ]*\) p99=\([^ ]*\) max=\([^ ]*\).*/\1,\2,\3/')
  echo "$name,$fn,$d,$ok,$bad,$pin,$c,$b,$chi,$ge,$gx,$cchi,$cb,$cc,$ar,$ac,$fl" | tee -a $OUT
done
rm -f /tmp/genus_$$.obj
