#!/bin/bash
# E4: curves of degree 3 to 6 -- Fermat curves and Kostlan random curves
# (functions_cp2.hpp) -- at depths 12, 14, 16, in base mode and with
# tangency threshold 0.05 + two repair rounds, each closed (--close).
# Output: degree.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=degree.csv
echo "curve,function,degree,depth,mode,leaves,touched,bad,genus_raw,genus_closed,closed_chi,closed_b,closed_c,area_raw,area_closed,F_max,refine_s,extract_s,peak_mib" > $OUT
for spec in "fermat_cubic 1 3" "fermat_quartic 3 4" "fermat_5 6 5" "fermat_6 7 6" \
            "kostlan_3 8 3" "kostlan_4 9 4" "kostlan_5 10 5" "kostlan_6 11 6"; do
  set -- $spec; name=$1; fn=$2; deg=$3
  for d in 12 14 16; do
    for mode in base both; do
      extra=""; [ $mode = both ] && extra="--tangency-threshold 0.05 --repair-rounds 2"
      log=$($P --function $fn --generic --depth $d $extra --close --timing --obj /tmp/deg_$$.obj 2>&1)
      leaves=$(echo "$log" | grep "final leaf count" | tail -1 | awk '{print $4}')
      ok=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*ok=\([0-9]*\).*/\1/')
      bad=$(echo "$log" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
      g0=$(echo "$log" | grep "^genus estimate" | head -1 | awk '{print $5}')
      g1=$(echo "$log" | grep "^genus estimate" | tail -1 | awk '{print $5}')
      last=$(echo "$log" | grep "^V=" | tail -1)
      chi=$(echo "$last" | sed 's/.* chi=\([-0-9]*\).*/\1/')
      b=$(echo "$last" | sed 's/.*loops b=\([0-9]*\).*/\1/')
      c=$(echo "$log" | grep "components c=" | tail -1 | sed 's/.*components c=\([0-9]*\).*/\1/')
      a0=$(echo "$log" | grep "^FS area" | head -1 | sed 's/.*ratio \([0-9.]*\).*/\1/')
      a1=$(echo "$log" | grep "^FS area" | tail -1 | sed 's/.*ratio \([0-9.]*\).*/\1/')
      fmax=$(echo "$log" | grep "^|F| at" | head -1 | sed 's/.*max=\([^ ]*\).*/\1/')
      t=$(echo "$log" | grep "^timing" | sed 's/timing: refinement \([0-9.]*\) s, extraction (incl. repair) + post-processing \([0-9.]*\) s, peak RSS \([0-9.]*\) MB/\1,\2,\3/')
      echo "$name,$fn,$deg,$d,$mode,$leaves,$((ok+bad)),$bad,$g0,$g1,$chi,$b,$c,$a0,$a1,$fmax,$t" | tee -a $OUT
    done
  done
done
rm -f /tmp/deg_$$.obj
