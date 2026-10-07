#!/bin/bash
# CP^1 x CP^1 with the same pipeline (examples/top/riemann_s2s2_glpt): the
# elliptic curve w^2=z^3-z, whose closure has a cusp at (inf,inf), against
# the parabola w^2=z, smooth there. Base mode, geodesic cones, Kahler filter
# 0.9, --close. "near" = cells with a vertex within 0.1 of (inf,inf) in the
# product Fubini-Study metric. Depth sweep with the default rotation, then
# five rotations at depth 14.
# Output: s2s2.csv
P=../../../examples/top/riemann_s2s2_glpt/riemann_s2s2_glpt
OUT=s2s2.csv
echo "curve,function,depth,seed,leaves,ok,bad,kahler,near_ok,near_bad,area_raw,area_closed,closed_chi,closed_b,closed_c,orientable,genus_est" > $OUT
run() { # name fn depth seed
  log=$($P --function $2 --depth $3 --generic-seed $4 --close --kahler-filter 0.9 2>&1)
  s=$(echo "$log" | grep "^SUMMARY")
  get() { echo "$s" | sed "s/.* $1=\([^ ]*\).*/\1/"; }
  last=$(echo "$log" | grep "^V=" | tail -1)
  chi=$(echo "$last" | sed 's/.* chi=\([-0-9]*\).*/\1/')
  b=$(echo "$last" | sed 's/.*loops b=\([0-9]*\).*/\1/')
  c=$(echo "$log" | grep "components c=" | tail -1 | sed 's/.*components c=\([0-9]*\).*/\1/')
  orient=$(echo "$log" | grep "orientation conflicts" | tail -1 | grep -q "NOT" && echo no || echo yes)
  g=$(echo "$log" | grep "^genus estimate" | tail -1 | awk '{print $5}')
  echo "$1,$2,$3,$4,$(get leaves),$(get ok),$(get bad),$(get kahler),$(get near_ok),$(get near_bad),$(get area_raw),$(get area_out),$chi,$b,$c,$orient,$g" | tee -a $OUT
}
for spec in "elliptic 0" "parabola 1"; do
  set -- $spec
  for d in 10 12 14 16; do run $1 $2 $d 12345; done
  for seed in 1 2 3 4 5; do run $1 $2 14 $seed; done
done
