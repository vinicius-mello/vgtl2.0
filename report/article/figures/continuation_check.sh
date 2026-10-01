#!/bin/bash
# Continuation from one root vs. from all 108 roots (--global-scan):
# compares the final leaf sets by an order-independent hash, and the mesh
# size against the earlier priority-queue criterion (--legacy-priority-refine).
# Output: continuation.csv
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
OUT=continuation.csv
echo "curve,function,depth,leaves,bad,hash,global_leaves,global_bad,global_hash,same,legacy_leaves" > $OUT
for spec in "conic 0" "elliptic 5" "fermat_cubic 1" "fermat_quartic 3"; do
  set -- $spec; name=$1; fn=$2
  for d in 8 12 14; do
    a=$($P --function $fn --generic --depth $d --leaf-hash --obj /tmp/cc_$$.obj 2>&1)
    b=$($P --function $fn --generic --depth $d --leaf-hash --global-scan --obj /tmp/cc_$$.obj 2>&1)
    c=$($P --function $fn --generic --depth $d --legacy-priority-refine --obj /tmp/cc_$$.obj 2>&1)
    la=$(echo "$a" | grep "^leaf set" | awk '{print $3}'); ha=$(echo "$a" | grep "^leaf set" | awk '{print $NF}')
    lb=$(echo "$b" | grep "^leaf set" | awk '{print $3}'); hb=$(echo "$b" | grep "^leaf set" | awk '{print $NF}')
    ba=$(echo "$a" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
    bb=$(echo "$b" | grep "extracted cells" | tail -1 | sed 's/.*bad=\([0-9]*\).*/\1/')
    lc=$(echo "$c" | grep "final leaf count" | tail -1 | awk '{print $4}')
    same=$([ "$ha" = "$hb" ] && echo yes || echo no)
    echo "$name,$fn,$d,$la,$ba,$ha,$lb,$bb,$hb,$same,$lc" | tee -a $OUT
  done
done
rm -f /tmp/cc_$$.obj
