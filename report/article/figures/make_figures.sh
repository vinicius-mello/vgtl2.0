#!/bin/bash
# Regenerates every figure of the article from scratch.
#   seed_crossing.pdf: fig_seed_crossing.py (needs matplotlib)
#   faces.pdf:         fig_faces.py (needs matplotlib)
#   elliptic_flat.png, quartic_alpha.png: riemann_cp2_glpt --close, then
#             Blender (views chosen interactively, recorded below as
#             CAMERA=x,y,z,rx,ry,rz,ortho)
#   badrate.pdf:       badrate_sweep.sh + fig_badrate.py
#   Tables:   causes_sweep.sh, genus_sweep.sh, rotation_sweep.sh, modes_table.sh,
#             continuation_check.sh, cost_sweep.sh, area_sweep.sh, affine_sweep.sh,
#             s2s2_sweep.sh (CP^1 x CP^1, examples/top/riemann_s2s2_glpt),
#             degree_sweep.sh (degree 3-6: Fermat and Kostlan random curves)
# (badrate_chartflat.csv / rotations_chartflat.csv keep the measurements
#  made with the old chart-flat faces, before geodesic faces became default.)
# Also prints the genus check of both rendered meshes.
set -e
cd "$(dirname "$0")"
P=../../../examples/top/riemann_cp2_glpt/riemann_cp2_glpt
PY=${PYTHON:-python3}

$PY fig_seed_crossing.py
$PY fig_faces.py

$P --function 5 --generic --depth 18 --tangency-threshold 0.05 --repair-rounds 2 \
   --chart 2 --flat-swap --cutoff 8 --close --obj elliptic_flat.obj | sed -n '/--- genus check/,$p'
$P --function 3 --generic --depth 17 --tangency-threshold 0.05 --repair-rounds 2 \
   --alpha-projection --close --obj quartic_alpha.obj | sed -n '/--- genus check/,$p'

ASPECT=0.8 CLIP=-2,2,-2,2,-3,3 \
CAMERA=-2.455743,-3.676821,4.041286,0.830780,0.000001,-0.593414,3.000000 \
  blender -b --python render_obj.py -- elliptic_flat.obj "$PWD/elliptic_flat.png" grey 30 -60 3.0 1800
ASPECT=1.0 \
CAMERA=2.969030,0.543071,-5.164394,2.614504,0.000000,1.755801,3.600000 \
  blender -b --python render_obj.py -- quartic_alpha.obj "$PWD/quartic_alpha.png" alpha 25 -45 3.6 1800 0 0 0.55
# crop the white margins
$PY - <<'PYEOF'
from PIL import Image
for f, th in [("elliptic_flat.png", 245), ("quartic_alpha.png", 248)]:
    im = Image.open(f).convert("RGB")
    bb = im.convert("L").point(lambda x: 255 if x < th else 0).getbbox(); pad = 20
    bb = (max(0, bb[0]-pad), max(0, bb[1]-pad), min(im.width, bb[2]+pad), min(im.height, bb[3]+pad))
    im.crop(bb).save(f)
PYEOF

./badrate_sweep.sh
$PY fig_badrate.py
./causes_sweep.sh
./genus_sweep.sh
./rotation_sweep.sh
./modes_table.sh
./continuation_check.sh
./cost_sweep.sh
./area_sweep.sh
./affine_sweep.sh
./s2s2_sweep.sh
./degree_sweep.sh
