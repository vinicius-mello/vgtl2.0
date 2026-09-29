"""Render an extracted Riemann-surface OBJ with Blender, in a clean
line-drawing-like style (white background, light grey surface, dark
outline), or with transparency and colour for nested alpha-projection
sheets.

Usage:
  blender -b --python render_obj.py -- IN.obj OUT.png MODE ELEV AZIM [DIST] [RES]
MODE: 'grey' or 'alpha'
"""
import sys, math, os
import bpy, bmesh
from mathutils import Vector

argv = sys.argv[sys.argv.index("--") + 1:]
inp, out, mode = argv[0], argv[1], argv[2]
elev, azim = float(argv[3]), float(argv[4])
dist = float(argv[5]) if len(argv) > 5 else 2.6
res = int(argv[6]) if len(argv) > 6 else 1400

bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.wm.obj_import(filepath=inp)
obj = bpy.context.selected_objects[0]
bpy.context.view_layer.objects.active = obj

# weld the independent per-polygon vertices, smooth normals
me = obj.data
bm = bmesh.new(); bm.from_mesh(me)
bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-6)
bmesh.ops.triangulate(bm, faces=bm.faces[:])
# optional box clip (in the OBJ's own coordinates): CLIP=xmin,xmax,ymin,ymax,zmin,zmax
import os
if os.environ.get("FILL"):
    # close the small holes left by bad cells (disclosed in the caption);
    # done before clipping so the clip boundary is never capped
    r = bmesh.ops.holes_fill(bm, edges=[e for e in bm.edges if e.is_boundary], sides=int(os.environ.get("FILL")))
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    print("holes filled:", len(r["faces"]))
    # holes_fill skips some loops (twisted ones near branch points): close
    # every remaining small boundary loop by a fan from its centroid
    maxs = int(os.environ.get("FILL")); seen = set(); nfan = 0
    for e in [e for e in bm.edges if e.is_boundary]:
        if e in seen or not e.is_valid: continue
        stack = [e]; comp = []
        while stack:
            x = stack.pop()
            if x in seen: continue
            seen.add(x); comp.append(x)
            for v in x.verts:
                for y in v.link_edges:
                    if y.is_boundary and y not in seen: stack.append(y)
        if len(comp) > maxs: continue
        vs = set(v for x in comp for v in x.verts)
        c = bm.verts.new(sum((v.co for v in vs), Vector()) / len(vs))
        for x in comp:
            bm.faces.new((x.verts[0], x.verts[1], c))
        nfan += 1
    print("fans:", nfan)
clip = os.environ.get("CLIP")
if clip:
    lim = [float(t) for t in clip.split(",")]
    for axis in range(3):
        for k, sign in ((0, -1), (1, 1)):
            no = Vector([0, 0, 0]); no[axis] = sign
            co = Vector([0, 0, 0]); co[axis] = lim[2 * axis + k]
            geom = bm.verts[:] + bm.edges[:] + bm.faces[:]
            bmesh.ops.bisect_plane(bm, geom=geom, plane_co=co, plane_no=no, clear_outer=True)
bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])
bm.to_mesh(me); bm.free()
me.update()
bpy.context.view_layer.update()
for p in me.polygons:
    p.use_smooth = True

# centre and normalise
vs = [obj.matrix_world @ v.co for v in me.vertices]
bb = [Vector((min(v[0] for v in vs), min(v[1] for v in vs), min(v[2] for v in vs))),
      Vector((max(v[0] for v in vs), max(v[1] for v in vs), max(v[2] for v in vs)))]
ctr = (bb[0] + bb[1]) / 2
size = max((max(v[i] for v in bb) - min(v[i] for v in bb)) for i in range(3))
obj.location = -ctr
bpy.context.view_layer.update()
bpy.ops.object.transform_apply(location=True)
s = 2.0 / size
obj.scale = (s, s, s)
bpy.ops.object.transform_apply(scale=True)

scene = bpy.context.scene
scene.render.resolution_x = res
scene.render.resolution_y = int(res * float(os.environ.get("ASPECT", "0.8")))
scene.render.film_transparent = False
scene.world = bpy.data.worlds.new("w")
scene.world.color = (1, 1, 1)
scene.view_settings.view_transform = "Standard"

cam_data = bpy.data.cameras.new("cam")
cam_data.type = "ORTHO"
cam_data.ortho_scale = dist
cam = bpy.data.objects.new("cam", cam_data)
scene.collection.objects.link(cam)
e, a = math.radians(elev), math.radians(azim)
cam.location = (6 * math.cos(e) * math.cos(a), 6 * math.cos(e) * math.sin(a), 6 * math.sin(e))
d = -cam.location
cam.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()
scene.camera = cam
if os.environ.get("CAMERA"):
    # CAMERA=x,y,z,rx,ry,rz,ortho_scale (the views chosen for the article
    # are recorded in make_figures.sh)
    c = [float(t) for t in os.environ["CAMERA"].split(",")]
    cam.location = c[0:3]; cam.rotation_euler = c[3:6]; cam_data.ortho_scale = c[6]
if os.environ.get("CAMERA_FROM"):
    # reuse a camera chosen interactively and saved in a .blend
    with bpy.data.libraries.load(os.environ["CAMERA_FROM"]) as (src, dst):
        dst.objects = [n for n in src.objects if n == "cam"]
    c0 = dst.objects[0]
    cam.location = c0.location.copy()
    cam.rotation_mode = c0.rotation_mode
    cam.rotation_euler = c0.rotation_euler.copy()
    cam_data.type = c0.data.type
    cam_data.ortho_scale = c0.data.ortho_scale
    print("camera from", os.environ["CAMERA_FROM"])

if mode == "grey":
    scene.render.engine = "BLENDER_WORKBENCH"
    sh = scene.display.shading
    sh.light = "STUDIO"
    sh.show_backface_culling = False  # open surfaces: both sides visible
    sh.color_type = "SINGLE"
    sh.single_color = (0.82, 0.82, 0.82)
    sh.show_object_outline = True
    sh.object_outline_color = (0, 0, 0)
    sh.show_cavity = False
    sh.show_specular_highlight = False
    sh.background_type = "WORLD"
    scene.display.render_aa = "32"
else:
    # nested alpha-projection sheets: Workbench X-ray (uniform transparency,
    # noise-free), per-vertex colour by distance from the centre
    import colorsys
    d = [v.co.length for v in me.vertices]
    lo, hi = min(d), max(d)
    q = sorted(d); lo, hi = q[int(0.02 * len(q))], q[int(0.98 * len(q))]
    col = me.color_attributes.new(name="Col", type="FLOAT_COLOR", domain="POINT")
    stops = [(0.0, (0.15, 0.30, 0.75)), (0.5, (0.95, 0.85, 0.30)), (1.0, (0.85, 0.25, 0.10))]
    def ramp(t):
        t = min(max(t, 0.0), 1.0)
        for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
            if t <= t1:
                u = (t - t0) / (t1 - t0)
                return [c0[k] + u * (c1[k] - c0[k]) for k in range(3)]
        return list(stops[-1][1])
    for i, v in enumerate(me.vertices):
        c = ramp((d[i] - lo) / (hi - lo))
        col.data[i].color = (c[0], c[1], c[2], 1.0)
    scene.render.engine = "BLENDER_WORKBENCH"
    sh = scene.display.shading
    sh.light = "STUDIO"
    sh.color_type = "VERTEX"
    sh.show_backface_culling = False
    sh.show_xray = True
    sh.xray_alpha = float(argv[9]) if len(argv) > 9 else 0.45
    sh.show_object_outline = False
    sh.show_specular_highlight = False
    sh.background_type = "WORLD"
    scene.display.render_aa = "32"
    print("distance range", lo, hi)

scene.render.filepath = out
if os.environ.get("SAVE_BLEND"):
    # save the prepared scene for interactive choice of the view; the
    # figure is later rendered from the saved camera (render_blend.py)
    bpy.ops.wm.save_as_mainfile(filepath=os.environ["SAVE_BLEND"])
    print("saved", os.environ["SAVE_BLEND"])
else:
    bpy.ops.render.render(write_still=True)
    print("rendered", out)
