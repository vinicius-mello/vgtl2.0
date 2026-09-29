"""Render a .blend saved by render_obj.py (SAVE_BLEND) from its own camera.
Usage: blender -b FILE.blend --python render_blend.py -- OUT.png [RES]"""
import sys, bpy
argv = sys.argv[sys.argv.index("--") + 1:]
sc = bpy.context.scene
if len(argv) > 1:
    r = int(argv[1]); ar = sc.render.resolution_y / sc.render.resolution_x
    sc.render.resolution_x = r; sc.render.resolution_y = int(r * ar)
sc.display.shading.show_backface_culling = False  # open surface: both sides
sc.render.filepath = argv[0]
bpy.ops.render.render(write_still=True)
print("rendered", argv[0])
