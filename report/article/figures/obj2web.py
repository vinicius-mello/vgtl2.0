#!/usr/bin/env python3
"""Convert an OBJ written by riemann_cp2_glpt / riemann_s2s2_glpt into the
compact binary mesh read by the web viewer (vinicius-mello.github.io/riemann).

  obj2web.py OBJ LOG OUTDIR ID [--merge] [--cutoff R] [--meta key=value ...]

Writes OUTDIR/ID.bin -- little-endian uint32 nv, nt, np, then nv*3 float32
positions, nt*3 uint32 indices and np uint8 polygon sizes (polygon k is
triangles sum(size[<k])-2k ... in fan order, so the viewer can draw the
polygon edges without the fan diagonals) -- and adds or replaces the entry
ID in OUTDIR/meshes.json. Polygons are fanned from their first vertex, which
keeps the orientation they were written with (the sign of omega). With
--merge, vertices with identical coordinates are merged: flat mode writes
them per polygon (use it only with a fixed --chart; the OBJ has 6 significant
digits, so distinct but very close points would merge too). With
--cutoff R, polygons with a vertex farther than R from the origin are
dropped. chi and the boundary count are those of the merged polygon mesh;
the numbers parsed from LOG (leaves, bad cells, genus, area) are those of
the run itself, in CP^2 or CP^1 x CP^1.
"""
import json, math, os, re, struct, sys


def parse_log(path):
    s = open(path).read()
    m = {}
    def last(rx, conv=float):
        r = re.findall(rx, s)
        return conv(r[-1]) if r else None
    m['leaves'] = last(r'final leaf count: (\d+)', int) or last(r'SUMMARY leaves=(\d+)', int)
    ok = last(r'extracted cells: ok=(\d+)', int)
    bad = last(r'extracted cells: ok=\d+ bad=(\d+)', int)
    if ok is None:
        ok = last(r'SUMMARY .* ok=(\d+)', int); bad = last(r'SUMMARY .* bad=(\d+)', int)
    m['ok'], m['bad'] = ok, bad
    m['genus_est'] = last(r'genus estimate \(2c-b-chi\)/2 = (-?[\d.]+)')
    m['genus'] = last(r'expected (\d+)', int)
    m['area_ratio'] = last(r'ratio (\d+\.\d+)') or last(r'area_out=(\d+\.\d+)')
    r = re.findall(r'^V=\d+ E=\d+ F=\d+ chi=(-?\d+)', s, re.M)
    m['chi_cp2'] = int(r[-1]) if r else None
    return {k: v for k, v in m.items() if v is not None}


def main():
    a = sys.argv[1:]
    obj, log, outdir, mid = a[:4]
    cutoff, meta, merge = None, {}, False
    i = 4
    while i < len(a):
        if a[i] == '--merge': merge = True; i += 1
        elif a[i] == '--cutoff': cutoff = float(a[i + 1]); i += 2
        elif a[i] == '--meta':
            i += 1
            while i < len(a) and not a[i].startswith('--'):
                k, v = a[i].split('=', 1)
                try: v = json.loads(v)
                except ValueError: pass
                meta[k] = v; i += 1
        else: sys.exit('unknown option ' + a[i])

    V, F = [], []
    for line in open(obj):
        t = line.split()
        if not t: continue
        if t[0] == 'v': V.append((float(t[1]), float(t[2]), float(t[3])))
        elif t[0] == 'f': F.append([int(x.split('/')[0]) - 1 for x in t[1:]])

    if cutoff is not None:
        far = [math.sqrt(x * x + y * y + z * z) > cutoff for x, y, z in V]
        F = [f for f in F if not any(far[v] for v in f)]

    index, P, remap = {}, [], {}
    for f in F:
        for v in f:
            if v in remap: continue
            key = V[v] if merge else v
            if key not in index: index[key] = len(P); P.append(V[v])
            remap[v] = index[key]
    F = [[remap[v] for v in f] for f in F]

    edges = {}
    for f in F:
        for k in range(len(f)):
            e = (min(f[k], f[(k + 1) % len(f)]), max(f[k], f[(k + 1) % len(f)]))
            edges[e] = edges.get(e, 0) + 1
    T = [(f[0], f[k], f[k + 1]) for f in F for k in range(1, len(f) - 1)]

    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, mid + '.bin'), 'wb') as out:
        out.write(struct.pack('<III', len(P), len(T), len(F)))
        out.write(struct.pack('<%df' % (3 * len(P)), *[c for p in P for c in p]))
        out.write(struct.pack('<%dI' % (3 * len(T)), *[v for t in T for v in t]))
        out.write(bytes(len(f) for f in F))

    entry = dict(meta)
    entry.update(parse_log(log))
    entry.update(id=mid, file=mid + '.bin', vertices=len(P), polygons=len(F),
                 triangles=len(T), chi=len(P) - len(edges) + len(F),
                 boundary_edges=sum(1 for c in edges.values() if c == 1),
                 bytes=os.path.getsize(os.path.join(outdir, mid + '.bin')))
    if cutoff is not None: entry['cutoff'] = cutoff
    man = os.path.join(outdir, 'meshes.json')
    L = json.load(open(man)) if os.path.exists(man) else []
    L = [e for e in L if e['id'] != mid] + [entry]
    json.dump(L, open(man, 'w'), indent=1, ensure_ascii=False)
    print(f"{mid}: {len(P)} vertices, {len(F)} polygons, {len(T)} triangles, chi={entry['chi']}, "
          f"boundary edges={entry['boundary_edges']}, {entry['bytes']/1e6:.2f} MB")


if __name__ == '__main__':
    main()
