#!/usr/bin/env python3
"""Simulator for stellar subdivision schemes that split 2-faces of 4-simplices.

A scheme is given by:
  face(t)   -> sorted tuple of 3 positions split in a cell of type t
  pos(t)    -> position where the new vertex is inserted
  child(t,q)-> type of child that deleted the q-th face vertex (q indexes face positions)
  ntypes
"""
import itertools, random, sys
from collections import defaultdict

D = 4

def kuhn_torus(n=4):
    """Reflected Kuhn triangulation of the 4-torus Z_n^4 (n even); color = #odd coords."""
    cells = set()
    for x in itertools.product(range(n), repeat=D):
        # cube with lower corner x; even corner / odd corner per coordinate
        ev = tuple(xi if xi % 2 == 0 else (xi + 1) % n for xi in x)
        od = tuple(xi if xi % 2 == 1 else (xi + 1) % n for xi in x)
        for perm in itertools.permutations(range(D)):
            p = list(ev); path = [tuple(p)]
            for i in perm:
                p[i] = od[i]; path.append(tuple(p))
            cells.add(tuple(path))
    return list(cells)

class Mesh:
    def __init__(self, scheme, seed_cells):
        self.S = scheme
        self.cells = {}          # id -> (verts tuple, level, type)
        self.vstar = defaultdict(set)
        self.nid = 0
        self.nv = 0
        for c in seed_cells:
            self.add(c, 0, scheme['t0'])
    def add(self, verts, lev, t):
        i = self.nid; self.nid += 1
        self.cells[i] = (verts, lev, t)
        for v in verts: self.vstar[v].add(i)
        return i
    def remove(self, i):
        verts = self.cells.pop(i)[0]
        for v in verts: self.vstar[v].discard(i)
    def split_face(self, i):
        verts, lev, t = self.cells[i]
        return frozenset(verts[p] for p in self.S['face'][t])
    def star(self, fs):
        it = iter(fs); s = set(self.vstar[next(it)])
        for v in it: s &= self.vstar[v]
        return s
    def subdivide(self, fs, cells):
        v = ('n', self.nv); self.nv += 1
        new = []
        for i in cells:
            verts, lev, t = self.cells[i]
            P = self.S['face'][t]; l = self.S['pos'][t]
            assert frozenset(verts[p] for p in P) == fs
            self.remove(i)
            for q, p in enumerate(P):
                rest = [verts[r] for r in range(D + 1) if r != p]
                rest.insert(l, v)
                new.append(self.add(tuple(rest), lev + 1, self.S['child'](t, q)))
        return new
    def compatible_star(self, i):
        fs = self.split_face(i); lev = self.cells[i][1]
        st = self.star(fs)
        bad = [j for j in st if self.split_face(j) != fs or self.cells[j][1] != lev]
        return fs, st, bad

def uniform(scheme, seed, rounds):
    M = Mesh(scheme, seed)
    for r in range(rounds):
        ids = list(M.cells)
        nbad = 0; ex = None
        for i in ids:
            fs, st, bad = M.compatible_star(i)
            if bad:
                nbad += 1
                if ex is None: ex = (M.cells[i], [M.cells[j] for j in bad][:2])
        print(f"  round {r}: {len(ids)} cells, incompatible cells: {nbad}")
        if nbad:
            print("   example:", ex); return False
        done = set()
        for i in ids:
            if i in M.cells:
                fs = M.split_face(i)
                M.subdivide(fs, M.star(fs))
    return True

class Fail(Exception): pass

def refine(M, i, depth=0, cap=200):
    if depth > cap: raise Fail("recursion cap")
    while i in M.cells:
        fs, st, bad = M.compatible_star(i)
        if not bad:
            M.subdivide(fs, st); return
        lev = M.cells[i][1]
        j = min(bad, key=lambda j: M.cells[j][1])
        if M.cells[j][1] >= lev:
            raise Fail(f"same/higher-level incompatibility: {M.cells[i]} vs {M.cells[j]}")
        refine(M, j, depth + 1, cap)

def adaptive(scheme, seed, steps, rng):
    M = Mesh(scheme, seed)
    target = (0.1, 0.2, 0.3, 0.1)
    for s in range(steps):
        # refine the finest cells near a moving point, emulate local refinement
        ids = list(M.cells)
        i = rng.choice(ids)
        try:
            refine(M, i)
        except Fail as e:
            print(f"  adaptive FAIL at step {s}: {e}"); return False
    levs = [c[1] for c in M.cells.values()]
    print(f"  adaptive ok: {steps} steps, {len(M.cells)} cells, max level {max(levs)}")
    return True

def run(name, S, rounds=3, steps=3000):
    print(name)
    seed = kuhn_torus(4)
    ok = uniform(S, seed, rounds)
    ok2 = adaptive(S, seed, steps, random.Random(1))
    return ok, ok2

if __name__ == '__main__':
    user = dict(face={0: (0, 1, 4), 1: (0, 1, 3), 2: (0, 1, 2)}, pos={0: 4, 1: 3, 2: 2}, t0=0)
    for k in (0, 1, 2):
        S = dict(user, child=(lambda k: lambda t, q: (t + k) % 3)(k))
        run(f"user scheme, child type = t+{k} mod 3", S)
