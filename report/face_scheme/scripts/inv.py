import random, itertools
from collections import Counter, defaultdict
from facesplit import *
user = dict(face={0:(0,1,4),1:(0,1,3),2:(0,1,2)}, pos={0:4,1:3,2:2}, t0=0, child=lambda t,q:(t+1)%3)

def pairs(M):
    fac = defaultdict(list)
    for i,(v,l,t) in M.cells.items():
        for p in range(5):
            fac[frozenset(v[:p]+v[p+1:])].append((i,p))
    C = Counter()
    for f,lst in fac.items():
        if len(lst)==2:
            (a,i),(b,j) = lst
            va,la,ta = M.cells[a]; vb,lb,tb = M.cells[b]
            # check ordered equality
            ordered = (va[:i]+va[i+1:]) == (vb[:j]+vb[j+1:])
            if la<lb or (la==lb and (i,ta)>(j,tb)): a,b,i,j,la,lb,ta,tb=b,a,j,i,lb,la,tb,ta
            C[(lb-la, ta, i, tb, j, ordered)] += 1
    return C

seed = kuhn_torus(4)
M = Mesh(user, seed)
for r in range(3):
    ids=list(M.cells)
    for i in ids:
        if i in M.cells:
            fs=M.split_face(i); M.subdivide(fs,M.star(fs))
    C=pairs(M)
    print("uniform round",r+1, "pairs (dlev, t_coarse, i, t_fine, j, ordered-equal):")
    for k,v in sorted(C.items()): print("   ",k,v)
# adaptive
M = Mesh(user, seed); rng=random.Random(7)
for s in range(6000):
    refine(M, rng.choice(list(M.cells)))
C=pairs(M)
print("adaptive:", len(M.cells), "max lev", max(c[1] for c in M.cells.values()))
for k,v in sorted(C.items()): print("   ",k,v)
