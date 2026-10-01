import itertools
src = open("hybrid.py").read().split("types = [")[0]
exec(src)                     # defines d, children, pair, check
def pos(p, q, l):             # position of omega[q] in child deleting omega[p], v inserted at l
    r = q - (p < q); return r + (r >= l)
def rule(T):
    n = len(T)
    for t in range(n):
        F, l = T[t]; Fn = T[(t+1) % n][0]
        if min(F) < l < max(F): return False                      # (C0)
        if l in Fn: return False                                   # (C1)
        for p, q in itertools.combinations(sorted(F), 2):          # (C2)
            a, b = pos(p, q, l), pos(q, p, l)
            if a != b and not {a, b} <= Fn: return False
    return True
types = [(frozenset(F), l) for s in (2, 3, 4) for F in itertools.combinations(range(d+1), s) for l in range(d+1)]
for n in (1, 2):
    agree = dis = 0
    for T in itertools.product(types, repeat=n):
        if check(T) == rule(T): agree += 1
        else: dis += 1; print("DISAGREE", T) if dis < 5 else None
    print(f"cycle length {n}: agree {agree}, disagree {dis}")
import random
rng = random.Random(3); agree = dis = 0
for _ in range(200000):
    T = tuple(rng.choice(types) for _ in range(rng.choice((3, 4, 5))))
    if check(T) == rule(T): agree += 1
    else: dis += 1
print(f"random cycles of length 3-5: agree {agree}, disagree {dis}")
