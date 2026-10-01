import itertools
def pos(p, q, l):
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
