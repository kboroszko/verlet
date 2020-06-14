# Patryk Czajka
# refactored to use python3 Kajetan Boroszko

import sys

p = int(input())

def prv(i, p):
    return i - 1 if i > 0 else p - 1

def nxt(i, p):
    return (i + 1) % p

S = set()

def compute(a, b, c):
    global S
    print( a, b, c)
    if a > b: a, b = b, a
    if c < b: c, b = b, c
    if a > b: a, b = b, a
    if (a, b, c) in S:
        print( "REDUNDANT")
        sys.exit(1)
    else:
        S.add((a, b, c))


for r in range(p):
    print("RANK", r)
    b = [prv(r, p), r, nxt(r, p)]
    i = 0

    for s in range(p - 3, -1, -3):
        for move in range(s):
            if move != 0 or s != p - 3:
                b[i] = prv(b[i], p)
            else:
                compute( b[1], b[1], b[1] )
                compute( b[1], b[1], b[2] )
                compute( b[0], b[0], b[2] )
            if s == p - 3:
                compute( b[0], b[1], b[1] )
            compute( b[0], b[1], b[2] )

        i = (i + 1) % 3

    if p % 3 == 0:
        i = prv(i, 3)
        b[i] = prv(b[i], p)
        if (r // (p // 3)) == 0:
            compute(b[0], b[1], b[2])


for i in range(p):
    for j in range(i, p):
        for k in range(j, p):
            if not (i, j, k) in S:
                print( "NO", (i, j, k))