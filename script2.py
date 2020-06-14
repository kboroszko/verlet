trips = []

with open("cmake-build-debug/out.txt", "r") as f:
    for line in f.readlines():
        tri = [int(x) for x in line.split()]
        trips.append(tri)
n = 6

found = []
not_found = []

for i in range(n):
    for j in range(n):
        for k in range(n):
            s = {i, j, k}
            print('checking', s)
            if s in trips :
                found.append(s)
                print('FOUND')
            else :
                not_found.append(s)

print()
print('FOUND:')
print(found)
print()
print('not found')
print(not_found)
print()
print('trips len', len(trips))
print('found len', len(found))

unique_tris = []
unique_found = []

for s in found:
    if s not in unique_found :
        unique_found.append(s)

for s in trips :
    if s not in unique_tris:
        unique_tris.append(s)

print()
print('unique trips len', len(unique_tris))
print('unique found len', len(unique_found))
print()
print(unique_found)