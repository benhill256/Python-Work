# This is original code developed to present my findings regarding patterns in the vertices of n-dimensional cubes.
# This program is good enough to visualize dimensions 1-16 in my experience, but beyond that the calculations are too intensive.

import numpy as np
import matplotlib.pyplot as plt

n = int(input("Input the desired number of dimensions: "))
d = {}
list_dict = {}
for i in range(n):
    list_dict["list_"+str(i)] = []


for i in range(2 ** n):
    b = bin(i)[2:]
    l = len(b)
    b = str(0) * (n - l) + b

    u = []
    for j in range(len(b)):
        u.append(int(b[j]))

    d[i] = [u, []]


keys = d.keys()


for i in d:
    m = d[i][0]

    # Create a list of n zeroes.
    g = [0 for i in range(n)]

    # Turn g into a copy of m, keeping m and g independent.
    for l in range(n):
        g[l] = m[l]
    
    # Compare m and g entrywise, changing one entry in g to be different from m.
    for j in range(n):
        if m[j] == 0:
            g[j] = 1
        elif m[j] == 1:
            g[j] = 0
        
        # Find the decimal number that the binary number g now corresponds to
        # and add that number to the original entry in d, because those vertices are connected.
        for k in keys:
            if d[k][0] == g:
                d[i][1].append(k)
        
        # Reset g to be equal to m, so that it can be compared to other numbers again.
        if g[j] == 0:
            g[j] = 1
        elif g[j] == 1:
            g[j] = 0

    # Sort numerically the numbers of the adjacent vertices.
    d[i][1] = sorted(d[i][1])

    for c in range(n):
        list_dict["list_"+str(c)].append(d[i][1][c])


for i in range(n):
    plt.plot(keys, list_dict["list_"+str(i)])

plt.show()