#!/usr/bin/env python3
# Plots the solid (block 1) of a facet_auto thermal stress result: temperature,
# p (in units of E; sigma = -p I +mu*(grad(u)+grad(u)^T)) and the x displacement
# usage: plot_stress.py Results/dataN_b1.dat
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mt

fname = sys.argv[1]

L = [l for l in open(fname).read().split('\n') if l.strip() and not l.startswith('#')]
N = int(L[0].split('N = ')[1].split(',')[0])
E = int(L[0].split('E = ')[1].split(',')[0])
P = np.array([[float(v) for v in l.split()] for l in L[1:1+N]])
T = np.array([[int(v)-1 for v in l.split()[:3]] for l in L[1+N:1+N+E]])
x, y, u, v, temp, p = P.T[:6]

tr = mt.Triangulation(x, y, T)
print("max |u|,|v|:", abs(u).max(), abs(v).max())
print("p range:", p.min(), p.max())

fig, ax = plt.subplots(3, 1, figsize=(10, 7))
for a, (val, name) in zip(ax, [(temp, 'T'), (p, 'p / E'), (u, 'u (x displacement)')]):
	c = a.tripcolor(tr, val, shading='gouraud', cmap='viridis')
	a.set_title(name, fontsize=10); fig.colorbar(c, ax=a)  # y is stretched: the ribbon is thin
	a.set_xlim(x.min(), x.max())
fig.suptitle('facet_auto solid with thermal stress: ' + fname)
fig.tight_layout()
fig.savefig(fname.replace('.dat', '_stress.png'), dpi=130)
