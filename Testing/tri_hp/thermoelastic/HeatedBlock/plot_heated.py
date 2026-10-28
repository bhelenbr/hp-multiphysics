#!/usr/bin/env python3
# Plots the deformed heated block (Results/data1_b0.dat, tecplot) colored by T and p
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mt

os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
L = [l for l in open('Results/data1_b0.dat').read().split('\n') if l.strip() and not l.startswith('#')]
N = int(L[0].split('N = ')[1].split(',')[0])
E = int(L[0].split('E = ')[1].split(',')[0])
P = np.array([[float(v) for v in l.split()] for l in L[1:1+N]])
T = np.array([[int(v)-1 for v in l.split()[:3]] for l in L[1+N:1+N+E]])
x, y, u, v, temp, p = P.T
print("top corner displacement: u(1,1)", u[np.isclose(x,1) & np.isclose(y,1)][0], " v(1,1)", v[np.isclose(x,1) & np.isclose(y,1)][0])
print("p range", p.min(), p.max())

fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
for a, (val, name) in zip(ax, [(temp, 'T'), (p, 'p (= -lambda div u +3K alpha (T-T_ref))')]):
	tr = mt.Triangulation(x+u, y+v, T)
	c = a.tripcolor(tr, val, shading='gouraud', cmap='viridis')
	a.triplot(tr, lw=0.3, color='k')
	a.plot([0,1,1,0,0], [0,0,1,1,0], 'r--', lw=0.8)
	a.set_aspect('equal'); a.set_title(name, fontsize=10); fig.colorbar(c, ax=a)
fig.suptitle('Thermoelastic plane strain: clamped base T=0, free top T=1, E=1, nu=0.3, alpha=0.1, p=2, 8x8')
fig.tight_layout()
fig.savefig('heated.png', dpi=130)
