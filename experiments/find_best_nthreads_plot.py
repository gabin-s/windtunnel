from tkinter import N
from turtle import width
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

d = np.loadtxt('find_best_nthreads_narrow_50runs.results', delimiter=' ', comments='#')

nrows, nthreads, t = d.T

fig = plt.figure()
ax = plt.axes(projection='3d', computed_zorder=False)

nrows    = np.unique(nrows)
nthreads = np.unique(nthreads)
t = t.reshape(len(nrows), -1)

"""# filter out data where nthreads >= 128
mask = nthreads < 128
nthreads = nthreads[mask]
t        = t[:, mask]"""

X, Y = np.meshgrid(nthreads, nrows)
ax.plot_surface(X, Y, t, cmap=cm.coolwarm, zorder=0)
ax.set_xlabel('threads')
ax.set_ylabel('rows')

m = np.argmin(t, axis=1)
a = np.arange(len(nrows))

ax.scatter(nthreads[m], nrows, t[a, m]+1e-3, c='black', depthshade=False, zorder=1, label='optimal for $N_{rows}$')
ax.set_zlabel('Runtime (s)')

plt.legend()
plt.tight_layout()

plt.figure()

plt.plot(nrows, nthreads[m], 'x--')
plt.grid()
plt.ylabel('number of threads')
plt.xlabel('number of rows')

plt.show()