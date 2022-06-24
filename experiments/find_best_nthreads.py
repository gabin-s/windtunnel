#!/usr/bin/env python3
import itertools
import subprocess, os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as cm

N_RUNS=5
EXEC_PATH='.././wind_omp'

max_iter=1000
threshold=-1 # threshold=-1 to prevent early stopping
inlet_pos=0
inlet_size=100

cols=1000

args_t=f"{{rows}} {cols} {max_iter} {threshold} {inlet_pos} {inlet_size} 0 0 0 0 0 0 123 456 789"  # no particles

print(f"# running {EXEC_PATH}, with <nthreads> thread(s)")
print(f"# args: {args_t}")

env = os.environ.copy()

rows_l     = np.arange(1000, 10000, 1000)
nthreads_l = np.arange(50, 90, 5)

r = []

coldrun = True

for rows, nthreads in itertools.product(rows_l, nthreads_l):
    # create the argument list
    args = [EXEC_PATH] + args_t.format(rows=rows).split(' ')

    # set the number of threads
    env['OMP_NUM_THREADS'] = str(nthreads)

    if coldrun:
        p = subprocess.Popen(args, stdout=subprocess.PIPE, env=env)
        lines = p.communicate()
        coldrun = False

    # execute command
    t = 0
    for _ in range(N_RUNS):
        p = subprocess.Popen(args, stdout=subprocess.PIPE, env=env)
        lines = p.communicate()[0].decode('utf-8').split('\n')
    
        t += float(next(filter(lambda e: e.startswith('Time: '), lines))[6:])

    t /= N_RUNS

    print(rows, nthreads, t)
    r.append(t)
    
r = np.array(r).reshape(len(nthreads_l), len(rows_l))

fig = plt.figure()
ax = plt.axes(projection='3d')

X, Y = np.meshgrid(rows_l, nthreads_l)
ax.plot_surface(X, Y, r, cmap=cm.coolwarm)
ax.xaxis.set_label('number of rows')
ax.yaxis.set_label('number of threads')

plt.show()