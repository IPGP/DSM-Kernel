#!/usr/bin/env python
"""
very simple visualization script for the DSM-Kernel example. For more advanced
visualizations, check out the folder /util/kernelviewer/
"""

import numpy as np
import matplotlib.pyplot as plt

fname_synthetics = 'output/eq.Explosion.Zs.dat.100s10s'
data = np.loadtxt(fname_synthetics)
time = data[:,0]
raw = data[:,1]
filtered = data[:,2]

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
ax.plot(time, filtered, label='filtered synthetics')
ax.set(xlabel='time [s]', ylabel='displacement [?]')
ax.legend()

fig.savefig('synthetics.png')
