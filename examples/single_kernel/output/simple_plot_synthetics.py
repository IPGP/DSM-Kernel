#!/usr/bin/env python
"""
very simple visualization script for the DSM-Kernel example. For more advanced
visualizations, check out the folder /util/kernelviewer/
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('eq.Explosion.Zs.dat.100s10s')
time = data[:,0]
raw = data[:,1]
filtered = data[:,2]
plt.plot(time, raw)
plt.plot(time, filtered)
plt.show()
