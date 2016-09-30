#!/usr/bin/env python
"""
This is a very simple plot script that reads a 3D DSM-Kernel sensitivity kernel
file and make a plot of it.
"""

import matplotlib.pyplot as plt
import numpy as np

def read_fortran_record(binfile, count, dtype):
    """reads a sequential fortran binary file record"""
    rec_start = np.fromfile(binfile, count=1, dtype=np.int32)
    content = np.fromfile(binfile, count=count, dtype=dtype)
    rec_end = np.fromfile(binfile, count=1, dtype=np.int32)
    assert rec_start == rec_end, 'record incomplete. wrong count or datatype'
    return content

fname_kernel = 'output/eq.Explosion.Z.Z.100s10s.kernel'
fname_grid = 'output/eq.Explosion.Z.Z.grid'
fname_plot = 'kernel.png'

kernelfile = open(fname_kernel, 'rb')
gridfile = open(fname_grid, 'rb')

# read grid info:
nr, nphi, ntheta, nktype = read_fortran_record(gridfile, count=4, dtype=np.int32)
nktype += 1 # starts counting from zero (should be changed in the code?)
radii = read_fortran_record(gridfile, count=nr, dtype=np.float32)
phis = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32)
thetas = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32)
phis = phis.reshape(ntheta, nphi)
thetas = phis.reshape(ntheta, nphi)

# read kernel data:
npoints = nr * nphi * ntheta
kernel = read_fortran_record(kernelfile, count=npoints * nktype, dtype=np.float32)
kernel = kernel.reshape(nktype, ntheta, nphi, nr)

# write some output information
r_min, r_max = radii[0], radii[-1]
phi_min, phi_max = phis[0, 0], phis[0, -1]
theta_min, theta_max = thetas[0, 0], thetas[-1, 0]
print 'kernel dimensions (nr={}, nphi={}, ntheta={})'.format(nr, nphi, ntheta)
print 'kernel types: {}'.format(nktype)
print 'radius range = {} -> {}'.format(r_min, r_max)
print 'phis range = {} -> {}'.format(phi_min, phi_max)
print 'thetas range = {} -> {}'.format(theta_min, theta_max)

# plot a slice of the kernel:
ktype_labels = {1: 'vp phase',
                3: 'vp amplitude',
                4: r'$Q_\kappa$ amplitude',
                5: 'vs phase',
                7: 'vs amplitude',
                8: r'$Q_\mu$ amplitude'}
itheta = ntheta/2 # plot the greatcircle plane

ncols = 3
irow_max = (nktype - 1) / ncols
nrows = irow_max + 1

fig, axes = plt.subplots(nrows, ncols, figsize=(10, 5), sharex=True, sharey=True)
extent = (phi_min, phi_max, r_min, r_max)
for iktype, ax in zip(range(nktype), axes.flatten()):
    irow = iktype / ncols
    icol = iktype % ncols

    kernel_slice = kernel[iktype, itheta, :, :].T
    std = np.std(kernel_slice)
    ax.imshow(kernel_slice, vmin=-std, vmax=std, aspect='auto', origin='bottom',
              extent=extent)
    ax.set(title=ktype_labels.get(iktype, ''))
    if irow == irow_max:
        ax.set_xlabel('distance [arc degrees]')
    if icol == 0:
        ax.set_ylabel('radius [km]')

fig.tight_layout(rect=[0.01, 0.03, 0.95, 0.99])
fig.savefig(fname_plot)

try:
    plt.show()
except Exception:
    print 'plot in {}'.format(fname_plot)
