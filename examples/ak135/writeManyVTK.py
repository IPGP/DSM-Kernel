#!/usr/bin/env python
"""
This is a very simple plot script that reads a 3D DSM-Kernel sensitivity kernel
file and make a plot of it.
"""

import matplotlib.pyplot as plt
import numpy as np
from pyevtk.hl import gridToVTK

def read_fortran_record(binfile, count, dtype):
    """reads a sequential fortran binary file record"""
    rec_start = np.fromfile(binfile, count=1, dtype=np.int32)
    content = np.fromfile(binfile, count=count, dtype=dtype)
    rec_end = np.fromfile(binfile, count=1, dtype=np.int32)
    assert rec_start == rec_end, 'record incomplete. wrong count or datatype'
    return content

#fname_kernel='output/test90.TP.S.T.200s20s.kernel'
#fname_grid='output/test90.TP.S.T.beta.grid'

def test(fA,fB):
    print(fA,fB)
def write_vts2(kernel,fname_grid,fname_vts):
    # read grid info:
    nr, nphi, ntheta, nktype = read_fortran_record(gridfile, count=4, dtype=np.int32)
    nktype += 1 # starts counting from zero (should be changed in the code?)
    radii = read_fortran_record(gridfile, count=nr, dtype=np.float32)
    phis = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32)
    thetas = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32)
    # write vtk file
    xgrid = np.outer(np.sin(np.radians(thetas)) * np.cos(np.radians(phis)), radii)
    ygrid = np.outer(np.sin(np.radians(thetas)) * np.sin(np.radians(phis)), radii)
    zgrid = np.outer(np.cos(np.radians(thetas)), radii)

    xgrid = xgrid.reshape(ntheta, nphi, nr)
    ygrid = ygrid.reshape(ntheta, nphi, nr)
    zgrid = zgrid.reshape(ntheta, nphi, nr)

    point_data = {'s'+'{:d}'.format(name): data for name, data in zip(range(nktype), kernel)}

    gridToVTK(fname_vts, xgrid, ygrid, zgrid, pointData=point_data)

    phis = phis.reshape(ntheta, nphi)
    thetas = phis.reshape(ntheta, nphi)

    # read kernel data:

    # write some output information
    r_min, r_max = radii[0], radii[-1]
    phi_min, phi_max = phis[0, 0], phis[0, -1]
    theta_min, theta_max = thetas[0, 0], thetas[-1, 0]
    print ('kernel dimensions (nr={}, nphi={}, ntheta={})'.format(nr, nphi, ntheta))
    print ('kernel types: {}'.format(nktype))
    print ('radius range = {} -> {}'.format(r_min, r_max))
    print ('phis range = {} -> {}'.format(phi_min, phi_max))
    print ('thetas range = {} -> {}'.format(theta_min, theta_max))

def write_vts(fname_kernel,fname_grid,fname_vts):
#fname_kernel='output/test90.TP.SS.T.200s20s.kernel'
#fname_grid='output/test90.TP.SS.T.beta.grid'
#fname_vts='SS90'
#fname_plot='SSkernel90.png'
    kernelfile = open(fname_kernel, 'rb')
    gridfile = open(fname_grid, 'rb')

    # read grid info:
    nr, nphi, ntheta, nktype = read_fortran_record(gridfile, count=4, dtype=np.int32)
    nktype += 1 # starts counting from zero (should be changed in the code?)
    radii = read_fortran_record(gridfile, count=nr, dtype=np.float32)
    phis = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32)
    thetas = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32)

    # read kernel
    npoints = nr * nphi * ntheta
    kernel = read_fortran_record(kernelfile, count=npoints * nktype, dtype=np.float32)
    kernel = kernel.reshape(nktype, ntheta, nphi, nr)
    kernel = kernel*1.e9

    # write vtk file
    xgrid = np.outer(np.sin(np.radians(thetas)) * np.cos(np.radians(phis)), radii)
    ygrid = np.outer(np.sin(np.radians(thetas)) * np.sin(np.radians(phis)), radii)
    zgrid = np.outer(np.cos(np.radians(thetas)), radii)

    xgrid = xgrid.reshape(ntheta, nphi, nr)
    ygrid = ygrid.reshape(ntheta, nphi, nr)
    zgrid = zgrid.reshape(ntheta, nphi, nr)

    point_data = {'s'+'{:d}'.format(name): data for name, data in zip(range(nktype), kernel)}

    gridToVTK(fname_vts, xgrid, ygrid, zgrid, pointData=point_data)

    phis = phis.reshape(ntheta, nphi)
    thetas = phis.reshape(ntheta, nphi)

    # read kernel data:
    
    # write some output information
    r_min, r_max = radii[0], radii[-1]
    phi_min, phi_max = phis[0, 0], phis[0, -1]
    theta_min, theta_max = thetas[0, 0], thetas[-1, 0]
    print ('kernel dimensions (nr={}, nphi={}, ntheta={})'.format(nr, nphi, ntheta))
    print ('kernel types: {}'.format(nktype))
    print ('radius range = {} -> {}'.format(r_min, r_max))
    print ('phis range = {} -> {}'.format(phi_min, phi_max))
    print ('thetas range = {} -> {}'.format(theta_min, theta_max))
    return kernel

f0='/home/nobuaki/scratch/Earth/ak135/kernel/'
distance=np.arange(50,100,10).tolist()
#print(distance)
for idistance in distance:
    Skernel=f0+'test'+'{:d}'.format(idistance)+'.TP.S.T.100s20s.kernel'
    Sgrid=f0+'test'+'{:d}'.format(idistance)+'.TP.S.T.beta.grid'
    Svtk=f0+'{:d}'.format(idistance)+'.S.vts'
    SSkernel=f0+'test'+'{:d}'.format(idistance)+'.TP.SS.T.100s20s.kernel' 
    SSgrid=f0+'test'+'{:d}'.format(idistance)+'.TP.SS.T.beta.grid'   
    SSvtk=f0+'{:d}'.format(idistance)+'.SS.vts'
    diffvtk=f0+'{:d}'.format(idistance)+'.diff.vts'
    #print(diffvtk)
    kernelS=write_vts(Skernel,Sgrid,Svtk)
    kernelSS=write_vts(SSkernel,SSgrid,SSvtk)
    kernelD=kernelSS-kernelS
    write_vts2(kernelD,Sgrid,diffvtk) 
