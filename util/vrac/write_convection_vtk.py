#!/Users/fujinobuaki/anaconda/bin/python
##!/tools/python


import matplotlib.pyplot as plt
import numpy as np
import sys
from pyevtk.hl import gridToVTK



def read_fortran_record(binfile, count, dtype, filetype):
    """reads a sequential fortran binary file record"""
    if filetype=='sequential':
        rec_start = np.fromfile(binfile, count=1, dtype=np.int32)
        content = np.fromfile(binfile, count=count, dtype=dtype)
        rec_end = np.fromfile(binfile, count=1, dtype=np.int32)
        assert rec_start == rec_end, 'record incomplete. wrong count or datatype'
    else:
        content = np.fromfile(binfile, count=count, dtype=dtype)
   
    return content


fname_kernel=sys.argv[1]
fname_vtk=sys.argv[2]
nktype=int(sys.argv[3])
radius_depth_type=sys.argv[4]







# read grid info:
nr=96
nphi=382
ntheta=192
print(nktype)

#radii=np.linspace(3479.5,6371.0,nr+1)
#local_phis=np.linspace(-180.0, 180.0, nphi+1)
#local_thetas=np.linspace(-90.0, 90.0, ntheta+1)
#phis=np.zeros(nphi*ntheta)
#thetas=np.zeros(nphi*ntheta
kernelfile = open(fname_kernel, 'rb')
xgrid=np.zeros((ntheta,nphi,nr))
ygrid=np.zeros((ntheta,nphi,nr))
zgrid=np.zeros((ntheta,nphi,nr))
kernel=np.zeros((nktype,ntheta,nphi,nr))
ipoint=0

lines=kernelfile.readlines()
for ir in range (0,nr):
    for itheta in range (0,ntheta):
        for iphi in range (0,nphi):
            #tmp=read_fortran_record(kernelfile,count=3+nktype,dtype=np.float32,filetype='')
            tmp=np.fromstring(lines[ipoint],dtype=float,sep=' ')
            #tmp=np.loadtxt(kernelfile, delimiter=' ')
            #print(tmp)
            #exit()
            rEarth=6371.0
            # if radius is input
            if(radius_depth_type=='radius'): radiusTmp=tmp[2]
            # if depth is input
            if(radius_depth_type=='depth'): radiusTmp=rEarth-tmp[2]
            zgrid[itheta,iphi,ir]=radiusTmp*np.sin(np.radians(tmp[1]))
            xgrid[itheta,iphi,ir]=radiusTmp*np.cos(np.radians(tmp[1]))*np.cos(np.radians(tmp[0]))
            ygrid[itheta,iphi,ir]=radiusTmp*np.cos(np.radians(tmp[1]))*np.sin(np.radians(tmp[0]))
            #vs[ipoint]=tmp[3]
            #vp[ipoint]=tmp[4]
            #vbulk[ipoint]=tmp[5]
            #rho[ipoint]=tmp[6]
            for iktype in range(0,nktype):
                kernel[iktype,itheta,iphi,ir]=tmp[2+iktype]
            ipoint=ipoint+1


        

kernelfile.close()
     
    # write vtk file
    #xgrid = np.outer(np.sin(np.radians(thetas)) * np.cos(np.radians(phis)), radii)
    #ygrid = np.outer(np.sin(np.radians(thetas)) * np.sin(np.radians(phis)), radii)
    #zgrid = np.outer(np.cos(np.radians(thetas)), radii)
     
xgrid = xgrid.reshape(ntheta, nphi, nr)
ygrid = ygrid.reshape(ntheta, nphi, nr)
zgrid = zgrid.reshape(ntheta, nphi, nr)
    
point_data = {'{:d}'.format(name): data for name, data in zip(range(nktype), kernel)}
     
gridToVTK(fname_vtk, xgrid, ygrid, zgrid, pointData=point_data)
     
    #phis = phis.reshape(ntheta, nphi)
    #thetas = phis.reshape(ntheta, nphi)
    

exit()
# write some output information
r_min, r_max = radii[0], radii[-1]
phi_min, phi_max = phis[0, 0], phis[0, -1]
theta_min, theta_max = thetas[0, 0], thetas[-1, 0]
print ('kernel dimensions (nr={}, nphi={}, ntheta={})'.format(nr, nphi, ntheta))
print ('kernel types: {}'.format(nktype))
print ('radius range = {} -> {}'.format(r_min, r_max))
print ('phis range = {} -> {}'.format(phi_min, phi_max))
print ('thetas range = {} -> {}'.format(theta_min, theta_max))

