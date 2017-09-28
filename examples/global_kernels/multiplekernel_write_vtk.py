#!/cm/shared/apps/intelpython/intelpython3/bin/python
##!/Users/fujinobuaki/anaconda/bin/python
##!/tools/python
"""
This is a very simple plot script that reads a 3D DSM-Kernel sensitivity kernel
file and make a plot of it.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
from pyevtk.hl import gridToVTK


head=("head","<i")
tail=("tail","<i")


def read_kernel_inffile(inffile):
    """reads kernel.inf files and generate all the variables needed"""
    iline=0
    kernelType=''
    Butterworth=''
    fastFFT=''
    timeWindow=''
    timeStart=''
    timeEnd=''
    timeIncrementSec=''
    for line in open(inffile):
        li=line.strip()
        if not li.startswith("#"):
            print(iline)
            tmpstring=line.rstrip()
            #print(line.rstrip())
            print(tmpstring)
            if(iline==1): outputDir=tmpstring
            if(iline==2): eventName=tmpstring
            if(iline==5): stationName=tmpstring
            if(iline==7): phaseName=tmpstring
            if(iline==8): component=tmpstring
            if(iline==9): kernelType=tmpstring
            if(any (x not in kernelType for x in ('V' and 'v'))):
                if(iline==11): filterName=tmpstring
                if(iline==9): timeIncrementSec=tmpstring
                if(iline==10): Butterworth=tmpstring
                if(iline==11): filterName=tmpstring
                if(Butterworth=='1'):
                    if(iline==17): fastFFT=tmpstring
                    if(fastFFT=='1'):
                        if(iline==22): timeWindow=tmpstring
                    else:
                        if(iline==21): timeWindow=tmpstring
                else:
                    if(iline==16): fastFFT=tmpstring
                    if(fastFFT=='1'):
                        if(iline==21): timeWindow=tmpstring
                    else:
                        if(iline==20): timeWindow=tmpstring


            else:
                if(iline==10): timeIncrementSec=tmpstring
                if(iline==11): Butterworth=tmpstring
                if(iline==12): filterName=tmpstring
                if(Butterworth=='1'): 
                    if(iline==18): fastFFT=tmpstring
                    if(fastFFT=='1'):
                        if(iline==23): timeWindow=tmpstring
                    else:
                        if(iline==22): timeWindow=tmpstring
                else:
                    if(iline==17): fastFFT=tmpstring
                    if(fastFFT=='1'):
                        if(iline==22): timeWindow=tmpstring
                    else:
                        if(iline==21): timeWindow=tmpstring
                       
                                  
            iline=iline+1 
    print(timeWindow) 
    print (kernelType)
    timeStart=timeWindow.split(" ")[0]
    timeEnd=timeWindow.split(" ")[1]
    return {'outputDir':outputDir,'eventName':eventName,'stationName':stationName,'phaseName':phaseName,'component':component,'kernelType':kernelType,'filterName':filterName,'timeStart':timeStart,'timeEnd':timeEnd,'timeIncrementSec':timeIncrementSec}

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


kernel_inf=read_kernel_inffile(sys.argv[1])

fname_grid=kernel_inf['outputDir']+'/'+kernel_inf['stationName']+'.' \
    +kernel_inf['eventName']+'.'+kernel_inf['phaseName']+'.'\
    +kernel_inf['component']+'.'+kernel_inf['kernelType']+'.grid'


if any (x not in kernel_inf['kernelType'] for x in ('V' and 'v')):
    print("this is not a video")
    fname_kernel=kernel_inf['outputDir']+'/'+kernel_inf['stationName']+'.' \
    +kernel_inf['eventName']+'.'+kernel_inf['phaseName']+'.'\
    +kernel_inf['component']+'.'+kernel_inf['filterName']+'.kernel'
    print(fname_kernel)

    kernelfile = open(fname_kernel, 'rb')
    gridfile = open(fname_grid, 'rb')
    
    fname_vtk=kernel_inf['stationName']+'.' \
    +kernel_inf['eventName']+'.'+kernel_inf['phaseName']+'.'\
    +kernel_inf['component']+'.'+kernel_inf['filterName']

     # read grid info:
    nr, nphi, ntheta, nktype = read_fortran_record(gridfile, count=4, dtype=np.int32,filetype='sequential')
    nktype += 1 # starts counting from zero (should be changed in the code?)
    radii = read_fortran_record(gridfile, count=nr, dtype=np.float32,filetype='sequential')
    phis = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32,filetype='sequential')
    thetas = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32,filetype='sequential')
    gridfile.close()
     
    # read kernel
    npoints = nr * nphi * ntheta
    kernel = read_fortran_record(kernelfile, count=npoints * nktype, dtype=np.float32,filetype='direct')
    kernel = kernel.reshape(nktype, ntheta, nphi, nr)
    kernelfile.close()
     
    # write vtk file
    xgrid = np.outer(np.sin(np.radians(thetas)) * np.cos(np.radians(phis)), radii)
    ygrid = np.outer(np.sin(np.radians(thetas)) * np.sin(np.radians(phis)), radii)
    zgrid = np.outer(np.cos(np.radians(thetas)), radii)
     
    xgrid = xgrid.reshape(ntheta, nphi, nr)
    ygrid = ygrid.reshape(ntheta, nphi, nr)
    zgrid = zgrid.reshape(ntheta, nphi, nr)
    
    point_data = {'{:d}'.format(name): data for name, data in zip(range(nktype), kernel)}
     
    gridToVTK(fname_vtk, xgrid, ygrid, zgrid, pointData=point_data)
     
    phis = phis.reshape(ntheta, nphi)
    thetas = phis.reshape(ntheta, nphi)
    
    

else:

    print("this is a video")
    print(kernel_inf['timeStart'])
    print(kernel_inf['timeEnd'])
    numSnaps=int((float(kernel_inf['timeEnd'].replace("d","e"))-float(kernel_inf['timeStart'].replace("d","e")))/float(kernel_inf['timeIncrementSec'].replace("d","e")))+1
    print(numSnaps)

    # il me faut calculer le nombre de snapshots pour video mode
    for itime in range (1,numSnaps+1):
        
        num_snap=str(itime).zfill(7)
        fname_kernel=kernel_inf['outputDir']+'/tmpvideo/'\
            +kernel_inf['stationName']+'.' \
            +kernel_inf['eventName']+'.'+kernel_inf['phaseName']+'.'\
            +kernel_inf['component']+'.'+kernel_inf['filterName']+'.'\
            +num_snap+'.video'
        #print(fname_kernel)
   
     

        fname_vtk=kernel_inf['stationName']+'.' \
            +kernel_inf['eventName']+'.'+kernel_inf['phaseName']+'.'\
            +kernel_inf['component']+'.'+kernel_inf['filterName']+'.'\
            +num_snap
        #print(fname_vtk)
       
        kernelfile = open(fname_kernel, 'rb')
        gridfile = open(fname_grid, 'rb')
        
        # read grid info:
        nr, nphi, ntheta, nktype = read_fortran_record(gridfile, count=4, dtype=np.int32,filetype='sequential')
        nktype += 1 # starts counting from zero (should be changed in the code?)
        radii = read_fortran_record(gridfile, count=nr, dtype=np.float32,filetype='sequential')
        phis = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32,filetype='sequential')
        thetas = read_fortran_record(gridfile, count=nphi * ntheta, dtype=np.float32,filetype='sequential')
        gridfile.close()
        
    # read kernel
        npoints = nr * nphi * ntheta
        kernel = read_fortran_record(kernelfile, count=npoints * nktype, dtype=np.float32,filetype='direct')
        kernel = kernel.reshape(nktype, ntheta, nphi, nr)
        kernelfile.close()
        
    # write vtk file
        xgrid = np.outer(np.sin(np.radians(thetas)) * np.cos(np.radians(phis)), radii)
        ygrid = np.outer(np.sin(np.radians(thetas)) * np.sin(np.radians(phis)), radii)
        zgrid = np.outer(np.cos(np.radians(thetas)), radii)
        
        xgrid = xgrid.reshape(ntheta, nphi, nr)
        ygrid = ygrid.reshape(ntheta, nphi, nr)
        zgrid = zgrid.reshape(ntheta, nphi, nr)
        
        point_data = {'{:d}'.format(name): data for name, data in zip(range(nktype), kernel)}
        
        gridToVTK(fname_vtk, xgrid, ygrid, zgrid, pointData=point_data)
        
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

