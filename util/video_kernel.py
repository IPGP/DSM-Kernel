#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import ipdb
import mydebug

def read_record(binfile,dtype='<i4'):
    nbytes1 = np.fromfile(binfile,'<i4',1)[0]
    record = np.fromfile(binfile,dtype,nbytes1/int(dtype[2]))
    nbytes2 = np.fromfile(binfile,'<i4',1)[0]
    assert (nbytes1==nbytes2),'record incomplete'
    return record

class KernelSeries(object):
    def __init__(s,gridfile,kernelfiles):
        s.gridfile = gridfile
        s.kernelfiles = kernelfiles
        s.nkernels = len(s.kernelfiles)

        infile = open(gridfile,'rb')
        
        s.nr,s.nphi,s.ntheta,\
        s.nkvtype,s.nfilt,s.nt_seismo,\
        s.nt_video = read_record(infile,dtype='<i4')
        s.tstart,s.tend = read_record(infile,dtype='<f4')
        s.mtensor  = read_record(infile,dtype='<f4')
        s.radii    = read_record(infile,dtype='<f4')
        s.phis     = read_record(infile,dtype='<f4').reshape(s.ntheta,s.nphi)
        s.thetas   = read_record(infile,dtype='<f4').reshape(s.ntheta,s.nphi)
        s.seismo   = read_record(infile,dtype='<f4')

        s.ntot = s.nr*s.ntheta*s.nphi

        print 'reading grid with dimension:',s.nr,s.ntheta,s.nphi

    def get_kernel(s,itime):
        fname = s.kernelfiles[itime]
        binfile = open(fname,'rb')
        kernel = np.fromfile(binfile,'<f4',s.ntot).reshape(s.nr,s.ntheta,s.nphi)
        binfile.close()
        return kernel

    def plot_seismo(s):
        plt.figure()
        plt.plot(s.seismo)

    def plot_kernel(s,itime):
        igreatcircle = s.ntheta/2
        kernel = np.abs(s.get_kernel(itime)[:,igreatcircle,:])
        kmax = np.max(kernel)
        plt.figure()
        norm = mpl.colors.LogNorm(kmax*1e-7,kmax)
        im = plt.imshow(kernel,norm=norm,origin='lower')
        plt.colorbar(im)

    def writevtk(s, folder='./'):
        from evtk.hl import gridToVTK
 
        xgrid = np.outer(np.sin(np.radians(s.thetas))*np.cos(np.radians(s.phis)),s.radii)
        ygrid = np.outer(np.sin(np.radians(s.thetas))*np.sin(np.radians(s.phis)),s.radii)
        zgrid = np.outer(np.cos(np.radians(s.thetas)),s.radii)

        xgrid = xgrid.reshape(s.ntheta,s.nphi,s.nr)
        ygrid = ygrid.reshape(s.ntheta,s.nphi,s.nr)
        zgrid = zgrid.reshape(s.ntheta,s.nphi,s.nr)

        for ikernel in range(s.nkernels):
            print 'writing kernel {:d}'.format(ikernel)
            kernel = np.transpose(s.get_kernel(ikernel),(1,2,0)).copy()
            gridToVTK(folder+'kernel_{:03d}'.format(ikernel),xgrid,ygrid,zgrid,pointData={'kernel':kernel})

#===== MAIN FUNCTION ====
def main():
    print 'main function'
    gridfile = sys.argv[1]
    fnames   = sorted(sys.argv[2:])
    print fnames
    kernels = KernelSeries(gridfile,fnames)
    kernels.writevtk(folder='vtk/')
    #plt.show()

#==== SCRIPT EXECUTION ====
if __name__ == "__main__":
    main()
