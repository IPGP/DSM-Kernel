#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import ipdb

class DSM_Kernel(object):
    def __init__(s,fname):
        """
        reads combined kernel file from Nobuakis KernelMaker code
        """
        infile = open(fname)

        #---- read header and save some values in the class ----
        headerstruct = np.dtype([('rec1s','<i4'),('ndims','<i4',3),('nwin','<i4'),('kc','<i4'),
                                 ('mtype','<i4'),('nktype','<i4'),('sym','<f4'),('nfilter','<i4'),
                                 ('ifilter','<i4'),('rec1e','<i4'),

                                 ('rec2s','<i4'),('dtn','<f4'),('twin','<f4',4),('trange','<f4',2),
                                 ('rec2e','<i4'),

                                 ('rec3s','<i4'),('rec3_data','<f4',5),('rec3e','<i4'),

                                 ('rec4s','<i4'),('latlondep_so','<f4',3),('sdep','<f4',1), #?
                                 ('latlon_re','<f4',2),('distance','<f4',1),('rec4e','<i4'),
                                 
                                 ('rec5s','<i4'),('rec5_data','<f4',6),('rec5e','<i4'),

                                 ('rec6s','<i4'),('moment_tensor','<f4',6),('rec6e','<i4'),

                                 ('rec7s','<i4'),('r0D','<f4',1),('rec7e','<i4')
                                 ])
        header = np.fromfile(infile,headerstruct,1)[0]
        for i in range(7):
            print 'record {:d}:'.format(i+1),
            print header['rec{:d}s'.format(i+1)],
            print header['rec{:d}e'.format(i+1)]

        s.nr,s.nphi,s.ntheta = header['ndims']
        s.nwin               = header['nwin']
        s.nktype             = header['nktype']
        s.nfilter            = header['nfilter']
        s.trange             = header['trange']
        s.slat,s.slon,s.sdep = header['latlondep_so']
        s.rlat,s.rlon        = header['latlon_re']
        s.distance           = header['distance']
        s.moment_tensor      = header['moment_tensor']
        s.r0D                = header['r0D']

        print 'grid dimension:',s.nr,s.nphi,s.ntheta
        print 'number of time samples:',s.nwin
        print 'window range:',s.trange

        #---- read data arrays ----
        nbytes1 = np.fromfile(infile,'<i4',1)[0]
        s.radii = np.fromfile(infile,'<f4',s.nr)
        nbytes2 = np.fromfile(infile,'<i4',1)[0]
        print 'record 8:',nbytes1,nbytes2

        nbytes1    = np.fromfile(infile,'<i4',1)[0]
        s.phitheta = np.fromfile(infile,'<f4',nbytes1/4)
        nbytes2    = np.fromfile(infile,'<i4',1)[0]
        print 'record 9:',nbytes1,nbytes2

        nbytes1    = np.fromfile(infile,'<i4',1)[0]
        s.thetaphi = np.fromfile(infile,'<f4',nbytes1/4)
        nbytes2    = np.fromfile(infile,'<i4',1)[0]
        print 'record 10:',nbytes1,nbytes2

        nbytes1  = np.fromfile(infile,'<i4',1)[0]
        s.kernel = np.fromfile(infile,'<f4',s.nr*s.nphi*s.ntheta*(s.nktype+1))
        s.kernel = s.kernel.reshape(s.nktype+1,s.ntheta,s.nphi,s.nr)
        nbytes2  = np.fromfile(infile,'<i4',1)[0]
        print 'record 11:',nbytes1,nbytes2
        for itype in range(s.nktype):
            print np.sum(s.kernel[itype])

        nbytes1  = np.fromfile(infile,'<i4',1)[0]
        s.displu = np.fromfile(infile,'<f4',s.nwin)
        nbytes2  = np.fromfile(infile,'<i4',1)[0]
        print 'record 12:',nbytes1,nbytes2

        #nbytes1  = np.fromfile(infile,'<i4',1)[0]
        #s.displv = np.fromfile(infile,'<f4',s.nwin)
        #nbytes2  = np.fromfile(infile,'<i4',1)[0]
        #print 'record 13:',nbytes1,nbytes2

        #nbytes1  = np.fromfile(infile,'<i4',1)[0]
        #s.displh = np.fromfile(infile,'<f4',s.nwin)
        #nbytes2  = np.fromfile(infile,'<i4',1)[0]
        #print 'record 14:',nbytes1,nbytes2

        infile.close()

    def plot_seismo(s):
        plt.figure()
        times = np.linspace(s.trange[0],s.trange[1],s.nwin)
        plt.plot(times,s.displu)
        #plt.plot(times,s.displv)
        #plt.plot(times,s.displh)

    def plot_slices(s,itype):
        limit = np.mean(np.abs(s.kernel[itype]))
        norm=plt.Normalize(-limit,limit)
        plt.figure()
        plt.imshow(np.transpose(s.kernel[itype,s.ntheta/2,:,:]),norm=norm,origin='bottom')

    def write_vtk(s):
        from evtk.hl import gridToVTK

        xgrid = np.outer(np.sin(np.radians(s.thetaphi))*np.cos(np.radians(s.phitheta)),s.radii)
        ygrid = np.outer(np.sin(np.radians(s.thetaphi))*np.sin(np.radians(s.phitheta)),s.radii)
        zgrid = np.outer(np.cos(np.radians(s.thetaphi)),s.radii)

        xgrid = xgrid.reshape(s.ntheta,s.nphi,s.nr)
        ygrid = ygrid.reshape(s.ntheta,s.nphi,s.nr)
        zgrid = zgrid.reshape(s.ntheta,s.nphi,s.nr)

        print s.kernel[5].shape
        print xgrid.shape
        print ygrid.shape
        print zgrid.shape

        gridToVTK('kernel',xgrid,ygrid,zgrid,pointData={'kernel':s.kernel[8,:,:,:]})

def main():
    fname = sys.argv[1]
    kernel = DSM_Kernel(fname)
    kernel.plot_seismo()
    kernel.plot_slices(7)
    kernel.write_vtk()
    plt.show()

if __name__ == "__main__":
    main()
