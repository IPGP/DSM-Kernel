<img src="etc/images/DSM_logo.png"  width="250" align="middle">

### Welcome to DSM Kernel Suite!
This suite of programmes is to calculate 3D finite frequency Fréchet sensitivity kernels (or 3D waveform partial derivatives) for 1D reference Earth models using Direct Solution Method. It consists of forward and back propagated strain Green's function calculation in a 2D plane, cross correlations of forward and back propagated wavefields to obtain sensitivity kernels. We developed also visualisation facilities of kernels. Those who would like to get kernels by comparing observed and synthetic waveforms, we developed python interface for that purpose as well. All the main motors are parallelised and you just have to submit prepared scripts. 

## Programmes inside DSM Kernel Suite
DSM Kernel Suite have five independent programmes for each step. If you fetch the application's source with git, this dependency will be fetched automatically. 
* SGTpsv : Green's function calculation for PSV mode (MPI fortran, NF) 
* SGTsh  : Green's function calculation for SH  mode (MPI fortran, NF)
* wave2kernel : visualisation of SAC files in order to compare them with synthetics and to decide which phase (python, HJ)
* KernelMaker : 3D Fréchet derivatives calculation (MPI fortran, NF)
* KernelViewer : 3D Kernel visualisation facilities (python, MM)
* ModelDrawer : 1D model ASCII file generator from polynomial model file (NF)

recommended requirements:
* ifort

optional requirements:
* obspy
* VTK

## Build and Install
* Precise compiler options using `etc/config_calcul/config.h` (compiler options etc.) 
* Now you just have to generate executables by `./configure; make'

## Examples
 * [single station and event](examples/single_kernel/README.md)
 * [multiple stations and events](examples/multiple_kernels/README.md)

## Kernel Gallery
<img src="etc/images/kernel1.png" width="300">

## Papers to be cited
Fuji, N., Chevrot, S., Zhao, L., Geller, R.J., Kawai, K. (2012) [Finite-frequency structural sensitivities of short-period compressional body waves](https://gji.oxfordjournals.org/content/190/1/522.full), Geophys. J. Int., 190, 522-540.

## Authors and Contributors
All rights reserved, 2016, Nobuaki Fuji (@seismobassoon), Matthias Meschede
(@MMesch), Kensuke Konishi (@kensuke1984), Hugo Jaegler (@hJaegler), Kenji
Kawai, Li Zhao, Sebastien Chevrot, Robert J. Geller, Vadim Monteiller, Dimitri
Komatitsch, Marie Calvet, Hiromitsu Mizutani

## Support or Contact
Having trouble with DSM Kernel Suite? Check out our
[documentation](http://ipgp.github.io/DSM-Kernel/) or [contact
support](email:nobuaki@ipgp.fr) and we will help you sort it out.

## Legal info
DSM Kernel Suite is Free/Libre/Open Source software and available under the
GPLv3 license ([License](LICENSE.txt)). A list of contributors can be found in
the about dialog.
