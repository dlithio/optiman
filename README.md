Generates stable and unstable manifolds

# Installation
## Prerequisites
You should have gfortran, python, numpy, scipy, and mayavi installed to run the program that finds the manifolds. There is also a program included to help find the fixed point and the unstable eigenvectors. To run this program you should also have LAPACK.

### Linux (Ubuntu or Debian)
```
sudo apt-get install gfortran liblapack3 liblapack-dev liblapack-doc python python-numpy python-scipy mayavi2
```

### Mac
* Install gfortran from https://gcc.gnu.org/wiki/GFortranBinaries. Be sure to read the detailed instructions, there are several steps you need to take before installing the .dmg.
* Install Enthought Canopy Express from https://store.enthought.com/downloads/. 
* Install LAPACK. It may come preinstalled for Mac 10.9 and above, but I am unable to test it. The following site and its comments could be helpful - https://pheiter.wordpress.com/2012/09/04/howto-installing-lapack-and-blas-on-mac-os/.

### Windows 
* Install gfortran from https://gcc.gnu.org/wiki/GFortranBinaries. Be sure to read the detailed instructions, there are several steps you need to take before installing the .dmg.
* Install Enthought Canopy Express from https://store.enthought.com/downloads/. 
* Install LAPACK. This may get a bit tricky. It'd probably be easier to find a Linux machine or install a virtual machine/dual boot ubuntu on your current system.

## Program installation
If you have git installed, simply clone the repository by 
```
git clone https://github.com/dlithio/optiman.git
```
Alternatively, you can also use the "Download Zip" button on the right side of the page and then manually unzip the folder. This could be done on the command line by 
```
unzip optiman-master.zip
```

# Computing Your First Manifold

There are several demos included in the demos folder. The demos lorenz, kse_origin, and kse_bottom_bimodal will likely be the closest to what you want to do and they are all done the same way.

1. Move to the base folder of the project and copy the contents of the demos/lorenz folder into it.
```
~/optiman(master ✗) ls
demos  LICENSE  makefile  obj  README.md  results  src
~/optiman(master ✗) cp demos/lorenz/* .
~/optiman(master ✗) ls
demos  fixed_point_input  initial_guess  LICENSE  lorenz.f90  makefile  obj  optiman_input  par  README.md  results  src
~/optiman(master ✗) 
```
