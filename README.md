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

 ```bash
~/optiman(master ✗) ls
demos  LICENSE  makefile  obj  README.md  results  src
~/optiman(master ✗) cp demos/lorenz/* .
~/optiman(master ✗) ls
demos  fixed_point_input  initial_guess  LICENSE  lorenz.f90  makefile  obj  optiman_input  par  README.md  results  src
~/optiman(master ✗) 
```

 This copies 4 important files into the directory. These are 4 files that you'll need to include with any manifold you try to compute.
 * **fixed_point_input** - Contains some parameters for finding the fixed point and eigenvectors of the fixed point that you'll be plotting the Manifold of. The newton iterations will not stop until the norm of the vector field is below *tolerance* or the program has completed *max_iterations* Newton steps.
 * **initial_guess** - Contains the initial guess for the fixed point as a text file. If you open and view this file, you'll see that the initial guess is (1,1,1). We could have specified (0,0,0) since we happen to know that's the exact fixed point we want, but this shows that the program can find the fixed point even when you're off the initial guess.
 * **lorenz.f90** - Needs to be a fortran module called *user_functions*. Must contain at least the 3 subroutines here - *setup*, *fcn*, and *get_jac*. *setup* is run once at the beginning of the program and can be used if you need to allocate arrays or do something more advanced (see the kse examples). *fcn* returns the vector field with the unstable manifold. *get_jac* is used for finding the jacobian. The *get_jac* routine included here can actually be used in any module you create: it's an approximate method that works for any field.
 * **par** - A file with 36 lines that will be passed to the programs and subroutines. The 36th line must be the dimension of the system. The rest are up to you. None besides par(36) are used here, but in the kse examples par(1) is used for a parameter value.
 * **optiman_input** - A file with all the settings for the manifold that will eventually be computed. These settings are described in a future section.

1. Create the fixed point program and run it. This will create the fixed_point and eigenvector files that are needed by the program.

 ```bash
~/optiman(master ✗) make clean
rm -f *.x
cd obj && rm -f *.mod *.o
~/optiman(master ✗) make fixed_point.x user_fcn=lorenz
gfortran -Jobj -O3 -march=native -ffree-form -ffree-line-length-none -c -Iobj lorenz.f90 -o obj/lorenz.o
gfortran -Jobj -O3 -march=native -ffree-form -ffree-line-length-none -c -Iobj -o obj/fixed_point.o src/fixed_point.f90 -L/usr/lib -L/usr/lib -lm -lblas -llapack
gfortran -Jobj -O3 -march=native -ffree-form -ffree-line-length-none -o fixed_point.x obj/lorenz.o obj/fixed_point.o -L/usr/lib -L/usr/lib -lm -lblas -llapack
~/optiman(master ✗) ./fixed_point.x 
 In            4  newton iterations we found a fixed point at
 
  -3.1005556992693429E-016  -3.1005557654437919E-016  -4.0577462739308000E-016
 
 the norm of the vector field at this point is    8.4411423845788430E-015
 it's has been saved in the file fixed_point for use by optiman.
 
 Here we list out all the eigenvalues of the jacobian at that fixed point
 that have positive real parts. If you do not see at least 2 positive eigenvalues,
 or see the ones you want, you should manually reverse the vector field simply
 by making the output of fcn negative
 If you select complex eigenvalues, you must select them in their conjugate pairs

 i            1
 real part of eigenvalue    22.827723451163457     
 imaginary part of eigenvalue    0.0000000000000000     

 i            3
 real part of eigenvalue    2.6666666666666665     
 imaginary part of eigenvalue    0.0000000000000000     
 
 Select the value of the first eigenvalue/vector you would like
 to use
1
 Select the value of the first eigenvalue/vector you would like
 to use
3
 your selections have been saved to eigval1,eigval2,eigvec1,eigvec2
 for use by optiman
~/optiman(master ✗) 
```

1. Create the manifold program and run it.

 ```bash
~/optiman(master ✗) make optiman.x user_fcn=lorenz
gfortran -Jobj -O3 -march=native -ffree-form -ffree-line-length-none -c -Iobj -o obj/ring.o src/ring.f90
gfortran -Jobj -O3 -march=native -ffree-form -ffree-line-length-none -c -Iobj -o obj/driver.o src/driver.f90
gfortran -Jobj -O3 -march=native -ffree-form -ffree-line-length-none -o optiman.x obj/lorenz.o obj/ring.o obj/driver.o -lm
~/optiman(master ✗) ./optiman.x 
   1.0000000000000000      %
   1.5000000000000000      %
...
   99.500000000000000      %
   100.00000000000000      %
~/optiman(master ✗) 
```

1. Store the results in the results folder for viewing later. This can be done with a simple command

 ```bash
~/optiman(master ✗) make store_results folder=my_first_manifold  
Results stored in results/my_first_manifold_150702095031
~/optiman(master ✗) ls
demos  LICENSE  makefile  obj  README.md  results  src
~/optiman(master ✗) 
```


# Visualizing your first manifold

1. Launch mayavi2. You may be able to launch this like any other program. You can launch it from the command line by typing *mayavi2*. In Ubuntu, the best way to launch it is by hitting Alt+F2 and then typing mayavi2 and hitting enter.

1. Click file->Run Python Script. Navigate to the project and then click the *src* folder. Double click on *Manifold.py*.

1. Use the built in python prompt, use the following commands to visualize your Manifold. You'll want to use your actual folder name rather than the one I have here.

 ```
my_first = Manifold('my_first_manifold_150702095031')
my_first.draw_manifold([0,1,2])
```

1. That's it, you should be able to look around at it! 
![alt text](https://raw.githubusercontent.com/dlithio/optiman/master/demos/images/lorenz.png "Mayavi Example
")
