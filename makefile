# User should set this directories. The vary by machine.
# If you've made sure to install these libraries, one way to try to find 
# the necessary files is to run the find commands listed before each
# variable. You could also contact the system administrator.
#
# find / -type f -name "*fft*.a" 2>/dev/null
fftw_lib_dir=/usr/lib/x86_64-linux-gnu
# find / -type f -name "*lapack*.a" 2>/dev/null
lapack_lib_dir=/usr/lib
# find / -type f -name "*blas*.a" 2>/dev/null
blas_lib_dir=/usr/lib
# find / -type f -name "fft*.f03" 2>/dev/null
fftw_include_dir=/usr/include

FFLAGS=$(debug_flags)
# Select the flags and libraries that were specified
ifeq (fast,$(flags))
FFLAGS=$(fast_flags)
endif

# As long as you've installed gfortran, nothing here should need to be
# changed.
FC = gfortran
fast_flags=-O3 -march=native -std=f2008
debug_flags=-Og -g -Wall -Wextra -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -pedantic -fimplicit-none -fbounds-check -fbacktrace -fcheck=all -std=f2008
libs=-lblas -llapack -lfftw3

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $@
	
optiman.x: clean $(user_fcn).o ring.o driver.o
	$(FC) $(FFLAGS) -I$(fftw_include_dir) -o optiman.x *.o -L$(fftw_lib_dir) -L$(lapack_lib_dir) -L$(blas_lib_dir) $(libs)
	
clean:
	rm -f *.exe *.mod *.o *.x
