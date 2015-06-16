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
FC = ifort
fast_flags=-O3 -xhost -ipo -fp-model strict -i4
debug_flags=-i4 -O0 -traceback -g -check all -check bounds -debug all -fp-stack-check -fpe0 -ftrapuv -warn all
libs=-lblas -llapack -lfftw3

FFTW3.o: FFTW3.f90
	$(FC) $(FFLAGS) -I$(fftw_include_dir) -c $*.f90 -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $@

ring.o: FFTW3.o

optiman.x: clean $(user_fcn).o ring.o driver.o
	$(FC) $(FFLAGS) -o optiman.x *.o -L$(fftw_lib_dir) -L$(lapack_lib_dir) -L$(blas_lib_dir) $(libs)
	
nse: utility_mod.o projector_mod_phys.o nse_mod_physf.o auto_mod_physf.o status_mod.o example1.o ring.o driver.o mrgrnk.o
	$(FC) $(FFLAGS) -o nse_optiman.exe *.o $(libs)

clean:
	rm -f *.exe *.mod *.o *.x
