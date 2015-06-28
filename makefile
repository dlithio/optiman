# User should set this directories. The vary by machine.
# If you've made sure to install these libraries, one way to try to find 
# the necessary files is to run the find commands listed before each
# variable. You could also contact the system administrator.
#
# find / -type f -name "*lapack*.a" 2>/dev/null
lapack_lib_dir=/usr/lib
# find / -type f -name "*blas*.a" 2>/dev/null
blas_lib_dir=/usr/lib

# As long as you've installed ifort, nothing here should need to be
# changed.
FC = ifort
ifeq (ifort,$(FC))
fast_flags=-O3 -xhost -ipo -fp-model strict -i4
debug_flags=-i4 -O0 -traceback -g -check all -check bounds -debug all -fp-stack-check -fpe0 -ftrapuv -warn all
fixed_point_libs=-lblas -llapack
optiman_libs=
endif
ifeq (gfortran,$(FC))
fast_flags=
debug_flags=
fixed_point_libs=-lm -lblas -llapack
optiman_libs=-lm
endif

FFLAGS=$(debug_flags)
# Select the flags and libraries that were specified
ifeq (fast,$(flags))
FFLAGS=$(fast_flags)
endif

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $@
	
fixed_point.o: fixed_point.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $@ -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

fixed_point.x: clean $(user_fcn).o fixed_point.o
	$(FC) $(FFLAGS) -o fixed_point.x *.o  -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

optiman.x: clean $(user_fcn).o ring.o driver.o
	$(FC) $(FFLAGS) -o optiman.x *.o $(optiman_libs)

nse: clean utility_mod.o projector_mod_phys.o nse_mod_physf.o auto_mod_physf.o status_mod.o nse.o ring.o driver.o mrgrnk.o
	$(FC) $(FFLAGS) -o nse_optiman.x *.o $(optiman_libs)
	
nse_fixed_point.x: clean utility_mod.o projector_mod_phys.o nse_mod_physf.o auto_mod_physf.o status_mod.o nse.o fixed_point.o mrgrnk.o
	$(FC) $(FFLAGS) -o nse_fixed_point.x *.o -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

clean:
	rm -f *.exe *.mod *.o *.x
