# User should set this directories. The vary by machine.
# If you've made sure to install these libraries, one way to try to find 
# the necessary files is to run the find commands listed before each
# variable. You could also contact the system administrator.
#
# find / -type f -name "*lapack*.a" 2>/dev/null
lapack_lib_dir=/usr/lib
# find / -type f -name "*blas*.a" 2>/dev/null
blas_lib_dir=/usr/lib

ODIR = obj
SDIR = src
INC = -Iobj

_fpobjs = fixed_point.o
FPOBJS = $(patsubst %,$(ODIR)/%,$(_fpobjs))

_optobbjs = ring.o driver.o
OPTOBJS = $(patsubst %,$(ODIR)/%,$(_optobbjs))

# As long as you've installed ifort, nothing here should need to be
# changed.
FC = ifort
ifeq (ifort,$(FC))
fast_flags=-module $(ODIR) -O3 -xhost -ipo -fp-model strict -i4
debug_flags=-module $(ODIR) -i4 -O0 -traceback -g -check all -check bounds -debug all -fp-stack-check -fpe0 -ftrapuv -warn all
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

$(ODIR)/$(user_fcn).o: $(user_fcn).f90
	$(FC) $(FFLAGS) -c $(INC) $(user_fcn).f90 -o $@

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(FFLAGS) -c $(INC) -o $@ $<
	
$(ODIR)/fixed_point.o: $(SDIR)/fixed_point.f90
	$(FC) $(FFLAGS) -c $(INC) -o $@ $< -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

fixed_point.x: $(ODIR)/$(user_fcn).o $(FPOBJS)
	$(FC) $(FFLAGS) -o fixed_point.x $^ -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

optiman.x: $(ODIR)/$(user_fcn).o $(OPTOBJS)
	$(FC) $(FFLAGS) -o optiman.x $^ $(optiman_libs)

nse: clean utility_mod.o projector_mod_phys.o nse_mod_physf.o auto_mod_physf.o status_mod.o nse.o ring.o driver.o mrgrnk.o
	$(FC) $(FFLAGS) -o nse_optiman.x *.o $(optiman_libs)
	
nse_fixed_point.x: clean utility_mod.o projector_mod_phys.o nse_mod_physf.o auto_mod_physf.o status_mod.o nse.o fixed_point.o mrgrnk.o
	$(FC) $(FFLAGS) -o nse_fixed_point.x *.o -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

store_results:
ifdef folder
	mkdir results/$(folder)
	mv eig* results/$(folder)/
	mv *.x results/$(folder)/
	mv header results/$(folder)/
	mv fixed_point results/$(folder)/
	mv initial_guess results/$(folder)/
	mv par results/$(folder)/
	mv *.f90 results/$(folder)/
	mv output results/$(folder)/
	mv *_input results/$(folder)/
	mv fdot results/$(folder)/
	mv t_angle results/$(folder)/
	touch version_info
	echo "These results come from software version" >> version_info
	git rev-parse HEAD >> version_info
	mv version_info results/$(folder)/
else
	@echo "you need to specify a folder to put the files in using doing"
	@echo "make store_results folder=folder_name_to_stash_stuff_in"
endif
	
clean:
	rm -f *.exe *.mod *.o *.x
	cd obj && rm -f *.exe *.mod *.o *.x

