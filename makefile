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

_nseobbjs = utility_mod.o projector_mod_phys.o nse_mod_physf.o auto_mod_physf.o status_mod.o nse.o mrgrnk.o
NSEOBJS = $(patsubst %,$(ODIR)/%,$(_nseobbjs))

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
	
$(ODIR)/utility_mod.o: utility_mod.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/projector_mod_phys.o: projector_mod_phys.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/nse_mod_physf.o: nse_mod_physf.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/auto_mod_physf.o: auto_mod_physf.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/status_mod.o: status_mod.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/mrgrnk.o: mrgrnk.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(FFLAGS) -c $(INC) -o $@ $<
	
$(ODIR)/fixed_point.o: $(SDIR)/fixed_point.f90
	$(FC) $(FFLAGS) -c $(INC) -o $@ $< -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

fixed_point.x: $(ODIR)/$(user_fcn).o $(FPOBJS)
	$(FC) $(FFLAGS) -o fixed_point.x $^ -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

optiman.x: $(ODIR)/$(user_fcn).o $(OPTOBJS)
	$(FC) $(FFLAGS) -o optiman.x $^ $(optiman_libs)

nse_fixed_point.x: $(NSEOBJS) $(ODIR)/$(user_fcn).o $(FPOBJS)
	$(FC) $(FFLAGS) -o nse_fixed_point.x $^ -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

nse_optiman.x: $(NSEOBJS) $(ODIR)/$(user_fcn).o $(OPTOBJS)
	$(FC) $(FFLAGS) -o nse_optiman.x $^ $(optiman_libs)
	
timestamp=$(shell date +"%y%m%d%H%M%S")
store_results:
ifdef folder
	mkdir results/$(folder)_$(timestamp)
	mv eig* results/$(folder)_$(timestamp)/
	mv *.x results/$(folder)_$(timestamp)/
	mv header results/$(folder)_$(timestamp)/
	mv fixed_point results/$(folder)_$(timestamp)/
	mv initial_guess results/$(folder)_$(timestamp)/
	mv par results/$(folder)_$(timestamp)/
	mv *.f90 results/$(folder)_$(timestamp)/
	mv output results/$(folder)_$(timestamp)/
	mv *_input results/$(folder)_$(timestamp)/
	mv fdot results/$(folder)_$(timestamp)/
	mv t_angle results/$(folder)_$(timestamp)/
	mv old_par results/$(folder)_$(timestamp)/
	mv auto_ndim_key results/$(folder)_$(timestamp)/
	mv M_key results/$(folder)_$(timestamp)/
	mv N_key results/$(folder)_$(timestamp)/
	mv kx_projections results/$(folder)_$(timestamp)/
	mv ky_projections results/$(folder)_$(timestamp)/
	rm -f version_info
	touch version_info
	echo "These results come from software version" >> version_info
	git rev-parse HEAD >> version_info
	mv version_info results/$(folder)_$(timestamp)/
else
	@echo "you need to specify a folder to put the files in using doing"
	@echo "make store_results folder=folder_name_to_stash_stuff_in"
endif
	
clean:
	rm -f *.x
	cd obj && rm -f *.mod *.o

