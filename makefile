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
INC = -I$(ODIR)

_fpobjs = mrgrnk.o fixed_point.o
FPOBJS = $(patsubst %,$(ODIR)/%,$(_fpobjs))

_optobbjs = ring.o driver.o
OPTOBJS = $(patsubst %,$(ODIR)/%,$(_optobbjs))

_nseobbjs = utility_mod.o projector_mod_phys.o

FC = gfortran
ifeq (intel,$(compliler))
FC = ifort
endif

# As long as you've installed ifort, nothing here should need to be
# changed.

ifeq (ifort,$(FC))
fast_flags=-module $(ODIR) -O3 -xhost -ipo -fp-model strict -i4
debug_flags=-module $(ODIR) -i4 -O0 -traceback -g -check all -check bounds -debug all -fp-stack-check -fpe0 -ftrapuv -warn all
fixed_point_libs= -lblas -llapack
optiman_libs=-lblas
_nseobbjs+= nse_mod_physf_intel.o 
endif
ifeq (gfortran,$(FC))
fast_flags=-J$(ODIR) -O3 -march=native -ffree-form -ffree-line-length-none
debug_flags=-J$(ODIR) -ffpe-trap=invalid,zero,overflow -fimplicit-none -ffree-form -ffree-line-length-none -fbounds-check -O0 -g -Waliasing -Wall -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation -Wno-align-commons -Wno-tabs -Wreal-q-constant -Wsurprising -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wtarget-lifetime -fbacktrace
fixed_point_libs=-lm -lblas -llapack
optiman_libs=-lm -lblas
_nseobbjs+=nse_mod_physf_gnu.o 
endif
_nseobbjs += status_mod.o auto_mod_physf.o nse.o mrgrnk.o
NSEOBJS = $(patsubst %,$(ODIR)/%,$(_nseobbjs))

FFLAGS=$(fast_flags)
# Select the flags and libraries that were specified
ifeq (debug,$(flags))
FFLAGS=$(debug_flags)
endif

$(ODIR)/$(user_fcn).o: $(user_fcn).f90
	$(FC) $(FFLAGS) -c $(INC) $(user_fcn).f90 -o $@
	
$(ODIR)/utility_mod.o: utility_mod.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/projector_mod_phys.o: projector_mod_phys.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/auto_mod_physf.o: auto_mod_physf.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/nse_mod_physf_intel.o: nse_mod_physf_intel.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/nse_mod_physf_gnu.o: nse_mod_physf_gnu.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@
	
$(ODIR)/status_mod.o: status_mod.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $@

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(FFLAGS) -c $(INC) -o $@ $<

fixed_point.x: $(ODIR)/$(user_fcn).o $(FPOBJS)
	$(FC) $(FFLAGS) -o fixed_point.x $^ -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

optiman.x: $(ODIR)/$(user_fcn).o $(OPTOBJS)
	$(FC) $(FFLAGS) -o optiman.x $^ -L$(blas_lib_dir) $(optiman_libs)

nse_fixed_point.x: $(NSEOBJS) $(ODIR)/$(user_fcn).o $(FPOBJS)
	$(FC) $(FFLAGS) -o nse_fixed_point.x $^ -L$(lapack_lib_dir) -L$(blas_lib_dir) $(fixed_point_libs)

nse_optiman.x: $(NSEOBJS) $(ODIR)/$(user_fcn).o $(OPTOBJS)
	$(FC) $(FFLAGS) -o nse_optiman.x $^ -L$(blas_lib_dir) $(optiman_libs)
	
timestamp:=$(shell /bin/date +"%y%m%d%H%M%S")
store_results:
ifdef folder
	@mkdir results/$(folder)_$(timestamp)
	@mv eig* results/$(folder)_$(timestamp)/
	@mv q_matrix results/$(folder)_$(timestamp)/
	@mv *.x results/$(folder)_$(timestamp)/
	@mv header results/$(folder)_$(timestamp)/
	@mv fixed_point results/$(folder)_$(timestamp)/
	@mv initial_guess results/$(folder)_$(timestamp)/
	@mv par results/$(folder)_$(timestamp)/
	@mv *.f90 results/$(folder)_$(timestamp)/
	@mv output results/$(folder)_$(timestamp)/
	@mv *_input results/$(folder)_$(timestamp)/
	@mv fdot results/$(folder)_$(timestamp)/ 2>/dev/null
	@mv t_angle results/$(folder)_$(timestamp)/
	@mv interp_info results/$(folder)_$(timestamp)/
ifneq ("$(wildcard old_par)","")
	@mv old_par results/$(folder)_$(timestamp)/
	@mv auto_ndim_key results/$(folder)_$(timestamp)/
	@mv M_key results/$(folder)_$(timestamp)/
	@mv N_key results/$(folder)_$(timestamp)/
	@mv kx_projections results/$(folder)_$(timestamp)/
	@mv ky_projections results/$(folder)_$(timestamp)/
endif
	@rm -f version_info
	@touch version_info
	@echo "These results come from software version" >> version_info
	@git rev-parse HEAD >> version_info
	@mv version_info results/$(folder)_$(timestamp)/
	@echo "Results stored in results/$(folder)_$(timestamp)"
else
	@echo "you need to specify a folder to put the files in using doing"
	@echo "make store_results folder=folder_name_to_stash_stuff_in"
endif

clean_results:
	rm -f eig* 
	rm -f q_matrix
	rm -f *.x
	rm -f header
	rm -f fixed_point
	rm -f initial_guess 
	rm -f par 
	rm -f *.f90
	rm -f output 
	rm -f *_input 
	rm -f fdot 
	rm -f t_angle 
	rm -f old_par
	rm -f auto_ndim_key
	rm -f M_key 
	rm -f N_key 
	rm -f kx_projections 
	rm -f ky_projections 
	rm -f q_matrix 
	rm -f interp_info

clean:
	rm -f *.x
	cd obj && rm -f *.mod *.o
