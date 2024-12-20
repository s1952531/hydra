 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all':
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))
# Check if LAPACK or OpenBLAS is installed
lapack_installed := $(shell command -v lapack-config >/dev/null 2>&1 && echo true || echo false)
openblas_installed := $(shell ldconfig -p | grep -q libopenblas && echo true || echo false)

# Debugging output (optional, to check the values during the make process)
$(info LAPACK installed: $(lapack_installed))
$(info OpenBLAS installed: $(openblas_installed))

#---------------------------------------------------------------------------------
# Rules:
vertical: $(objects) $(sourcedir)/init/vertical.f90
ifeq ($(lapack_installed), true)
	$(f90) parameters.o $(sourcedir)/init/vertical.f90 -o vertical $(flags) -llapack
else ifeq ($(openblas_installed), true)
	$(f90) parameters.o $(sourcedir)/init/vertical.f90 -o vertical $(flags) -lopenblas
else
	$(error Either LAPACK or OpenBLAS must be installed.)
endif

eddy: $(objects) $(sourcedir)/init/eddy.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/eddy.f90 -o eddy $(flags)

rest: $(objects) $(sourcedir)/init/rest.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest.f90 -o rest $(flags)

spin_down: $(objects) $(sourcedir)/init/spin_down.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/spin_down.f90 -o spin_down $(flags)

ellipse_bath: $(objects) $(sourcedir)/init/ellipse_bath.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipse_bath.f90 -o ellipse_bath $(flags)

sine_bath: $(objects) $(sourcedir)/init/sine_bath.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/sine_bath.f90 -o sine_bath $(flags)

random_bath: $(fft_lib) $(objects) $(sourcedir)/init/random_bath.f90
	$(f90) $(fft_lib) parameters.o constants.o $(sourcedir)/init/random_bath.f90 -o random_bath $(flags)

random_km2_bath: $(fft_lib) $(objects) $(sourcedir)/init/random_km2_bath.f90
	$(f90) $(fft_lib) parameters.o constants.o $(sourcedir)/init/random_km2_bath.f90 -o random_km2_bath $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)
