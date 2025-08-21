 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
vstrip: $(objects) $(sourcedir)/init/vstrip.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/vstrip.f90 -o vstrip $(flags)

eddy: $(objects) $(sourcedir)/init/eddy.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/eddy.f90 -o eddy $(flags)

ranpv: $(objects) $(fft_lib) $(sourcedir)/init/ranpv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/ranpv.f90 -o ranpv $(flags)

vortex: $(objects) $(fft_lib) $(sourcedir)/init/vortex.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/vortex.f90 -o vortex $(flags)

init_tracer: $(objects) $(sourcedir)/init/init_tracer.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/init_tracer.f90 -o init_tracer $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)
