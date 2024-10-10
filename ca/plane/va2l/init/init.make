 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
ranpv: $(objects) $(fft_lib) $(sourcedir)/init/ranpv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/ranpv.f90 -o ranpv $(flags)

wave: parameters.o constants.o $(sourcedir)/init/wave.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/wave.f90 -o wave $(flags)

vstrip: parameters.o constants.o $(sourcedir)/init/vstrip.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/vstrip.f90 -o vstrip $(flags)

balinit: $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/balinit.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/balinit.f90 -o balinit $(flags)

swbalinit: $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/swbalinit.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/swbalinit.f90 -o swbalinit $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


