 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
energy: $(objects) $(fft_lib) $(sourcedir)/post/energy.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/energy.f90 -o energy $(flags)

slice: $(objects) $(sourcedir)/post/slice.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/slice.f90 -o slice $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

