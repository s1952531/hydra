 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
balance: $(objects) $(fft_lib) $(sourcedir)/post/balance.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/balance.f90 -o balance $(flags)

psi: $(objects) $(fft_lib) $(sourcedir)/post/psi.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/psi.f90 -o psi $(flags)

g2c: $(objects) $(sourcedir)/post/g2c.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/g2c.f90 -o g2c $(flags)

slice: $(objects) $(sourcedir)/post/slice.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/slice.f90 -o slice $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

