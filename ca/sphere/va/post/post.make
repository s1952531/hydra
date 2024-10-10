 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
r4toc2: $(objects) post/r4toc2.f90
	$(f90) parameters.o constants.o post/r4toc2.f90 -o r4toc2 $(flags)
flowprop: $(objects) post/flowprop.f90
	$(f90) parameters.o constants.o post/flowprop.f90 -o flowprop $(flags)
kinetic: $(objects) post/kinetic.f90
	$(f90) parameters.o constants.o post/kinetic.f90 -o kinetic $(flags)
zonal: $(objects) post/zonal.f90
	$(f90) parameters.o constants.o post/zonal.f90 -o zonal $(flags)
genfg: $(objects) $(fft_lib) post/genfg.f90
	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

