 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
slug: $(objects) $(sourcedir)/init/slug.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/slug.f90 -o slug $(flags)

straka: $(objects) $(sourcedir)/init/straka.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/straka.f90 -o straka $(flags)

robert: $(objects) $(sourcedir)/init/robert.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/robert.f90 -o robert $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


