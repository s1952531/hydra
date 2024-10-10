 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:

relax: $(objects) $(sourcedir)/init/relax.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/relax.f90 -o relax $(flags)

rest-state: $(objects) $(sourcedir)/init/rest-state.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest-state.f90 -o rest-state $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


