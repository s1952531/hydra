 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
vertical: $(objects) $(sourcedir)/init/vertical.f90
	$(f90) parameters.o $(sourcedir)/init/vertical.f90 -o vertical $(flags) -lopenblas

rod: $(objects) $(sourcedir)/init/rod.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rod.f90 -o rod $(flags)

ellipse_bath: $(objects) $(sourcedir)/init/ellipse_bath.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipse_bath.f90 -o ellipse_bath $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


