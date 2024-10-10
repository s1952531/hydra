 # Set existence of directory variable used in main makefile:
init_exists = true

 # Determine f90 codes existing in init directory for making with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#-----------------------------------------------------------------------------
 #Rules:
ring: $(objects) $(sourcedir)/init/ring.f90
	$(f90) $(objects) $(sourcedir)/init/ring.f90 -o ring $(flags)

monopoles: $(objects) $(sourcedir)/init/monopoles.f90
	$(f90) $(objects) $(sourcedir)/init/monopoles.f90 -o monopoles $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)
