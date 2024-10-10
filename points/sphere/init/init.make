 # Set existence of directory variable used in main makefile:
init_exists = true

 # Determine f90 codes existing in init directory for making with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#-----------------------------------------------------------------------------
 #Rules:
dipoles: $(objects) $(sourcedir)/init/dipoles.f90
	$(f90) $(objects) $(sourcedir)/init/dipoles.f90 -o dipoles $(flags)

monopoles: $(objects) $(sourcedir)/init/monopoles.f90
	$(f90) $(objects) $(sourcedir)/init/monopoles.f90 -o monopoles $(flags)

pair: $(objects) $(sourcedir)/init/pair.f90
	$(f90) $(objects) $(sourcedir)/init/pair.f90 -o pair $(flags)

powdip: $(objects) $(sourcedir)/init/powdip.f90
	$(f90) $(objects) $(sourcedir)/init/powdip.f90 -o powdip $(flags)

powhex: $(objects) $(sourcedir)/init/powhex.f90
	$(f90) $(objects) $(sourcedir)/init/powhex.f90 -o powhex $(flags)

trihex: $(objects) $(sourcedir)/init/trihex.f90
	$(f90) $(objects) $(sourcedir)/init/trihex.f90 -o trihex $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)
