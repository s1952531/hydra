 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:

inigamma: $(objects) $(fft_lib) $(sourcedir)/init/inigamma.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/inigamma.f90 -o inigamma $(flags)

randompv: $(objects) $(fft_lib) $(sourcedir)/init/randompv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/randompv.f90 -o randompv $(flags)

rossby: $(sourcedir)/init/rossby.f90
	$(f90) $(sourcedir)/init/rossby.f90 -o rossby $(flags)

doublejet: $(objects) $(sourcedir)/init/doublejet.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/doublejet.f90 -o doublejet $(flags)

bickley: $(objects) $(sourcedir)/init/bickley.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/bickley.f90 -o bickley $(flags)

thermal-bickley: $(objects) $(sourcedir)/init/thermal-bickley.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/thermal-bickley.f90 -o thermal-bickley $(flags)

rest-state: $(objects) $(sourcedir)/init/rest-state.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest-state.f90 -o rest-state $(flags)

relax: $(objects) $(sourcedir)/init/relax.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/relax.f90 -o relax $(flags)

quadtopo: $(objects) $(sourcedir)/init/quadtopo.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/quadtopo.f90 -o quadtopo $(flags)

harmonic: $(objects) $(sourcedir)/init/harmonic.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/harmonic.f90 -o harmonic $(flags)

staircase: $(objects) $(sourcedir)/init/staircase.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/staircase.f90 -o staircase $(flags)

patches: $(objects) $(sourcedir)/init/patches.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/patches.f90 -o patches $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


