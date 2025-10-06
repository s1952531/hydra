 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
select: $(objects) $(sourcedir)/post/select.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/select.f90 -o select $(flags)

extend: $(objects) $(fft_lib) $(sourcedir)/post/extend.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/extend.f90 -o extend $(flags)

froude: $(objects) $(fft_lib) $(sourcedir)/post/froude.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/froude.f90 -o froude $(flags)

cross: $(objects) $(fft_lib) $(sourcedir)/post/cross.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/cross.f90 -o cross $(flags)

new_cross: $(objects) $(fft_lib) $(sourcedir)/post/new_cross.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/new_cross.f90 -o new_cross $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

