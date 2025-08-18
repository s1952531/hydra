 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
bci2l: $(objects) $(sourcedir)/post/bci2l.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/bci2l.f90 -o bci2l $(flags)

 # Suppress output apart from errors and warnings:
.SILENT:

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

