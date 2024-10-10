 # Set existence of directory variable used in main makefile:
post_exists = true

 # Determine f90 codes existing in post directory for making with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#-----------------------------------------------------------------------------
 #Rules:
r4toc2: $(objects) post/r4toc2.f90
	$(f90) parameters.o constants.o post/r4toc2.f90 -o r4toc2 $(flags)

p2g: $(objects) post/p2g.f90
	$(f90) parameters.o constants.o post/p2g.f90 -o p2g $(flags)

image: $(objects) post/image.f90
	$(f90) parameters.o constants.o post/image.f90 -o image $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

