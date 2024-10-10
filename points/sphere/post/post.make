 # Set existence of directory variable used in main makefile:
post_exists = true

 # Determine f90 codes existing in post directory for making with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#-----------------------------------------------------------------------------
 #Rules:
log-d2-pdf: $(objects) post/log-d2-pdf.f90
	$(f90) parameters.o constants.o post/log-d2-pdf.f90 -o log-d2-pdf $(flags)
epdf-dipoles: $(objects) post/epdf-dipoles.f90
	$(f90) parameters.o constants.o post/epdf-dipoles.f90 -o epdf-dipoles $(flags)
image: $(objects) post/image.f90
	$(f90) parameters.o constants.o post/image.f90 -o image $(flags)
p2g: $(objects) post/p2g.f90
	$(f90) parameters.o constants.o post/p2g.f90 -o p2g $(flags)
espec: post/espec.f90
	$(f90) post/espec.f90 -o espec $(flags) $(shflags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

