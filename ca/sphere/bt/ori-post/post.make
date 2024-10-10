 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
r4toc2: $(objects) post/r4toc2.f90
	$(f90) parameters.o constants.o post/r4toc2.f90 -o r4toc2 $(flags)
zonal: $(objects) post/zonal.f90
	$(f90) parameters.o constants.o post/zonal.f90 -o zonal $(flags)
measure: $(objects) post/measure.f90
	$(f90) parameters.o constants.o post/measure.f90 -o measure $(flags)
project: $(objects) post/project.f90
	$(f90) parameters.o constants.o post/project.f90 -o project $(flags)
enstrophy: $(objects) post/enstrophy.f90
	$(f90) parameters.o constants.o post/enstrophy.f90 -o enstrophy $(flags)
genfg: $(objects) $(fft_lib) post/genfg.f90
	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)
scatter: $(objects) $(fft_lib) post/scatter.f90
	$(f90) $(fft_lib) $(objects) post/scatter.f90 -o scatter $(flags)
unsteady: $(objects) $(fft_lib) post/unsteady.f90
	$(f90) $(fft_lib) $(objects) post/unsteady.f90 -o unsteady $(flags)
shape: $(objects) $(fft_lib) post/shape.f90
	$(f90) $(fft_lib) $(objects) post/shape.f90 -o shape $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

