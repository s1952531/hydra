 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
ispectra: $(objects) $(fft_lib) $(sourcedir)/post/ispectra.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/ispectra.f90 -o ispectra $(flags)
profile: $(objects) $(fft_lib) $(sourcedir)/post/profile.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/profile.f90 -o profile $(flags)
zonal: $(objects) $(fft_lib) $(sourcedir)/post/zonal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/zonal.f90 -o zonal $(flags)
ortho: $(objects) post/ortho.f90
	$(f90) parameters.o constants.o post/ortho.f90 -o ortho $(flags)
genfg: $(objects) $(fft_lib) post/genfg.f90
	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

