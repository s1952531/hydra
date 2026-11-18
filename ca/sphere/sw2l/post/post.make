 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#-----------------------------------------------------------------------------
 #Rules:
dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)

verify_zdg: $(objects) $(fft_lib) $(sourcedir)/post/verify_zdg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/verify_zdg.f90 -o verify_zdg $(flags)

check_h: $(objects) $(fft_lib) $(sourcedir)/post/check_h.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/check_h.f90 -o check_h $(flags)

check_inversion: $(objects) $(fft_lib) $(sourcedir)/post/check_inversion.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/check_inversion.f90 -o check_inversion $(flags)

test_inversion: $(objects) $(fft_lib) $(sourcedir)/post/test_inversion.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/test_inversion.f90 -o test_inversion $(flags)

test_project: $(objects) post/test_project.f90
	$(f90) parameters.o constants.o post/test_project.f90 -o test_project $(flags)

zonal: $(objects) $(fft_lib) $(sourcedir)/post/zonal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/zonal.f90 -o zonal $(flags)

hov: $(objects) $(fft_lib) post/hov.f90
	$(f90) $(fft_lib) $(objects) post/hov.f90 -o hov $(flags)

ortho: $(objects) post/ortho.f90
	$(f90) parameters.o constants.o post/ortho.f90 -o ortho $(flags)

ini_ortho: $(objects) post/ini_ortho.f90
	$(f90) parameters.o constants.o post/ini_ortho.f90 -o ini_ortho $(flags)

avgdiag: $(objects) post/avgdiag.f90
	$(f90) parameters.o post/avgdiag.f90 -o avgdiag $(flags)

 #NOT WORKING YET:
 #genfg: $(objects) $(fft_lib) post/genfg.f90
 #	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

