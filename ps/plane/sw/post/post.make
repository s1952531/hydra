 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
rossby-froude: $(objects) $(fft_lib) $(sourcedir)/post/rossby-froude.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/rossby-froude.f90 -o rossby-froude $(flags)

gamma-tilde: $(objects) $(fft_lib) $(sourcedir)/post/gamma-tilde.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/gamma-tilde.f90 -o gamma-tilde $(flags)

velocity: $(objects) $(fft_lib) $(sourcedir)/post/velocity.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/velocity.f90 -o velocity $(flags)

dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)

hspectrum: $(objects) $(fft_lib) $(sourcedir)/post/hspectrum.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/hspectrum.f90 -o hspectrum $(flags)

ispectra: $(objects) $(fft_lib) $(sourcedir)/post/ispectra.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/ispectra.f90 -o ispectra $(flags)

norms: $(objects) $(sourcedir)/post/norms.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/norms.f90 -o norms $(flags)

vertical-velocity: $(objects) $(sourcedir)/post/vertical-velocity.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/vertical-velocity.f90 -o vertical-velocity $(flags)

imbal: $(objects) $(sourcedir)/post/imbal.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/imbal.f90 -o imbal $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

