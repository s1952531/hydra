 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
gamma-tilde: $(objects) $(fft_lib) $(sourcedir)/post/gamma-tilde.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/gamma-tilde.f90 -o gamma-tilde $(flags)

rossby-froude: $(objects) $(fft_lib) $(sourcedir)/post/rossby-froude.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/rossby-froude.f90 -o rossby-froude $(flags)

mean-flow: $(objects) $(fft_lib) $(sourcedir)/post/mean-flow.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/mean-flow.f90 -o mean-flow $(flags)

energy: $(objects) $(fft_lib) $(sourcedir)/post/energy.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/energy.f90 -o energy $(flags)

dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)

fopvbal: $(objects) $(fft_lib) $(sourcedir)/post/fopvbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/fopvbal.f90 -o fopvbal $(flags)

ispectra: $(objects) $(fft_lib) $(sourcedir)/post/ispectra.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/ispectra.f90 -o ispectra $(flags)

hspectrum: $(objects) $(fft_lib) $(sourcedir)/post/hspectrum.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/hspectrum.f90 -o hspectrum $(flags)

norms: $(objects) $(sourcedir)/post/norms.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/norms.f90 -o norms $(flags)

imbal: $(objects) $(sourcedir)/post/imbal.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/imbal.f90 -o imbal $(flags)

genfg: $(objects) $(fft_lib) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

