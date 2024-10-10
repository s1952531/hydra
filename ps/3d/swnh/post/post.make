 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
rossby-froude: $(objects) $(fft_lib) $(sourcedir)/post/rossby-froude.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/rossby-froude.f90 -o rossby-froude $(flags)

energy: $(objects) $(fft_lib) $(sourcedir)/post/energy.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/energy.f90 -o energy $(flags)

accel: $(objects) $(fft_lib) $(sourcedir)/post/accel.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/accel.f90 -o accel $(flags)

mass: $(objects) $(fft_lib) $(sourcedir)/post/mass.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/mass.f90 -o mass $(flags)

velocity: $(objects) $(fft_lib) $(sourcedir)/post/velocity.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/velocity.f90 -o velocity $(flags)

variability: $(objects) $(fft_lib) $(sourcedir)/post/variability.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/variability.f90 -o variability $(flags)

gamma-tilde: $(objects) $(fft_lib) $(sourcedir)/post/gamma-tilde.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/gamma-tilde.f90 -o gamma-tilde $(flags)

dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)

hspectrum: $(objects) $(fft_lib) $(sourcedir)/post/hspectrum.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/hspectrum.f90 -o hspectrum $(flags)

decomp_spec: $(objects) $(fft_lib) $(sourcedir)/post/decomp_spec.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/decomp_spec.f90 -o decomp_spec $(flags)

ispectra: $(objects) $(fft_lib) $(sourcedir)/post/ispectra.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/ispectra.f90 -o ispectra $(flags)

norms: $(objects) $(sourcedir)/post/norms.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/norms.f90 -o norms $(flags)

g2c: $(objects) $(sourcedir)/post/g2c.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/g2c.f90 -o g2c $(flags)

vertical-velocity: $(objects) $(sourcedir)/post/vertical-velocity.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/vertical-velocity.f90 -o vertical-velocity $(flags)

pressure: $(objects) $(sourcedir)/post/pressure.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/pressure.f90 -o pressure $(flags)

nh-pressure: $(objects) $(sourcedir)/post/nh-pressure.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/nh-pressure.f90 -o nh-pressure $(flags)

imbal: $(objects) $(sourcedir)/post/imbal.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/imbal.f90 -o imbal $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

