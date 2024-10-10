 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/spectral.f90
objects = parameters.o constants.o spectral.o

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/deriv1d.f90 ~/hydra/lib/stafft/sta2dfft.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(test_flags)

#----------------------------------------------------------------------
 # Rules:
mapgen: $(fft_lib) $(objects) $(sourcedir)/mapgen.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/mapgen.f90 -o mapgen $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

.PHONY: clean setup
clean:
	rm -f *.o *.mod
