 # Set main source and object files:
sources = parameters.f90 constants.f90
objects = parameters.o constants.o

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/deriv1d.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
#test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
diffuse: $(objects) $(fft_lib) diffuse.f90
	$(f90) $(fft_lib) $(objects) diffuse.f90 -o diffuse $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

 # Phony rules for setup etc.: 
.PHONY: clean
clean:
	rm -f *.o *.mod
