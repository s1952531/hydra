 # Set main source and object files:
sources = parameters.f90 constants.f90 #Main f90 sources
objects = parameters.o constants.o #Main .o files 

# Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: pam dambreak profile running hspectrum

pam: $(objects) pam.f90
	$(f90) $(objects) pam.f90 -o pam $(flags)

relax: $(objects) relax.f90
	$(f90) $(objects) relax.f90 -o relax $(flags)

dambreak: $(objects) dambreak.f90
	$(f90) $(objects) dambreak.f90 -o dambreak $(flags)

steadyic: $(objects) steadyic.f90
	$(f90) $(objects) steadyic.f90 -o steadyic $(flags)

profile: $(objects) profile.f90
	$(f90) $(objects) profile.f90 -o profile $(flags)

running: $(objects) running.f90
	$(f90) $(objects) running.f90 -o running $(flags)

equil: $(objects) equil.f90
	$(f90) $(objects) equil.f90 -o equil $(flags)

test_equil: $(objects) test_equil.f90
	$(f90) $(objects) test_equil.f90 -o test_equil $(flags)

hspectrum: $(objects) $(fft_lib) hspectrum.f90
	$(f90) $(fft_lib) $(objects) hspectrum.f90 -o hspectrum $(flags)

$(objects): $(sources) 
	$(f90) $(sources) -c $(flags)

.PHONY: clean
clean:
	rm -f *.o *.mod
