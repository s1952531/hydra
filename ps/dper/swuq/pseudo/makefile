 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .

 # Set main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/spectral.f90 \
          $(sourcedir)/common.f90 #Main f90 sources
objects = parameters.o constants.o spectral.o common.o #Main .o files 

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/sta2dfft.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(test_flags)

#----------------------------------------------------------------------
 # Rules:
all: swgw packet

swgw: $(objects) $(fft_lib) $(sourcedir)/swgw.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/swgw.f90 -o swgw $(flags)

packet: $(objects) $(sourcedir)/packet.f90
	$(f90) parameters.o constants.o $(sourcedir)/packet.f90 -o packet $(flags)

gausstopo: $(objects) $(sourcedir)/gausstopo.f90
	$(f90) parameters.o constants.o $(sourcedir)/gausstopo.f90 -o gausstopo $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

 # Phony rules for setup etc.: 
.PHONY: all clean
clean:
	rm -f *.o *.mod