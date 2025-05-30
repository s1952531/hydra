 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

 # Set main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/setup.f90
objects = parameters.o constants.o setup.o

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/sta2dfft.f90

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium -lopenblas #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace -lopenblas  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: epbeq extract install

epbeq: $(objects) $(fft_lib) $(sourcedir)/epbeq.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/epbeq.f90 -o epbeq $(flags)

extract: $(objects) $(fft_lib) $(sourcedir)/extract.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/extract.f90 -o extract $(flags)

setuptest: $(objects) $(fft_lib) $(sourcedir)/setuptest.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/setuptest.f90 -o setuptest $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

 # Phony rules for setup etc.: 
.PHONY: all clean setup install
clean:
	rm -f *.o *.mod
install:
	find -type f -perm /a+x -exec /bin/mv {} $(installdir)/ \;
