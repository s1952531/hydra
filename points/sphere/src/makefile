 # Set up source and install directories.  Note: these can be modified 
 # if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

 # Define main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/variables.f90 
objects = parameters.o constants.o variables.o

 # Set compiler options: 
norm_flags = -O3 -mcmodel=medium #Normal compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Testing flags

 # Define additional flags to use spherical harmonic and FFTW package:
shflags = -L/home/dgd/spherical_harmonics/SHTOOLS/lib -lSHTOOLS2.7 -I/home/dgd/spherical_harmonics/SHTOOLS/modules -L/home/dgd/fftw-3.3.2 -lfftw3 -lm -fno-second-underscore -s -m64

 # Choose fortran compiler and flags:
f90 = gfortran
flags = $(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: pvs proxy_init_all proxy_post_all install

pvs: $(objects) $(sourcedir)/pvs.f90
	$(f90) $(objects) $(sourcedir)/pvs.f90 -o pvs $(flags)

$(objects): $(sources) 
	$(f90) $(sources) -c $(flags)

 # Phony rules for clean up etc.: 
.PHONY: all clean install
clean:
	rm -f *.o *.mod
install:
	find -type f -perm /a+x -exec /bin/mv {} $(installdir)/ \;

 # Include intialisation and post-processing rules if these exist:
-include init/init.make
-include post/post.make

 # Proxy rules to ensure whatever exists in the subdirectories also gets built:
ifdef init_exists
proxy_init_all: init_all
else 
proxy_init_all:
	 @echo No init directory found 
endif

ifdef post_exists
proxy_post_all: post_all
else 
proxy_post_all:
	 @echo No post directory found 
endif
