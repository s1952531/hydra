 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

 # Set main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/variables.f90 $(sourcedir)/spectral.f90 \
          $(sourcedir)/contours.f90 $(sourcedir)/common.f90 $(sourcedir)/congen.f90 $(sourcedir)/evolution.f90 #Main f90 sources
objects = parameters.o constants.o variables.o spectral.o \
          contours.o common.o congen.o evolution.o #Main .o files 

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/deriv1d.f90 

spflags = -L/home/dgd/spherical_harmonics/SHTOOLS/lib -lSHTOOLS2.7 -I/home/dgd/spherical_harmonics/SHTOOLS/modules -L/home/dgd/fftw-3.3.2 -lfftw3 -lm -O3 -fno-second-underscore -s -m64 -mcmodel=medium

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: caps setup proxy_init_all proxy_post_all install

caps: $(objects) $(fft_lib) $(sourcedir)/caps.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/caps.f90 -o caps $(flags)

powerspec: $(objects) $(sourcedir)/powerspec.f90
	$(f90) parameters.o constants.o $(sourcedir)/powerspec.f90 -o powerspec $(spflags)

project: $(objects) $(sourcedir)/project.f90
	$(f90) parameters.o constants.o $(sourcedir)/project.f90 -o project $(flags)
genfg: $(objects) $(fft_lib) $(sourcedir)/genfg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/genfg.f90 -o genfg $(flags)
scatter: $(objects) $(fft_lib) $(sourcedir)/scatter.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/scatter.f90 -o scatter $(flags)
unsteady: $(objects) $(fft_lib) $(sourcedir)/unsteady.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/unsteady.f90 -o unsteady $(flags)
shape: $(objects) $(fft_lib) $(sourcedir)/shape.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/shape.f90 -o shape $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

 # Phony rules for setup etc.: 
.PHONY: all clean setup install
clean:
	rm -f *.o *.mod
setup:
	test -d $(installdir)/cont || mkdir $(installdir)/cont
	test -d $(installdir)/grid || mkdir $(installdir)/grid
	test -d $(installdir)/fine || mkdir $(installdir)/fine
install:
	find -type f -perm /a+x -exec /bin/mv {} $(installdir)/ \;

 # Include intialisation and post-processing rules 
 # if these exist:
-include init/init.make
-include post/post.make

 # Proxy rules to ensure whatever exists in the sub-directories
 # also gets built:
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
