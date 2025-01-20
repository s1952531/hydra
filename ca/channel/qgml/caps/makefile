 # Set up source and install directories; 
 # these can be edited if the source/install structure needed is different:
sourcedir = .
installdir = ..

sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 \
          $(sourcedir)/spectral.f90 $(sourcedir)/contours.f90 \
          $(sourcedir)/congen.f90 $(sourcedir)/common.f90 \
          $(sourcedir)/evolution.f90 #Main f90 sources
objects = parameters.o constants.o spectral.o contours.o \
          congen.o common.o evolution.o #Main object files 

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/deriv1d.f90 \
          ~/hydra/lib/stafft/sta2dfft.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=large #Compiler flags
test_flags = -O3 -mcmodel=large -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: caps genfg rest setup install

caps: $(fft_lib) $(objects) $(sourcedir)/caps.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/caps.f90 -o caps $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

 # Suppress output apart from errors and warnings:
.SILENT:

.PHONY: clean setup
clean:
	rm -f *.o *.mod
setup:
	test -d $(installdir)/cont || mkdir $(installdir)/cont
	test -d $(installdir)/fine || mkdir $(installdir)/fine
	test -d $(installdir)/evolution || mkdir $(installdir)/evolution
install:
	find -type f -perm /a+x -exec /bin/mv {} $(installdir)/ \;

 # Include intialisation and post-processing rules if these exist:
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
