
#
# Comment by DM (2007-07-04)
# The problem on ormen
# - there is a f77 (or g77) but it is based on gcc 3.3.5
# - when you try to use it with the cpp preprocessor (option
#   -x f77-cpp-input), it says :
#   f77: installation problem, cannot exec `cc1': No such file or directory
# - There are indeed some files in /usr/lib64/gcc-lib/x86_64-suse-linux/3.3.5
#   but not enough (there is no cc1 !)
# - In fact, the installed version of gcc is 4.0.2 (see gcc --version)
#   and indeed cc1 can be found in /usr/lib64/gcc/x86_64-suse-linux/4.0.2
# - This version of gcc is said to be installed with f95 langage support
#   (see gcc --v) but f95 or g95 is unknown
#
# So we have either :
# - a fotran compiler, but with no support of CPP directives, and no way to mix
#   it with C/C++ code
# - or a C/C++ working environment with no fortran support.
#
# We get fortran support with -x switch
# and we add the needed libraries
#

# Valid options
#    SAME_RANDOM_VALUES
#        A pseudo-random function is used in the simulation. If we want to get the
#        same results with the standard version and the 2 steps versions, some
#        care must be taken, this symbol activate the corresponding code

#LPJ_OPTIONS_STD
#LPJ_OPTIONS_STEP_1A
#LPJ_OPTIONS_STEP_1B
#LPJ_OPTIONS_STEP_2

#
# default options (main is in fortran)
#

# f77-cpp-input means pass the file through cpp

# NETCDF     = /opt/local/netcdf
# NETCDF     = /usr/local
# FC         = g77
# FC         =ifort

NETCDF     = /home/doug/Documents2/ncdf/netcdf-3.6.3/cxx
FC         = gfortran
FCOPTIONS  = -x f77-cpp-input
# FCOPTIONS  = -x f77-cpp-input -O3
# FCOPTIONS  = -fpp -O3
CPPOPTIONS = -c -L$(NETCDF)/include/
CXX        = c++
CXXOPTIONS = -O3
LD         = $(FC)
LDLIBS     = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -L/usr/lib/gcc/x86_64-linux-gnu/4.9/ -lstdc++
#LDLIBS     = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -L/usr/lib/gcc/darwin/3.3 -lstdc++
#LDLIBS     = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf #-L/lib/gcc/i686-pc-cygwin/4.5.0 -lstdc++ /usr/bin/x86_64-w64-mingw32-gcc
#
# Options for known enviroment
# Use of Intel 8.1 compiler seems to be frequent
#

ifeq ($(TARGET),INTEL_FC_81)
  FC        = ifort
  FCOPTIONS = -fpp -O3
  NETCDF    = /opt/local/intel_fc_81
  LDLIBS    = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -lstdc++
endif

ifeq ($(TARGET),ORMEN_GCC4)
  FC     = gcc
  NETCDF = /home/ggxyz/software/netcdf
  LD     = c++
  LDLIBS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -lfrtbegin -lg2c -lgfortran
endif

EXE = motif-lpj motif-lpj-step1a motif-lpj-step1b motif-lpj-step2

EXE1 = motif-lpj

EXE2 = motif-lpj-step1a motif-lpj-step2

EXE3 = motif-lpj-step1a motif-lpj-step1b motif-lpj-step2

#EXE = motif-lpj

OBJS = lpjmain-std.o lpjmain-step1a.o lpjmain-step1b.o lpjmain-step2.o \
       lpjio-std.o   lpjio-step1a.o   lpjio-step1b.0   lpjio-step2.o   \

#OBJS = lpjmain-std.o lpjio-std.o   Array.o

.PHONY: clean

all: $(EXE)

one_step : $(EXE1)

two_step : $(EXE2)

three_step : $(EXE3)

clean:
	rm -f $(OBJS) $(EXE)

motif-lpj: lpjmain-std.o lpjio-std.o Array.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

motif-lpj-step1a: lpjmain-step1a.o Array.o lpjio-step1a.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

motif-lpj-step1b: lpjmain-step1b.o lpjio-step1b.o Array.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

motif-lpj-step2: lpjmain-step2.o lpjio-step2.o Array.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

lpjio-std.o: lpjio.cpp
	$(CXX) $(LPJ_OPTIONS_STD) $(CPPOPTIONS) $(CXXOPTIONS) -o $@ $^

lpjio-step1a.o: lpjio.cpp
	$(CXX) $(LPJ_OPTIONS_STEP_1A) $(CPPOPTIONS) $(CXXOPTIONS) -DLPJ_STEP_1A -o $@ -c $^

lpjio-step1b.o: lpjio.cpp
	$(CXX) $(LPJ_OPTIONS_STEP_1B) $(CPPOPTIONS) $(CXXOPTIONS) -DLPJ_STEP_1B -o $@ -c $^

lpjio-step2.o: lpjio.cpp
	$(CXX) $(LPJ_OPTIONS_STEP_2) $(CPPOPTIONS) $(CXXOPTIONS) -DLPJ_STEP_2 -o $@ -c $^

lpjmain-std.o: lpjmain.f
	$(FC) $(LPJ_OPTIONS_STD) $(FCOPTIONS) -o $@ -c $^

lpjmain-step1a.o: lpjmain.f
	$(FC) $(LPJ_OPTIONS_STEP_1A) $(FCOPTIONS) -DLPJ_STEP_1A -o $@ -c $^

lpjmain-step1b.o: lpjmain.f
	$(FC) $(LPJ_OPTIONS_STEP_1B) $(FCOPTIONS) -DLPJ_STEP_1B -o $@ -c $^

lpjmain-step2.o: lpjmain.f
	$(FC) $(LPJ_OPTIONS_STEP_2) $(FCOPTIONS) -DLPJ_STEP_2 -o $@ -c $^

Array.o: array.cpp
	$(CXX) $(LPJ_OPTIONS_STD) $(CPPOPTIONS) $(CXXOPTIONS) -o $@ -c $^
