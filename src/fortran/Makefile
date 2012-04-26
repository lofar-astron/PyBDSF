# All the routines which are needed for libraries for PyBDSM are compiled here
# using gnu95, f2py and suitable flags depending on os. Include files are in ../includes
# $(FC) -c *.f -o *.o

DIR=../../fBDSM
CFLAGS = -I$(DIR)/includes
ODIR=./
LIBDIR=../../fBDSM/libs
IODIR=../../fBDSM/io
MDIR=../../fBDSM/wout_pgplot

ifeq ($(shell uname -m), x86_64)
 FC=gfortran -fPIC -fno-second-underscore $(CFLAGS)     
else   # ! i686
 FC=/usr/bin/g77 -g -fno-second-underscore -fbounds-check $(CFLAGS)
endif

ifeq ($(shell uname -s), Darwin)
 FC=gfortran -fbounds-check $(CFLAGS)
endif

OBJ1 = iland_mat2list.o writearray.o gaul2gaulbin.o writefitshead.o sub_prog.o \
 readarray.o get_keyword.o check.o 

OBJ1_f = iland_mat2list.f writearray.f gaul2gaulbin.f writefitshead.f

all: _pytesselate.so _py2fbdsm.so

LIBS=$(LIBDIR)/libfbdsm_woutpg.a

_pytesselate.pyf : 
	f2py2.5 pytess_simple.f pytess_roundness.f -m _pytesselate -h _pytesselate.pyf

_pytesselate.so: _pytesselate.pyf pytess_simple.o pytess_roundness.o
	f2py2.5 --fcompiler=gnu95 -c _pytesselate.pyf pytess_simple.o pytess_roundness.o dummy.f 

_py2fbdsm.pyf : 
	cp $(IODIR)/readarray.f $(IODIR)/writearray.f $(IODIR)/writefitshead.f ./
	cp $(MDIR)/check.f $(MDIR)/get_keyword.f $(MDIR)/iland_mat2list.f $(MDIR)/sub_prog.f ./
	f2py2.5 $(OBJ1_f) -m _py2fbdsm -h _py2fbdsm.pyf

_py2fbdsm.so: _py2fbdsm.pyf $(OBJ1)
	f2py2.5 --fcompiler=gnu95 -c _py2fbdsm.pyf $(OBJ1) $(LIBS) dummy.f


%.o: %.f
	$(FC) -c $< -o $@

clean: 
	rm -fr *.o 

cleanall: 
	rm -fr *.o *.so *.pyf


