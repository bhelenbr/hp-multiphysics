OBJ = allocate.o copy.o tobasis.o gtol.o output.o curvinit.o ptprobe.o mg_alloc.o \
rsdl.o minvrt.o movco.o movfin.o tstep.o nstage.o bdry.o tadvance.o length.o adapt.o \
surface.o drag.o block.o blocks.o rblocks.o main.o initfunc.o l2error.o

VIZ = pv3routines.o

LIB = -L$(HOME)/Codes/lib -lutil -lmyblas -lmesh -lquad -lbasis -lm
DEFINES = -Df2cFortran
CPPFLAGS = -I$(HOME)/Codes/include/

### Begin machine dependent

## Darwin:
## Machine specific for OS X
#CPPFLAGS += -O3 -Wno-long-double 
#LIB += -framework veclib
#DEFINES += -DNOUNDERSCORE -DPOINTERLONG -DNOUNS
#PV3HOME = ${HOME}/Packages/pV3/
#FORTRTL = /opt/ibmcmp/lib
#FORTLIB = -lxlf90_r.A -lxlfmath.A
#FORTMAIN =

## Linux:
##  Machine specific for Linux
CXX=icpc #(cares)
#CXX=g++3 #(shelob)
CPPFLAGS += -O3
CPPFLAGS += -cxxlib-icc #(cares)
#LIB += -llapack -lblas -lg2c #(shelob)
LIB += -L$(HOME) -L/usr/lib/gcc-lib/i386-pc-linux/3.2.3/ -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
#PV3HOME = /usr/local/pV3/
#DEFINES += -DUNDERSCORE -DPOINTERLONG
#FORTRTL = /usr/local/pgi/linux86/lib#(shelob)
#FORTMAIN = $(FORTRTL)/pgfmain.o #(shelob)
#FORTLIB = -lpgftnrtl -lpgc #(shelob)


## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)
#LIB += -lcomplib.sgimath

## Aix:
## Machine specific for Jupiter
#CPPFLAGS += -O3
#LIB += -lessl \
/afs/cu/software/Development/lapack-3.0/rs_aix52/LAPACK/lapack_RS6K.a -lxlf -lxlf90
#DEFINES = -DIBMR2Fortran -DNOUNDERSCORE -DPOINTERLONG -DNOUNS

### End machine dependent
CPPFLAGS += $(DEFINES)
PVLIB = -L$(PV3HOME)/clients/$(PVM_ARCH) -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3 -L$(FORTRTL) $(FORTLIB)

hp: $(OBJ)
	${CXX} -o $@ $(CPPFLAGS) $(OBJ) $(LIB)
	
hp.pV3: $(OBJ) $(VIZ)
	${CXX} -o $@ $(OBJ) $(FORTMAIN) $(LIB) $(VIZ) $(PVLIB)

$(OBJ): spectral_hp.h hp_mgrid.h blocks.h defines.h surface.h
