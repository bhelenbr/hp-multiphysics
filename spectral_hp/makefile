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
CPPFLAGS += -O3
LIB += -framework veclib
DEFINES += -DNOUNDERSCORE -DPOINTERLONG -DNOUNS

## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
#CXX=g++3 #(shelob)
#CPPFLAGS += -O3
#LIB += -llapack -lblas -lg2c
#PV3HOME = /usr/local/pV3/
#DEFINES += -DUNDERSCORE -DPOINTERLONG
#LD = /usr/local/pgi/linux86/bin/pgf77 #(shelob)

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
PGIHOME = /usr/local/pgi/linux86/lib
PVLIB = -L$(PV3HOME)/clients/$(PVM_ARCH) -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3 -L$(PGIHOME) -lpgftnrtl -lpgc

hp: $(OBJ)
	${CXX} -o $@ $(OBJ) $(LIB)
	
hp.pV3: $(OBJ) $(VIZ)
	${CXX} -o $@ $(OBJ) $(VIZ) $(PGIHOME)/pgfmain.o $(LIB) $(PVLIB)
	
$(OBJ): spectral_hp.h hp_mgrid.h blocks.h defines.h surface.h
