#!/bin/bash

SYS=$(uname)

echo '
OBJ = allocate.o copy.o tobasis.o gtol.o output.o curvinit.o ptprobe.o mg_alloc.o \
rsdl.o minvrt.o movco.o movfin.o tstep.o nstage.o bdry.o tadvance.o length.o adapt.o \
surface.o drag.o block.o blocks.o rblocks.o main.o initfunc.o
VIZ = pv3routines.o
MYLIB = -L$(HOME)/Codes/lib -lutil -lmyblas -lmesh -lquad -lbasis
' > makefile

if [ $SYS = "linux-gnu" ]; then

echo '
PVLIB = -L/usr/local/pV3/clients/$(PVM_ARCH) -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
LIB = -lm -llapack -lblas
DEFINES = -DPV3 -DUNDERSCORE -DPOINTERLONG 
CPPFLAGS = -O3 -Df2cFortran $(DEFINES) -I$(HOME)/Codes/include/

hp: $(OBJ) $(VIZ)
	/usr/local/pgi/linux86/bin/pgf77 ${CPPFLAGS} -o $@ $(OBJ) $(VIZ) $(MYLIB) $(PVLIB) $(LIB) 
' >> makefile

elif [ $SYS = "Darwin" ]; then

echo '
CXX = c++
PVLIB = -L$(HOME)/Packages/pV3 -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
LIB = -lm -framework veclib
DEFINES = -DNO_PV3 -DNOUNDERSCORE -DPOINTERLONG -DNOUNS
CPPFLAGS = -O3 -Df2cFortran $(DEFINES) -I$(HOME)/Codes/include/

hp: $(OBJ)
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) $(MYLIB) $(LIB)
' >> makefile

elif [ $SYS = "Irix" ]; then

echo '
LIB = -lm -lcomplib.sgimath 
CPPFLAGS = -Ofast=IP30
#CPPFLAGS = -Ofast=IP27
#CPPFLAGS = -Ofast=IPrk5

hp: $(OBJ) $(VIZ)
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) $(VIZ) $(MYLIB) $(PVLIB) $(LIB)
' >> makefile

elif [ $SYS = "Aix" ]; then

echo '
PVLIB = -L$(HOME)/Packages/pV3 -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
LIB = -lm -lessl -lxlf -lxlf90
DEFINES = -DNOUNDERSCORE -DPOINTERLONG -DNOUNS
CCFLAGS = -O3 -DIBMR2Fortran $(DEFINES) -I$(HOME)/Codes/include/
hp: $(OBJ)
	${CCC} -O3 -o $@ $(OBJ) $(MYLIB) $(PVLIB) /afs/cu/software/Development/lapack-3.0/rs_aix52/LAPACK/lapack_RS6K.a $(LIB) 
' >> makefile
fi