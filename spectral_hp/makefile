OBJ = allocate.o copy.o tobasis.o gtol.o output.o curvinit.o ptprobe.o mg_alloc.o \
rsdl.o minvrt.o movco.o movfin.o tstep.o nstage.o bdry.o tadvance.o length.o adapt.o\
 surface.o drag.o block.o blocks.o rblocks.o main.o initfunc.o
VIZ = pv3routines.o

MYLIB = -L$(HOME)/Codes/lib -lutil -lmyblas -lmesh -lquad -lbasis 

# LINUX 
PVLIB = \
-L/usr/local/pV3/clients/$(PVM_ARCH) -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
LIB = -lm -llapack -lblas
DEFINES = -DPV3 -DUNDERSCORE -DPOINTERLONG 
CPPFLAGS = -O3 -Df2cFortran $(DEFINES) -I$(HOME)/Codes/include/

# OS X
#CXX = c++
#PVLIB = -L$(HOME)/Packages/pV3 -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
#LIB = -lm -framework veclib
#DEFINES = -DNO_PV3 -DNOUNDERSCORE -DPOINTERLONG -DNOUNS
#CPPFLAGS = -O1 -Df2cFortran $(DEFINES) -I$(HOME)/Codes/include/

# UNCOMMENT FOR SGI 
#LIB = -lm -lcomplib.sgimath 
#CPPFLAGS = -Ofast=IP30 -64 -mips4 $(OPT)
#CPPFLAGS = -Ofast=IP27
#CPPFLAGS = -Ofast=IPrk5

hp: $(OBJ) $(VIZ)
	/usr/local/pgi/linux86/bin/pgf77 ${CPPFLAGS} -o $@ $(OBJ) $(VIZ) $(MYLIB) $(PVLIB) $(LIB) 

hp.osx: $(OBJ)
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) $(MYLIB) $(PVLIB) $(LIB) 

hp.sgi: $(OBJ) $(VIZ)
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) $(VIZ) $(MYLIB) $(PVLIB) $(LIB) 
