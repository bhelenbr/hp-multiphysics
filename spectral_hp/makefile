SPHP = allocate.o copy.o tobasis.o gtol.o output.o curvinit.o ptprobe.o
HPMG = mg_alloc.o rsdl.o minvrt.o movco.o movfin.o tstep.o nstage.o bdry.o tadvance.o length.o adapt.o surface.o drag.o
BLCKS 	= block.o blocks.o rblocks.o
PLANAR = main.o initfunc.o
VIZ = pv3routines.o

#CXX = g++3


# LINUX 
MYLIB = -L$(HOME)/Codes/lib -lutil -lmyblas -lmesh -lquad -lbasis 
PVLIB = \
-L/usr/local/pV3/clients/$(PVM_ARCH) -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
#LIB = -lm -llapack -lf77blas -lcblas -latlas -lg2c
LIB = -lm -llapack -lblas
DEFINES = -DPV3 -DUNDERSCORE -DPOINTERLONG 
CPPFLAGS = -O3 -Df2cFortran $(DEFINES) -I$(HOME)/Codes/include/

# OS X
#CXX = c++
#MYLIB = -L$(HOME)/Codes/lib -lutil -lmyblas -lmesh -lquad 
#PVLIB = \
#-L$(HOME)/Codes/packages/pV3 -L$(PVM_ROOT)/lib/$(PVM_ARCH) -lpV3 -lgpvm3 -lpvm3
#LIB = -L/sw/lib -lm -llapack -lcblas -lf77blas -latlas -lg2c
#DEFINES = -DPV3 -DNOUNDERSCORE -DPOINTERLONG -DNOUNS
#CPPFLAGS = -O1 $(DEFINES) -I$(HOME)/Codes/include/

# UNCOMMENT FOR SGI 
#LIB = $(MYLIB) -lm -lcomplib.sgimath 
#CPPFLAGS = -Ofast=IP30 -64 -mips4 $(OPT)
#CPPFLAGS = -Ofast=IP27
#CPPFLAGS = -Ofast=IPrk5

hp: $(PLANAR) $(VIZ) $(BLCKS) $(SPHP) hp_mgrid.h surface.h $(HPMG) blocks.h $(BLCKS)
# OSX
#	${CXX} ${CPPFLAGS} -o $@ $(BLCKS) $(HPMG) $(SPHP) $(HPBASE) $(PLANAR) $(VIZ) $(MYLIB) $(PVLIB) $(LIB) 
# LINUX
	/usr/local/pgi/linux86/bin/pgf77 ${CPPFLAGS} -o $@ $(BLCKS) $(HPMG) $(SPHP) $(HPBASE) $(PLANAR) $(VIZ) $(MYLIB) $(PVLIB) $(LIB) 
#
	rm *.o
