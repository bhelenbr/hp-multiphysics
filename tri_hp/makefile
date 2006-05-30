BASE = adapt.o allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o minvrt.o movco.o movfin.o \
nstage.o output.o ptprobe.o tadvance.o tobasis.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o ins/surface.o
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/length.o cd/rsdl.o cd/tstep.o
OBJPS = planestrain/allocate.o planestrain/bdry.o planestrain/getnewbdry.o planestrain/length.o planestrain/rsdl.o planestrain/tstep.o
OBJSWIRL = swirl/allocate.o swirl/bdry.o swirl/getnewbdry.o swirl/getnewibc.o swirl/length.o swirl/rsdl.o swirl/tstep.o

OBJ = $(BASE)
OBJ += $(OBJCD) $(OBJINS) $(OBJPS) $(OBJSWIRL)

VIZ = pv3routines.o

LIB = -L$(HOME)/Codes/lib -lmesh -lquad -lbasis -lmyblas -lutil -linput_map -lm
DEFINES = -Df2cFortran
CPPFLAGS = -I$(HOME)/Codes/include/ -I. -Iins -Icd -Iplanestrain -Iswirl 

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
CPPFLAGS += -O3 -cxxlib-icc -DNO_MPISRC -I/opt/mpich/myrinet/intel/include #(cares)
#CPPFLAGS += -cxxlib-icc #(cares)
#LIB += -llapack -lblas -lg2c #(shelob)
LIB += -L$(HOME) -L/usr/lib/gcc-lib/i386-pc-linux/3.2.3/ -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
LIB +=   -L/opt/mpich/myrinet/intel/lib -lpmpich++ -lmpich -L/opt/gm/lib -lgm -lpthread #(cares)
LDFLAGS += -static-libcxa #(cares)
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
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) $(LIB)
	
hp.pV3: $(OBJ) $(VIZ)
	${CXX} -o $@ $(OBJ) $(FORTMAIN) $(LIB) $(VIZ) $(PVLIB)

$(OBJ): *.h

