OBJ = initialize.o intgrt.o intgrt1d.o proj.o \
proj1d.o ptvalues.o hpprobe.o

CPPFLAGS = -I$(HOME)/Codes/include/ 
DEFINES = -Df2cFortran
LIB = -L$(HOME)/Codes/lib -lmyblas -lm

### Begin machine dependent

## Darwin:
## Machine specific for OS X
#AR = libtool
#ARFLAGS = -o
#CPPFLAGS += -O3 -funroll-loops -falign-loops=16
#LIB += -framework veclib
 
## Linux:
##  Machine specific for Linux
CXX=icpc #(cares)
#CXX=g++3 #(shelob)
#CPPFLAGS += -O3
CPPFLAGS += -O3 -cxxlib-icc
#LIB += -llapack -lblas -lg2c #(shelob)
LIB += -L. -latlas -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
#DEFINES = -DpgiFortran

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
#DEFINES = -DIBMR2Fortran

### End machine dependent

CPPFLAGS += $(DEFINES)

all: tri_basis libbasis.a 

tri_basis: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIB)

libbasis.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libbasis.a $(HOME)/Codes/lib/
	cp -f tet_basis.h $(HOME)/Codes/include/
	rm *.o

$(OBJ): tet_basis.h
